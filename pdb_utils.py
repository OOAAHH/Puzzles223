import copy
import math
import os,sys
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Any
from pathlib import Path

from Bio.PDB import *
import pandas as pd

import msgs,math,numpy
import mcannotate
import rnaview
import utils

@dataclass
class PDBConfig:
	"""PDB工具配置类"""
	BIN_DIR: str = os.getcwd()
	BACKBONE_ATOMS: List[str] = ["C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'","O4'", "O5'", "OP1", "OP2", "P"]
	HEAVY_ATOMS: List[str] = ["C2", "C4", "C5", "C6", "C8", "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"]
	ALL_ATOMS: List[str] = BACKBONE_ATOMS + HEAVY_ATOMS
	RMSDD_ATOMS: List[str] = ["C4", "C8", "P", "C1'"]
	MAX_ERRORS: int = 5

	@property
	def all_atoms(self) -> List[str]:
		"""获取所有原子列表"""
		return self.BACKBONE_ATOMS + self.HEAVY_ATOMS

# 创建全局配置对象
config = PDBConfig()

# from: http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
# Implements the Gauss error function.
#   erf(z) = 2 / sqrt(pi) * integral(exp(-t*t), t = 0..z)
#
# fractional error in math formula less than 1.2 * 10 ^ -7.
# although subject to catastrophic cancellation when z in very close to 0
# from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
def erf(z):
		t = 1.0 / (1.0 + 0.5 * abs(z))
		# use Horner's method
		ans = 1 - t * math.exp( -z*z -  1.26551223 +
								  t * ( 1.00002368 +
								  t * ( 0.37409196 + 
								  t * ( 0.09678418 + 
								  t * (-0.18628806 + 
								  t * ( 0.27886807 + 
								  t * (-1.13520398 + 
								  t * ( 1.48851587 + 
								  t * (-0.82215223 + 
								  t * ( 0.17087277))))))))))
		if z >= 0.0:
				return ans
		else:
				return -ans

class PDBNormalizer:
	"""PDB文件标准化处理类"""
	
	def __init__(self, fres_list: str, fatoms_list: str):
		"""初始化标准化器
		
		Args:
			fres_list: 残基列表文件路径
			fatoms_list: 原子列表文件路径
		"""
		self._res_list: Dict[str, str] = {}
		self._atom_list: Dict[str, str] = {}
		self._in_model: bool = False
		self._in_atom: bool = False
		self._chain_found: bool = False
		self._row_count: int = 0
		self._ok: bool = True
		
		self._load_res_list(fres_list)
		self._load_atom_list(fatoms_list)
	
	def parse(self, finput: str, foutput: str) -> bool:
		"""解析PDB文件
		
		Args:
			finput: 输入文件路径
			foutput: 输出文件路径
			
		Returns:
			是否解析成功
		"""
		out_txt = []
		
		try:
			with open(finput) as fi:
				for row in fi:
					self._row_count += 1
					row = row.strip()
					rec_name = row[:6]
					
					if rec_name == "MODEL ":
						row = self.parse_model(row)
						row = ""
					elif rec_name == "ENDMDL":
						row = self.parse_endmdl(row)
						row = ""
					elif rec_name[:3] == "TER":
						row = self.parse_ter(row)
					elif rec_name in ("ATOM  ", "HETATM"):
						row = self.parse_atom(row)
					else:
						continue
					
					if row:
						out_txt.append(row)
			
			if self._in_atom:
				out_txt.append("TER")
			
			if self._ok:
				with open(foutput, "w") as f:
					f.write("\n".join(out_txt))
			
			return self._ok
			
		except Exception as e:
			msgs.show("ERROR", f"Error parsing PDB file: {str(e)}")
			return False
	
	def parse_model(self, row: str) -> str:
		"""解析MODEL记录"""
		if self._in_model:
			self.show_err("'ENDMDL' not found.")
		
		if self._in_atom:
			self.show_err("Missing 'MODEL' before 'ATOM' declaration.")
		
		self._in_model = True
		return row
	
	def parse_endmdl(self, row: str) -> str:
		"""解析ENDMDL记录"""
		if not self._in_model:
			msgs.show("Warning", "Missing 'MODEL' declaration.")
		
		if self._in_atom:
			msgs.show("Warning", "Missing 'TER' declaration.")
			
		self._in_model = False
		self._in_atom = False
		return "ENDMDL"
	
	def parse_ter(self, row: str) -> str:
		"""解析TER记录"""
		result = "TER" if self._in_atom else ""
		self._in_atom = False
		return result
	
	def parse_atom(self, row: str) -> str:
		"""解析ATOM记录"""
		# 获取字段
		serial = row[6:11]
		name = row[12:16].strip()
		altLoc = row[16]
		resName = row[17:20].strip()
		chainID = row[21]
		resSeq = int(row[22:26])
		iCode = row[26]
		x = row[30:38]
		y = row[38:46]
		z = row[46:54]
		occupancy = row[54:60]
		tempFactor = row[60:66]
		element = row[76:78]
		charge = row[78:80]
		
		if element.strip() == 'H':
			return ""
		
		# 检查残基名
		resName_norm = self._res_list.get(resName)
		if resName_norm is None:
			self.show_err(f"Unknown residue name: '{resName}'.")
			return ""
		elif resName_norm == "-":
			return ""
		
		resName = resName_norm
		
		# 检查原子名
		name_norm = self._atom_list.get(name)
		if name_norm is None:
			self.show_err(f"Unknown atom name: '{name}' in residue '{resName}'")
			return ""
		elif name_norm == "-":
			return ""
		
		name = name_norm.ljust(3)
		
		# 检查链ID
		if chainID == " ":
			if self._chain_found:
				self.show_err("One of the chains is missing!")
			else:
				chainID = "A"
		else:
			self._chain_found = True
		
		# 标准化字段
		occupancy = "  1.00"  # MolProbity要求
		tempFactor = "  0.00" if tempFactor == "" else tempFactor
		element = name[0] if element == "" else element
		
		self._in_atom = True
		return f"ATOM  {serial:5}  {name:3}{altLoc}{resName:3} {chainID}{resSeq:4d}{iCode}   {x:8}{y:8}{z:8}{occupancy:6}{tempFactor:6}          {element:2}{charge:2}"
	
	def show_err(self, msg: str) -> None:
		"""显示错误信息"""
		msgs.show("ERROR", f"Line {self._row_count}: {msg}\n")
		self._ok = False
	
	def _load_res_list(self, fres_list: str) -> None:
		"""加载残基列表"""
		try:
			with open(fres_list) as f:
				pairs = [line.split() for line in f if not line.startswith("#")]
				self._res_list = {name: nt for name, nt in pairs}
		except Exception as e:
			msgs.show("ERROR", f"Error loading residue list: {str(e)}")
	
	def _load_atom_list(self, fatoms_list: str) -> None:
		"""加载原子列表"""
		try:
			with open(fatoms_list) as f:
				pairs = [line.split() for line in f if not line.startswith("#")]
				self._atom_list = {name: name_norm for name, name_norm in pairs}
		except Exception as e:
			msgs.show("ERROR", f"Error loading atom list: {str(e)}")

#		
# get the sequence list from a pdb either raw or indexed
#
@dataclass
class Residue:
	"""残基类"""
	chain: str
	pos: int
	nt: str
	res: Any  # Bio.PDB.Residue类型
	
	def key(self) -> str:
		"""生成残基键值"""
		return f"{self.chain}:{self.pos}"
	
	def __str__(self) -> str:
		"""字符串表示"""
		return f"{self.chain}:{self.pos}:{self.nt} > {self.res}"

class PDBStruct:
	"""PDB结构类"""
	
	def __init__(self):
		"""初始化PDB结构"""
		self._pdb_file: Optional[str] = None
		self._struct: Optional[Any] = None  # Bio.PDB.Structure类型
		self._res_list: List[Residue] = []
		self._res_seq: List[int] = []
		self._res_index: Dict[str, List[Optional[int]]] = {}
		self._interactions: List[Tuple[str, int, int, str]] = []
	
	def load(self, pdb_file: str, index_name: Optional[str] = None) -> bool:
		"""加载PDB结构
		
		Args:
			pdb_file: PDB文件路径
			index_name: 索引文件路径
			
		Returns:
			是否加载成功
		"""
		self._pdb_file = pdb_file
		
		try:
			ok = self._load_struct()
			
			if ok and index_name is not None:
				ok = self._load_index(index_name)
			else:
				ok = self._load_index2()
				
			if ok:
				ok = self._load_annotations_3D()
				
			return ok
			
		except Exception as e:
			msgs.show("ERROR", f"Error loading PDB structure: {str(e)}")
			return False
	
	def raw_sequence(self) -> str:
		"""获取原始序列"""
		return "".join(self._res_list[ndx].nt for ndx in self._res_seq)
	
	def res_sequence(self) -> List[Any]:
		"""获取残基序列"""
		return [self._res_list[ndx].res for ndx in self._res_seq]
	
	def get_interactions(self, type: str = "ALL") -> List[Tuple[str, int, int, str]]:
		"""获取相互作用
		
		Args:
			type: 相互作用类型
			
		Returns:
			相互作用列表
		"""
		if type == "ALL":
			return self._interactions
		elif type in ("PAIR"):
			return [x for x in self._interactions if x[0] in ("PAIR_2D", "PAIR_3D")]
		elif type in ("PAIR_2D", "PAIR_3D", "STACK"):
			return [x for x in self._interactions if x[0] == type]
		else:
			msgs.show("FATAL", f"Wrong interaction type '{type}' expected: 'ALL', 'PAIR', 'PAIR_2D', 'PAIR_3D' or 'STACK'")
			return []
	
	@property
	def struct(self) -> Optional[Any]:
		"""获取结构对象"""
		return self._struct
	
	@property
	def res_seq(self) -> List[int]:
		"""获取残基序列索引"""
		return self._res_seq
	
	@property
	def res_list(self) -> List[Residue]:
		"""获取残基列表"""
		return self._res_list
	
	@property
	def pdb_file(self) -> Optional[str]:
		"""获取PDB文件路径"""
		return self._pdb_file
	
	def rad_gir(self) -> float:
		"""计算半径"""
		rmean = Vector([0.0, 0.0, 0.0])
		count = 0
		
		for res in self._res_list:
			for a in res.res:
				rmean += a.get_vector()
				count += 1
		
		rmean = rmean/float(count)
		
		rsum = 0.0
		for res in self._res_list:
			for a in res.res:
				rsum += (a.get_vector() - rmean) * (a.get_vector() - rmean)
		
		return math.sqrt(rsum / count)
	
	def _load_struct(self) -> bool:
		"""加载结构"""
		try:
			parser = PDBParser()
			self._struct = parser.get_structure("struct", self._pdb_file)
			
			if len(self._struct) > 1:
				msgs.show("WARNING", f"{len(self._struct)} models found. Only the first will be used!")
			
			self._res_list = []
			self._res_seq = []
			self._res_index = {}
			
			# 只获取第一个模型
			model = self._struct[0]
			count = 0
			for chain in model.child_list:
				for res in chain.child_list:
					new_residue = Residue(chain.id, res.id[1], res.resname.strip(), res)
					self._res_list.append(new_residue)
					self._res_seq.append(count)
					self._res_index[new_residue.key()] = [count, None]
					count += 1
			
			return True
			
		except Exception as e:
			msgs.show("ERROR", f"Error loading structure: {str(e)}")
			return False
	
	def _load_index(self, index_name: str) -> bool:
		"""加载索引"""
		try:
			self._res_seq = []
			entries = []
			
			with open(index_name) as f:
				for row in f:
					row = row.strip()
					if not row.startswith("#") and row:
						entries.extend(row.split(","))
			
			for entry in entries:
				parts = entry.split(":")
				if len(parts) != 3:
					msgs.show("ERROR", f"Bad index entry: '{entry}'")
					return False
				
				chain = parts[0]
				pos = int(parts[1])
				count = int(parts[2])
				
				# 获取索引位置
				ndx = self._get_index(chain, pos, 0)
				if ndx is None:
					return False
				
				# 获取位置
				for i in range(ndx, ndx + count):
					if i >= len(self._res_list):
						msgs.show("ERROR", f"Bad count {count} in index entry: '{entry}'")
						return False
					
					if self._res_list[i].chain != chain:
						msgs.show("ERROR", f"Position {i} in index entry: '{entry}' is outside the chain")
						return False
					
					self._res_seq.append(i)
					self._res_index[self._res_list[i].key()][1] = len(self._res_seq) - 1
			
			return True
			
		except Exception as e:
			msgs.show("ERROR", f"Error loading index: {str(e)}")
			return False
	
	def _load_index2(self) -> bool:
		"""加载默认索引"""
		self._res_seq = []
		for i in range(len(self._res_list)):
			self._res_seq.append(i)
			self._res_index[self._res_list[i].key()][1] = len(self._res_seq) - 1
		return True
	
	def _load_annotations_3D(self) -> bool:
		"""加载3D注释"""
		try:
			self._interactions = []
			mca = mcannotate.MCAnnotate()
			mca.load(self._pdb_file, os.path.dirname(self._pdb_file))
			
			for (type, chain_a, pos_a, nt_a, chain_b, pos_b, nt_b, extra1, extra2, extra3) in mca.interactions:
				# 获取配对的第一个位置的排名
				rank_a = self._get_index(chain_a, pos_a, 1)
				rank_b = self._get_index(chain_b, pos_b, 1)
				
				if rank_a is None or rank_b is None:
					continue
				
				extra = extra1 if type == "STACK" else f"{extra1}{extra2}"
				self._interactions.append((type, min(rank_a, rank_b), max(rank_a, rank_b), extra))
			
			return True
			
		except Exception as e:
			msgs.show("ERROR", f"Error loading 3D annotations: {str(e)}")
			return False
	
	def _get_index(self, chain: str, pos: int, field: int) -> Optional[int]:
		"""获取索引
		
		Args:
			chain: 链ID
			pos: 位置
			field: 字段索引
			
		Returns:
			索引值
		"""
		key = f"{chain}:{pos}"
		idx = self._res_index.get(key)
		
		if idx is None:
			sys.stderr.write(f"WARNING\tBad index key: {idx} '{key}' {field}\n")
			return None
		 
		data = idx[field]
		if data is None and field == 0:
			sys.stderr.write(f"WARNING\tBad index key: {idx} '{key}' {field}\n")
		
		return data

class PDBComparer:
	"""PDB结构比较类"""
	
	def __init__(self):
		"""初始化比较器"""
		pass
	
	def mcq(self, f1: str, f2: str) -> float:
		"""计算MCQ分数
		
		Args:
			f1: 第一个PDB文件
			f2: 第二个PDB文件
			
		Returns:
			MCQ分数
		"""
		cmd = f'/home/chichau/RNA-Puzzles/bin/mcq-cli-2019-03-18/mcq-cli/mcq-global {f1} {f2} >mcq.log'
		utils.command(cmd)
		
		try:
			ss = open('mcq.log').read().strip().split()[-1]
			df = pd.read_csv(f'{ss}/matrix.csv', index_col=0, sep=',')
			v = df.iloc[0,1]
			print(f"MCQ: {v}")
		except Exception:
			print("mcq error")
			v = 0
		return v
	
	def tm(self, f1: str, f2: str) -> float:
		"""计算TM-score
		
		Args:
			f1: 第一个PDB文件
			f2: 第二个PDB文件
			
		Returns:
			TM-score
		"""
		cmd = f'/home/chichau/RNA-Puzzles/bin/RNAalign/RNAalign {f1} {f2} >tm.log'
		utils.command(cmd)
		
		try:
			v = float(open('tm.log').read().split('\n')[14][9:17])
			print(f"TMscore: {v}")
		except Exception:
			print("TM-score error")
			v = 0
		return v
	
	def gdt(self, f1: str, f2: str) -> float:
		"""计算GDT分数
		
		Args:
			f1: 第一个PDB文件
			f2: 第二个PDB文件
			
		Returns:
			GDT分数
		"""
		cmd = f'java -jar {config.BIN_DIR}/gdt.jar {f2} {f1} >gdt.log'
		os.system(cmd)
		
		try:
			x = open('gdt.log').read().strip().split('\n')[1].split(',')[-1]
			if x == 'NaN':
				return 0
			v = float(x)
		except Exception:
			v = 0
		return v
	
	def rmsd(self, src_struct: PDBStruct, trg_struct: PDBStruct, fit_pdb: Optional[str] = None) -> Optional[float]:
		"""计算RMSD
		
		Args:
			src_struct: 源结构
			trg_struct: 目标结构
			fit_pdb: 拟合后的PDB文件路径
			
		Returns:
			RMSD值
		"""
		src_residues = src_struct.res_sequence()
		trg_residues = trg_struct.res_sequence()
		
		atoms = self._get_atoms_struct(config.ALL_ATOMS, src_residues, trg_residues)
		if atoms is None:
			return None
			
		src_atoms, trg_atoms = atoms
		
		# 计算RMSD值并应用到目标结构
		sup = Superimposer()
		print(f"src_atoms {len(src_atoms)} trg_atoms {len(trg_atoms)}")
		sup.set_atoms(src_atoms, trg_atoms)
		
		# 复制结构以保持目标结构不变
		fit_struct = copy.deepcopy(trg_struct.struct)
		sup.apply(fit_struct.get_atoms())
		
		# 保存拟合后的结构
		if fit_pdb is not None:
			io = PDBIO()
			io.set_structure(fit_struct)
			io.save(fit_pdb)
		
		return sup.rms
	
	def rmsd2(self, src_struct: PDBStruct, trg_struct: PDBStruct, fit_pdb: Optional[str] = None) -> Optional[List[float]]:
		"""计算RMSD并返回距离列表
		
		Args:
			src_struct: 源结构
			trg_struct: 目标结构
			fit_pdb: 拟合后的PDB文件路径
			
		Returns:
			距离列表
		"""
		src_residues = src_struct.res_sequence()
		trg_residues = trg_struct.res_sequence()
		
		atoms = self._get_atoms_struct(config.ALL_ATOMS, src_residues, trg_residues)
		if atoms is None:
			return None
			
		src_atoms, trg_atoms = atoms
		
		# 计算RMSD值并应用到目标结构
		sup = Superimposer()
		print(f"src_atoms {len(src_atoms)} trg_atoms {len(trg_atoms)}")
		sup.set_atoms(src_atoms, trg_atoms)
		
		# 复制结构以保持目标结构不变
		fit_struct = copy.deepcopy(trg_struct.struct)
		sup.apply(fit_struct.get_atoms())
		
		fit_atoms = [i for i in fit_struct.get_atoms()]
		
		return [i-j for i,j in zip(fit_atoms, src_atoms)]
	
	def pvalue(self, m: float, N: int, param: str) -> float:
		"""计算p值
		
		Args:
			m: RMSD值
			N: 残基数
			param: 参数类型
			
		Returns:
			p值
		"""
		if param == "+":
			a = 5.1
			b = 15.8
		elif param == "-":
			a = 6.4
			b = 12.7
		else:
			msgs.show("FATAL", f"Wrong p-value parameter '{param}'. Expected '+' or '-'")
			return 0.0
			
		RMSD = a * (N ** 0.41) - b
		Z = (m - RMSD) / 1.8
		return (1.0 + erf(Z / (2**0.5))) / 2.0
	
	def INF(self, src_struct: PDBStruct, trg_struct: PDBStruct, type: str) -> float:
		"""计算INF分数
		
		Args:
			src_struct: 源结构
			trg_struct: 目标结构
			type: 相互作用类型
			
		Returns:
			INF分数
		"""
		P = TP = FP = FN = 0
		
		for (stype, sb1, sb2, sextra) in src_struct.get_interactions(type):
			P += 1
			found = False
			for (ttype, tb1, tb2, textra) in trg_struct.get_interactions(type):
				if (stype == ttype and sb1 == tb1 and sb2 == tb2 and sextra == textra):
					found = True
					break
			
			if found:
				TP += 1
			else:
				FN += 1
		
		for (ttype, tb1, tb2, textra) in trg_struct.get_interactions(type):
			found = False
			for (stype, sb1, sb2, sextra) in src_struct.get_interactions(type):
				if (stype == ttype and sb1 == tb1 and sb2 == tb2 and sextra == textra):
					found = True
					break
			
			if not found:
				FP += 1
		
		if TP == 0 and (FP == 0 or FN == 0):
			return -1.0
			
		PPV = float(TP) / (float(TP) + float(FP))
		STY = float(TP) / (float(TP) + float(FN))
		return (PPV * STY) ** 0.5
	
	def DP(self, src_struct: PDBStruct, trg_struct: PDBStruct, template_txt: str, dname: str, dp_script: str) -> None:
		"""生成DP配置
		
		Args:
			src_struct: 源结构
			trg_struct: 目标结构
			template_txt: 模板文本
			dname: 输出目录
			dp_script: DP脚本路径
		"""
		# 准备配置文件
		txt = []
		txt.append("matrix=True")
		txt.append("quiet_err = True")
		txt.append(f"out_dir = '{dname}'")
		txt.append(f"ref_model = ('{src_struct.pdb_file}', 0)")
		txt.append(f"cmp_model = [('{trg_struct.pdb_file}', 0)]")
		
		aligns = self._build_dp_alignments(src_struct, trg_struct)
		aligns_txt = [f"('{align[0]}', {align[1]}, '{align[2]}', {align[3]}, {align[4]})" for align in aligns]
		txt.append(f"aligns = [{', '.join(aligns_txt)}]")
		txt.append(template_txt)
		
		fname_cfg = f"{trg_struct.pdb_file}.cfg"
		fname_log = f"{trg_struct.pdb_file}.log"
		
		with open(fname_cfg, "w") as f:
			f.write("\n".join(txt))
		
		# 运行DP生成器
		os.system(f"python {dp_script} -c {fname_cfg} > {fname_log}")
	
	def VARNA(self, src_struct: PDBStruct, trg_struct: PDBStruct, algorithm: str = "radiate") -> Dict[str, Any]:
		"""生成VARNA配置
		
		Args:
			src_struct: 源结构
			trg_struct: 目标结构
			algorithm: 算法类型
			
		Returns:
			VARNA配置字典
		"""
		edges = {"W": "wc", "S": "s", "H": "h"}
		
		data = {
			"sequenceDBN": src_struct.raw_sequence(),
			"structureDBN": "." * len(src_struct.raw_sequence()),
			"algorithm": algorithm
		}
		
		aux_bps = []
		for (stype, sb1, sb2, sextra) in src_struct.get_interactions("PAIR"):
			color = "#FF0000"
			for (ttype, tb1, tb2, textra) in trg_struct.get_interactions("PAIR"):
				if (stype == ttype and sb1 == tb1 and sb2 == tb2 and sextra == textra):
					color = "#00FF00"
					break
			aux_bps.append(f"({sb1+1},{sb2+1}):color={color},edge5={edges[sextra[0]]},edge3={edges[sextra[1]]},stericity={sextra[2:]}")
		
		data["auxBPs"] = ";".join(aux_bps)
		return data
	
	def _get_atoms_residue(self, atom_list: List[str], src_res: Any, trg_res: Any) -> Tuple[List[Any], List[Any]]:
		"""获取残基中的原子
		
		Args:
			atom_list: 原子列表
			src_res: 源残基
			trg_res: 目标残基
			
		Returns:
			源原子列表和目标原子列表
		"""
		src_atom_list = []
		trg_atom_list = []
		
		src_atom_list_tmp = [a for a in src_res if a.get_name() in atom_list]
		trg_atom_list_tmp = [a for a in trg_res if a.get_name() in atom_list]
		
		# 对每个参考原子
		for src_atom in src_atom_list_tmp:
			found = False
			src_name = src_atom.get_full_id()[4][0]
			
			# 在比较结构中搜索同名原子
			for trg_atom in trg_atom_list_tmp:
				trg_name = trg_atom.get_full_id()[4][0]
				
				# 如果找到原子则保存并继续下一个
				if src_name == trg_name:
					src_atom_list.append(src_atom)
					trg_atom_list.append(trg_atom)
					found = True
					break
			
			if not found:
				msgs.show("WARNING", f"Atom {src_name} from residue {src_res.id} not found in target atom list")
		
		return src_atom_list, trg_atom_list
	
	def _get_atoms_struct(self, atom_list: List[str], src_residues: List[Any], trg_residues: List[Any]) -> Optional[Tuple[List[Any], List[Any]]]:
		"""获取结构中的原子
		
		Args:
			atom_list: 原子列表
			src_residues: 源残基列表
			trg_residues: 目标残基列表
			
		Returns:
			源原子列表和目标原子列表
		"""
		if len(src_residues) != len(trg_residues):
			msgs.show("ERROR", "Different number of residues!")
			return None
		
		src_atoms = []
		trg_atoms = []
		
		for src_res, trg_res in zip(src_residues, trg_residues):
			sa, ta = self._get_atoms_residue(atom_list, src_res, trg_res)
			src_atoms.extend(sa)
			trg_atoms.extend(ta)
		
		return src_atoms, trg_atoms
	
	def _build_dp_alignments(self, src_struct: PDBStruct, trg_struct: PDBStruct) -> List[List[Any]]:
		"""构建DP对齐
		
		Args:
			src_struct: 源结构
			trg_struct: 目标结构
			
		Returns:
			对齐列表
		"""
		aligns = []
		
		schain = tchain = ""
		spos = tpos = -1
		count = 0
		item = None
		
		for i, j in zip(src_struct.res_seq, trg_struct.res_seq):
			sres, tres = src_struct.res_list[i], trg_struct.res_list[j]
			
			# 如果编号或链发生变化
			if (sres.pos != (spos+1) or tres.pos != (tpos+1) or 
				sres.chain != schain or tres.chain != tchain):
				if count > 0:
					item[4] = count
					aligns.append(item)
				
				schain, tchain = sres.chain, tres.chain
				item = [sres.chain, sres.pos, tres.chain, tres.pos, None]
				count = 1
			else:
				count += 1
			
			spos, tpos = sres.pos, tres.pos
		
		# 添加最后一个对齐
		if count > 0:
			item[4] = count
			aligns.append(item)
		
		return aligns
