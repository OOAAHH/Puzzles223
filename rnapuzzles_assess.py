#!/usr/bin/python

import sys,os
from operator import attrgetter
from Bio.PDB import *
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
import pandas as pd

import pdb_utils
import utils
import extract


class RNAPuzzleConfig:
	"""RNA Puzzle评估配置类"""
	def __init__(self, puzzle_name, work_dir, pvalue=None):
		self.puzzle_name = puzzle_name
		self.work_dir = work_dir
		self.pvalue = pvalue
		self.setup_paths()
		
	def setup_paths(self):
		"""设置所有路径配置"""
		# MolProbity配置
		self.molprobity_bin_dir = '/home/chichau/Software/MolProbity/cmdline'
		
		# 基础路径配置
		self.dp_script = f"{self.work_dir}/bin/DPS_1.0.2/dp.py"
		self.config_dir = f"{self.work_dir}/bin/config"
		self.residues_list = f"{self.config_dir}/residues.list"
		self.atoms_list = f"{self.config_dir}/atoms.list"
		
		# 主要目录配置
		self.puzzle_dir = f"{self.work_dir}/Puzzles/{self.puzzle_name}"
		self.solution_dir = f"{self.work_dir}/Puzzles/solution/{self.puzzle_name}"
		self.original_dir = f"{self.work_dir}/Puzzles/original/{self.puzzle_name}"
		self.index_dir = f"{self.puzzle_dir}/index"
		
		# 步骤目录配置
		self.step1_dir = f"{self.puzzle_dir}/step1"
		self.step2_dir = f"{self.puzzle_dir}/step2"
		self.step3_dir = f"{self.puzzle_dir}/step3"
		self.step4_dir = f"{self.puzzle_dir}/step4"
		self.step5_dir = f"{self.puzzle_dir}/step5"
		self.molprobity_dir = f"{self.puzzle_dir}/molprobity"
		
		# 输出文件配置
		self.evals_file = f"{self.puzzle_dir}/evals.txt"
		self.table_file = f"{self.work_dir}/Puzzles/table/{self.puzzle_name}.txt"
		
		# Markdown和数据目录配置
		self.markdown_file = f"/home/chichau/rnapuzzles/source/_posts/2000-01-01-{self.puzzle_name}-3d.markdown"
		self.data_dir = f"/home/chichau/rnapuzzles/source/data/{self.puzzle_name}"

class RNAPuzzleAssessor:
	"""RNA Puzzle评估类
	
	此类实现了RNA结构预测的完整评估流程，包括：
	1. 目录准备
	2. 解决方案标准化
	3. 预测结构标准化
	4. 结构测量与评估
	"""
	
	def __init__(self, config):
		"""初始化评估器
		
		Args:
			config: RNAPuzzleConfig实例，包含所有配置信息
		"""
		self.config = config
		self.pdb_normalizer = pdb_utils.PDBNormalizer(config.residues_list, config.atoms_list)
		self.comparer = pdb_utils.PDBComparer()
		
	def assess(self):
		"""执行完整的评估流程
		
		评估流程包括以下步骤：
		1. 创建工作目录 (Step 0)
		2. 标准化解决方案文件 (Step 1)
		3. 标准化预测文件 (Step 2)
		4. 测量和评估结构 (Step 3)
		"""
		self.make_directory()  # Step 0
		self.normalize_solutions()  # Step 1
		self.normalize_predictions()  # Step 2
		self.measure()  # Step 3
		
	def make_directory(self):
		"""创建评估所需的目录结构 (Step 0)
		
		创建以下目录:
		- puzzle_dir/: 主工作目录
		- solution_dir/: 解决方案文件目录
		- original_dir/: 原始预测文件目录
		- index_dir/: 索引文件目录
		- step1_dir/ 到 step5_dir/: 处理步骤目录
			* step1_dir/: 存储清理后的文件
			* step2_dir/: 存储分割后的文件
			* step3_dir/: 存储标准化后的文件
			* step4_dir/: 存储提取后的文件
			* step5_dir/: 存储最终评估结果
		- molprobity_dir/: MolProbity分析目录
		- data_dir/: 最终数据目录
			* data_dir/pdb/: 结构文件
			* data_dir/dp/: 结构对比图
		"""
		dirs = [self.config.puzzle_dir, self.config.solution_dir, self.config.original_dir, 
				self.config.index_dir, self.config.molprobity_dir, self.config.step1_dir, 
				self.config.step2_dir, self.config.step3_dir, self.config.step4_dir, 
				self.config.step5_dir, self.config.data_dir, 
				f'{self.config.data_dir}/pdb', f'{self.config.data_dir}/dp']
		for d in dirs:
			utils.MakeDIR(d)
			
	def normalize_solutions(self):
		"""标准化解决方案文件 (Step 1)
		
		处理流程:
		1. 读取split.txt文件获取解决方案信息
		2. 对每个解决方案文件:
			- 复制并重命名为标准格式
			- 清理文件格式
			- 标准化PDB文件
			- 如果需要，提取指定残基
			- 生成索引文件
		
		输入文件:
		- solution_dir/split.txt: 解决方案配置文件
		- solution_dir/*.pdb: 原始解决方案PDB文件
		
		输出文件:
		- solution_dir/{puzzle_name}_solution_{i}.pdb: 标准化的解决方案文件
		- solution_dir/{puzzle_name}_extract_{i}.pdb: 提取的结构文件
		- solution_dir/{puzzle_name}_solution_{i}.index: 残基索引文件
		
		错误处理:
		- 如果split.txt不存在，记录错误
		- 如果标准化失败，记录错误
		"""
		sys.stderr.write("Step1:	Normalize solution file!\n")
		os.chdir(self.config.solution_dir)
		SPLIT_FILE = 'split.txt'
		if os.path.isfile(SPLIT_FILE):
			splits = open(SPLIT_FILE).read().strip().split("\n")
			for (i, line) in enumerate(splits):
				xx = line.split('\t')
				SOLUTION_NAME = f'{xx[0]}.pdb'
				SOLUTION_NORMAL = f'{self.config.puzzle_name}_solution_{i}.pdb'
				os.system(f"cp {SOLUTION_NAME} {SOLUTION_NORMAL}")
				utils.CleanFormat(SOLUTION_NORMAL)
				ok = self.pdb_normalizer.parse(SOLUTION_NORMAL, SOLUTION_NORMAL)
				if not ok:
					sys.stderr.write("ERROR:	solution file not normalized!\n")
				if len(xx) > 1:
					coords = xx[1]
					extract.extract_PDB(SOLUTION_NORMAL, coords, f'{self.config.puzzle_name}_extract_{i}.pdb')
					open(f'{self.config.puzzle_name}_solution_{i}.index', 'w').write(coords)
					sys.stderr.write(f"INFO:	Solution {i} extracted\n")
				else:
					os.system(f"cp {SOLUTION_NORMAL} {self.config.puzzle_name}_extract_{i}.pdb")
		else:
			sys.stderr.write('ERROR	split.txt not found!')
			
	def normalize_predictions(self):
		"""标准化预测文件 (Step 2)
		
		处理流程:
		1. 读取预测文件列表
		2. 对每个预测文件:
			- 清理文件格式
			- 重命名文件
			- 如果是多模型文件，进行分割
			- 标准化PDB格式
			- 提取指定残基(如果有索引文件)
		
		输入文件:
		- original_dir/*.pdb: 原始预测PDB文件
		- work_dir/Puzzles/list/{puzzle_name}.list: 预测文件列表
		- index_dir/*.index: 残基索引文件
		
		输出文件:
		- step1_dir/*.pdb: 清理后的文件
		- step2_dir/*.pdb: 分割后的文件
		- step3_dir/*.pdb: 标准化后的文件
		- step4_dir/*.pdb: 提取后的文件
		
		错误处理:
		- 如果预测文件未注册，记录错误
		- 如果序列不匹配，记录错误
		- 如果标准化失败，记录错误
		
		Returns:
			list: 处理后的文件列表
		"""
		pdb_files = [f for f in os.listdir(self.config.original_dir) if f.endswith(".pdb")]
		result_list = pd.read_csv(f"{self.config.work_dir}/Puzzles/list/{self.config.puzzle_name}.list", 
								header=None, index_col=0, sep=' ')
		
		sol_struct = pdb_utils.PDBStruct()
		sol_index = f'{self.config.solution_dir}/{self.config.puzzle_name}_solution_0.index'
		if not os.path.isfile(sol_index):
			sol_index = None
		sol_struct.load(f'{self.config.solution_dir}/{self.config.puzzle_name}_solution_0.pdb', sol_index)
		sol_raw_seq = sol_struct.raw_sequence()
		
		outs = []
		for pdb_file in sorted(pdb_files):
			print(f"=============================================={pdb_file}")
			utils.CleanFormat(f'{self.config.original_dir}/{pdb_file}')
			sys.stderr.write(f"INFO	Original PDB: '{pdb_file}'.\n")
			
			if pdb_file not in result_list.index:
				sys.stderr.write(f"ERROR	Unregistered pdb file: {pdb_file}\n")
				continue
				
			name1 = result_list.loc[pdb_file][1]
			sys.stderr.write(f"INFO	Renamed PDB: '{name1}'.\n")
			os.system(f"cp {self.config.original_dir}/{pdb_file} {self.config.step1_dir}/{name1}.pdb")
			
			grp = name1.replace(self.config.puzzle_name,'').replace('.pdb','').split('_')[0]
			if name1.find('_') > 0:
				print('yes')
				os.system(f'cp {self.config.step1_dir}/{name1}.pdb {self.config.step2_dir}/{self.config.puzzle_name}_{name1}.pdb')
				outs.append(f'{self.config.puzzle_name}_{name1}.pdb')
			else:
				outs.extend(self._split_pdb(f'{self.config.step1_dir}/{name1}.pdb', grp, name1))
				
		return outs
		
	def _split_pdb(self, f, id, N):
		"""将PDB文件分割成多个文件"""
		p = PDBParser(PERMISSIVE=1)
		io = PDBIO()
		out = []
		struct = p.get_structure('SI', f)
		for i, mdl in enumerate(struct):
			ent = Structure(struct.id)
			ent.add(mdl)
			io.set_structure(ent)
			io.save(f'{self.config.step2_dir}/{self.config.puzzle_name}_{id}_{i+1}.pdb')
			out.append(f'{self.config.puzzle_name}_{id}_{i+1}.pdb')
		return out

	def measure(self):
		"""执行结构测量与评估 (Step 3)
		
		处理流程:
		1. 对每个预测结构:
			- 加载结构和索引
			- 与所有解决方案比较，选择最佳匹配
			- 计算评估指标:
				* RMSD
				* P-value
				* MCQ分数
				* TM-score
				* INF分数(ALL/WC/NWC/STACK)
			- 生成结构对比图
		2. 运行MolProbity分析
		3. 生成评估报告
		
		输入文件:
		- step4_dir/*.pdb: 预测结构文件
		- solution_dir/*_solution_*.pdb: 解决方案文件
		- work_dir/Puzzles/dp/{puzzle_name}.cfg: DP配置模板
		
		输出文件:
		- step5_dir/*.pdb: 叠合后的结构文件
		- step5_dir/*.svg: 结构对比图
		- step5_dir/*.png: 结构对比图片
		- molprobity_dir/molprobity.txt: MolProbity分析结果
		- puzzle_dir/evals.txt: 评估结果文件
		- table_dir/{puzzle_name}.txt: 排名表格
		- table_dir/{puzzle_name}.csv: CSV格式排名表格
		- markdown_file: 结果展示页面
		
		评估指标:
		- RMSD: 结构偏差
		- P-value: 统计显著性
		- MCQ: 模型质量
		- TM-score: 拓扑相似度
		- INF: 相互作用网络保真度
		- Clash score: 原子碰撞分数
		"""
		VARNA_ALGORITHM = 'radiate'
		DP_TEMPLATE = f"{self.config.work_dir}/Puzzles/dp/{self.config.puzzle_name}.cfg"
		dp_template_txt = open(f"{self.config.work_dir}/Puzzles/dp/empty.cfg").read()
		if os.path.isfile(DP_TEMPLATE):
			dp_template_txt = open(DP_TEMPLATE).read()
			
		pdb_files = [f for f in os.listdir(self.config.step4_dir) if f.endswith(".pdb")]
		solutions = sorted([f"{self.config.solution_dir}/{x}" for x in os.listdir(self.config.solution_dir) 
						  if f"{self.config.puzzle_name}_solution_" in x and x.endswith(".pdb")])
		eval_list = []
		
		for pdb_file in sorted(pdb_files):
			print(pdb_file)
			res_struct = pdb_utils.PDBStruct()
			index_file = f"{self.config.index_dir}/{pdb_file.replace('PZ','').replace('.pdb','')}.index"
			
			if not os.path.isfile(index_file):
				if pdb_file.find('_') > 0:
					index_file = f"{self.config.index_dir}/{pdb_file.split('_')[1]}.index"
					if not os.path.isfile(index_file):
						index_file = f"{self.config.index_dir}/xx.index"
						if not os.path.isfile(index_file):
							index_file = None
							
			res_struct.load(f"{self.config.step3_dir}/{pdb_file}", index_file)
			res_raw_seq = res_struct.raw_sequence()
			
			# 与所有可用解决方案比较并选择最佳方案
			best_rmsd = 1e100
			best_pvalue = 1e100
			best_sol_struct = None
			best_sol_ndx = -1
			best_sol_struct_file = None
			
			for i, solution_file in enumerate(solutions):
				sys.stderr.write(f"INFO	Get solution structure data from '{solution_file}'.\n")
				solution_index = solution_file.replace('.pdb','.index')
				if not os.path.isfile(solution_index):
					solution_index = None
					
				sol_struct = pdb_utils.PDBStruct()
				sol_struct.load(solution_file, solution_index)
				sol_raw_seq = sol_struct.raw_sequence()
				
				if sol_raw_seq != res_raw_seq:
					sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
					sys.stderr.write(f'{pdb_file}\n')
					sys.stderr.write(f"DATA Solution sequence --> '{sol_raw_seq}'\n")
					sys.stderr.write(f"DATA Result sequence   --> '{res_raw_seq}'\n")
					
					if len(sol_raw_seq) != len(res_raw_seq):
						sys.stderr.write(f"ERROR	Length different! {len(sol_raw_seq)} - {len(res_raw_seq)}\n")
						continue
						
				rmsd = self.comparer.rmsd(sol_struct, res_struct)
				print(f"INFO: {res_struct.pdb_file} {sol_struct.pdb_file} {index_file} {solution_index} RMSD: {rmsd}")
				sys.stderr.write(f"INFO Partial RMSD --> {rmsd}\n")
				pvalue = self.comparer.pvalue(rmsd, len(sol_raw_seq), self.config.pvalue)
				sys.stderr.write(f"INFO Partial P-Value --> {pvalue:e}\n")
				
				if rmsd < best_rmsd:
					best_rmsd = rmsd
					best_sol_struct = sol_struct
					best_sol_struct_file = solution_file
					best_sol_ndx = i
					
				if pvalue < best_pvalue:
					best_pvalue = pvalue
					
			if best_sol_struct is None:
				sys.stderr.write("ERROR	No solutions selected!\n")
				continue
				
			# 评估结果
			eval = utils.Eval(problem=self.config.puzzle_name, 
							lab=pdb_file.split('_')[1], 
							result=int(pdb_file.split('_')[-1].replace('.pdb','')))
			eval_list.append(eval)
			eval.best_sol_ndx = best_sol_ndx
			eval.result_fit = f"{self.config.step5_dir}/{pdb_file}"
			eval.rmsd = self.comparer.rmsd(best_sol_struct, res_struct, eval.result_fit)
			eval.pvalue = best_pvalue
			eval.mcq = self.comparer.mcq(f"{self.config.step4_dir}/{pdb_file}", best_sol_struct_file)
			print(f"mcq, {eval.mcq}")
			eval.tm = self.comparer.tm(f"{self.config.step4_dir}/{pdb_file}", best_sol_struct_file)
			print(f"TM-score, {eval.tm}")
			eval.gdt = 0
			
			# 计算INF（所有类型）
			INF = self.comparer.INF(best_sol_struct, res_struct, type="ALL")
			eval.INF_ALL = INF
			eval.DI_ALL = abs(eval.rmsd / INF)
			INF = self.comparer.INF(best_sol_struct, res_struct, type="PAIR_2D")
			eval.INF_WC = INF
			INF = self.comparer.INF(best_sol_struct, res_struct, type="PAIR_3D")
			eval.INF_NWC = INF
			INF = self.comparer.INF(best_sol_struct, res_struct, type="STACK")
			eval.INF_STACK = INF
			
			data = self.comparer.VARNA(best_sol_struct, res_struct, VARNA_ALGORITHM)
			
			svg_file = f"{self.config.step5_dir}/{eval.file_name()}.svg"
			png_file = f"{self.config.step5_dir}/{eval.file_name()}.png"
			
			if not os.path.isfile(svg_file):
				self.comparer.DP(best_sol_struct, res_struct, dp_template_txt, 
							   self.config.step5_dir, self.config.dp_script)
			else:
				sys.stderr.write("INFO	Skipping DP!\n")
				
			if not os.path.isfile(png_file):
				print("yes,no convert")
			else:
				sys.stderr.write("INFO Skipping convert!\n")
				
			eval.ok = True
			sys.stderr.write("OK Result was processed OK\n")
			
		# 使用Molprobity计算所有文件的clash score
		MOLPROBITY_FILE = f"{self.config.molprobity_dir}/molprobity.txt"
		if not os.path.isfile(MOLPROBITY_FILE):
			os.system(f"{self.config.molprobity_bin_dir}/reduce-build {self.config.solution_dir} {self.config.molprobity_dir}")
			os.system(f"{self.config.molprobity_bin_dir}/reduce-build {self.config.step4_dir} {self.config.molprobity_dir}")
			os.system(f"{self.config.molprobity_bin_dir}/oneline-analysis -nocbeta -norota -norama {self.config.molprobity_dir} > {MOLPROBITY_FILE}")
		else:
			sys.stderr.write(f"INFO Skipping MOLPROBITY file {MOLPROBITY_FILE}!\n")
			
		# 解析molprobity文件
		utils.molprobity_parse(MOLPROBITY_FILE, eval_list)
		# 保存完整结果到最终的"eval"文件
		utils.save_evals_list(eval_list, self.config.evals_file)
		
		evals = utils.load_evals_list(self.config.evals_file)
		utils.compute_evals_ranks(evals)
		print(evals)
		
		out = ''
		evals = sorted(evals, key=lambda x: x.rmsd_rank)
		for eval in evals:
			out += f'{eval.rankline()}\n'
		open(self.config.table_file, 'w').write(out)
		
		df = pd.read_csv(self.config.table_file, header=None, sep='\t', index_col=0)
		if df.shape[1] == 23:
			df.columns = ['model','RMSD','RMSD_rank','DI','DI_rank','INF_all','INF_all_rank',
						 'INF_wc','INF_wc_rank','INF_nwc','INF_nwc_rank','INF_stack',
						 'INF_stack_rank','clash','clash_rank','MCQ','MCQ_rank','TM',
						 'TM_rank','GDT','GDT_rank','p-value','solution']
		df.to_csv(self.config.table_file.replace('.txt','.csv'))
		
		utils.save_md(evals, self.config.markdown_file)
		
		os.system(f"cp {self.config.solution_dir}/*solution*.pdb {self.config.data_dir}/pdb/")
		os.system(f"cp {self.config.step5_dir}/*.pdb {self.config.data_dir}/pdb/")
		os.system(f"cp {self.config.step5_dir}/*.png {self.config.data_dir}/dp/")

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print("Usage: python rnapuzzles_assess.py <puzzle_name> <work_dir> [pvalue]")
		sys.exit(1)
		
	# 创建配置对象
	config = RNAPuzzleConfig(sys.argv[1], sys.argv[2])
	if len(sys.argv) > 3:
		config.pvalue = sys.argv[3]
		
	# 创建评估器并执行评估
	assessor = RNAPuzzleAssessor(config)
	assessor.assess()



