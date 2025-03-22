#!/usr/bin/python

from operator import attrgetter
import os
import msgs
import pandas as pd
from dataclasses import dataclass
from typing import List, Optional, Dict, Any
from pathlib import Path

@dataclass
class UtilsConfig:
    """工具类配置"""
    def __init__(self):
        self.rank_attributes = {
            "rmsd": False,
            "pvalue": False,
            "DI_ALL": False,
            "clashscore": False,
            "mcq": False,
            "gdt": False,
            "INF_ALL": True,
            "INF_WC": True,
            "INF_NWC": True,
            "INF_STACK": True,
            "tm": True
        }

class FileHandler:
    """文件处理类"""
    @staticmethod
    def make_dir(directory: str) -> None:
        """创建目录，如果已存在则提示"""
        if not os.path.isdir(directory):
            os.makedirs(directory)
        else:
            print(f'Directory "{directory}" already exists')

    @staticmethod
    def clean_format(file_path: str) -> None:
        """清理文件格式"""
        os.system(f"mac2unix -q {file_path}")
        os.system(f"dos2unix -q {file_path}")

    @staticmethod
    def execute_command(cmd: str) -> None:
        """执行命令并检查返回码"""
        ret_code = os.system(cmd)
        if ret_code != 0:
            print(f"FATAL: Command '{cmd}' ended with error code '{ret_code}'")

    @staticmethod
    def read_file(file_path: str) -> str:
        """读取文件内容"""
        with open(file_path, 'r') as f:
            return f.read().strip()

    @staticmethod
    def write_file(file_path: str, content: str) -> None:
        """写入文件内容"""
        with open(file_path, 'w') as f:
            f.write(content)

def MakeDIR(DIR):
    """创建目录（向后兼容）"""
    FileHandler.make_dir(DIR)

def CleanFormat(f):
    """清理文件格式（向后兼容）"""
    FileHandler.clean_format(f)

def command(cmd):
    """执行命令（向后兼容）"""
    FileHandler.execute_command(cmd)

class Result:
    """结果类（向后兼容）"""
    def __init__(self, problem, original_file, lab, result):
        self.problem = problem
        self.file_original = original_file
        self.lab = lab
        self.result = result
        self.file = f"{problem}_{lab}_{result}.pdb"

def read_results_list(problem, fname):
    """读取结果列表（向后兼容）"""
    return FileHandler.read_file(fname).strip().split("\n")

def get_index_file(pdb_file, pdb_result_file=""):
    """获取索引文件（向后兼容）"""
    index_file = f"{pdb_file.replace('.pdb', '')}.{pdb_result_file.replace('.pdb', '')}.index"
    
    if not pdb_result_file or not os.path.isfile(index_file):
        index_file = pdb_file.replace('.pdb', '.index')
    
    if not os.path.isfile(index_file):
        msgs.show("INFO", f"INDEX SKIPPED! '{index_file}' skipped for '{pdb_file}'.")
        return None
    else:
        msgs.show("INFO", f"INDEX FOUND! '{index_file}' for '{pdb_file}'.")
        return index_file

def molprobity_parse(fname: str, evals: List[Eval]) -> None:
    """解析MolProbity结果"""
    try:
        for line in FileHandler.read_file(fname).split('\n'):
            if line and "#" not in line:
                data = line.split(":")
                if len(data) < 9:
                    continue
                    
                id_parts = data[0].strip("FH.pdb").split("_")
                if len(id_parts) == 3 and id_parts[1] != 'solution' and not id_parts[2].startswith('model'):
                    eval = find_eval(evals, id_parts[0], id_parts[1], int(id_parts[2]))
                    if eval is None:
                        print(f'Cannot find {data[0]}')
                    else:
                        eval.clashscore = float(data[8])
                else:
                    print(f"skip: {line}")
    except Exception as e:
        print(f"FATAL: Error parsing MolProbity file '{fname}': {str(e)}")

@dataclass
class Eval:
    """评估结果类"""
    problem: str = ""
    original: str = ""
    lab: str = ""
    result: str = ""
    result_fit: str = ""
    
    # 评估指标
    rmsd: float = 1e100
    pvalue: float = 1e100
    DI_ALL: float = 1e100
    INF_ALL: float = 0.0
    INF_WC: float = 0.0
    INF_NWC: float = 0.0
    INF_STACK: float = 0.0
    clashscore: float = 1e100
    best_sol_ndx: int = -1
    mcq: float = 1e100
    tm: float = 1e100
    gdt: float = 1e100
    
    # 排名
    rmsd_rank: int = 0
    pvalue_rank: int = 0
    DI_ALL_rank: int = 0
    INF_ALL_rank: int = 0
    INF_WC_rank: int = 0
    INF_NWC_rank: int = 0
    INF_STACK_rank: int = 0
    clashscore_rank: int = 0
    mcq_rank: int = 0
    tm_rank: int = 0
    gdt_rank: int = 0
    
    ok: bool = False
    
    def parse(self, row: str) -> bool:
        """解析评估结果行"""
        data = row.split()
        if len(data) != 15:
            return False
            
        try:
            self.problem = data[0]
            self.lab = data[1]
            self.result = int(data[2])
            self.rmsd = float(data[3])
            self.pvalue = float(data[4])
            self.DI_ALL = float(data[5])
            self.INF_ALL = float(data[6])
            self.INF_WC = float(data[7])
            self.INF_NWC = float(data[8])
            self.INF_STACK = float(data[9])
            self.clashscore = float(data[10])
            self.mcq = float(data[11])
            self.tm = float(data[12])
            self.gdt = float(data[13])
            self.best_sol_ndx = int(data[14])
            
            # 初始化排名
            self._init_ranks()
            return True
        except (ValueError, IndexError):
            return False
    
    def _init_ranks(self) -> None:
        """初始化所有排名为0"""
        for attr in ['rmsd', 'pvalue', 'DI_ALL', 'INF_ALL', 'INF_WC', 
                    'INF_NWC', 'INF_STACK', 'clashscore', 'mcq', 'tm', 'gdt']:
            setattr(self, f"{attr}_rank", 0)
    
    def file_name(self) -> str:
        """生成文件名"""
        return f"{self.problem}_{self.lab}_{self.result}"
    
    def set_rank(self, attr: str, rank: int) -> None:
        """设置指定属性的排名"""
        if hasattr(self, f"{attr}_rank"):
            setattr(self, f"{attr}_rank", rank)
        else:
            msgs.show("FATAL", f"Can't set attribute '{attr}' in class 'Eval'")
    
    def __str__(self) -> str:
        """字符串表示"""
        return " ".join([
            str(self.problem),
            str(self.original),
            self.lab,
            str(self.result),
            f"{self.rmsd:7.3f}",
            f"{self.pvalue:.3e}",
            f"{self.DI_ALL:7.3f}",
            f"{self.INF_ALL:7.3f}",
            f"{self.INF_WC:7.3f}",
            f"{self.INF_NWC:7.3f}",
            f"{self.INF_STACK:7.3f}",
            f"{self.clashscore:7.3f}",
            f"{self.mcq:7.3f}",
            f"{self.tm:7.3f}",
            f"{self.gdt:7.3f}",
            str(self.best_sol_ndx)
        ])
    
    def rankline(self) -> str:
        """生成排名行"""
        return "\t".join([
            self.lab.replace('PostExp','_PostExp').replace('PreExp','_PreExp'),
            str(self.result),
            f"{self.rmsd:.3f}",
            str(self.rmsd_rank),
            f"{self.DI_ALL:.3f}",
            str(self.DI_ALL_rank),
            f"{self.INF_ALL:.3f}",
            str(self.INF_ALL_rank),
            f"{self.INF_WC:.3f}",
            str(self.INF_WC_rank),
            f"{self.INF_NWC:.3f}",
            str(self.INF_NWC_rank),
            f"{self.INF_STACK:.3f}",
            str(self.INF_STACK_rank),
            f"{self.clashscore:.3f}",
            str(self.clashscore_rank),
            f"{self.mcq:.3f}",
            str(self.mcq_rank),
            f"{self.tm:.3f}",
            str(self.tm_rank),
            f"{self.gdt:.3f}",
            str(self.gdt_rank),
            f"{self.pvalue:.3e}",
            str(self.best_sol_ndx)
        ])

def find_eval(evals: List[Eval], problem: str, lab: str, result: int) -> Optional[Eval]:
    """查找匹配的评估结果"""
    return next((eval for eval in evals 
                 if eval.problem == problem and eval.lab == lab and eval.result == result), None)

def load_evals_list(fname: str) -> List[Eval]:
    """加载评估结果列表"""
    evals = []
    try:
        rows = FileHandler.read_file(fname).split("\n")
        for row in rows:
            eval = Eval()
            if not eval.parse(row):
                print(f"FATAL: Syntax error in evals file '{fname}' row '{row}'")
                continue
            eval.ok = True
            evals.append(eval)
    except Exception as e:
        print(f"FATAL: Error loading evals file '{fname}': {str(e)}")
    return evals

def save_evals_list(evals: List[Eval], fname: str) -> None:
    """保存评估结果列表"""
    try:
        with open(fname, "w") as f:
            for eval in filter(lambda e: e.ok, evals):
                f.write(f"{eval}\n")
    except Exception as e:
        print(f"FATAL: Error saving evals file '{fname}': {str(e)}")

def compute_evals_ranks(evals: List[Eval]) -> None:
    """计算评估结果的排名"""
    config = UtilsConfig()
    
    # 使用pandas进行排名计算
    for attr, reverse in config.rank_attributes.items():
        try:
            df = pd.Series([getattr(eval, attr) for eval in evals])
            ranks = df.rank(method='min', ascending=(~reverse))
            for eval, rank in zip(evals, ranks):
                eval.set_rank(attr, int(rank))
        except Exception as e:
            print(f"Warning: Error computing rank for {attr}: {str(e)}")

def sort_evals( evals, attr ):
    evals.sort( key=attrgetter(attr) )

def save_md(evals: List[Eval], fname: str) -> None:
    """保存Markdown格式的评估结果"""
    try:
        table = []
        for eval in evals:
            row = [
                f"<tr><td>{eval.lab}</td><td>{eval.result}</td>",
                f"<td>{eval.rmsd:.3f}</td>",
                f"<td>{eval.DI_ALL:.3f}</td>",
                f"<td>{eval.INF_ALL:.3f}</td>",
                f"<td>{eval.INF_WC:.3f if eval.INF_WC >= 0.0 else '-'}</td>",
                f"<td>{eval.INF_NWC:.3f if eval.INF_NWC >= 0.0 else '-'}</td>",
                f"<td>{eval.INF_STACK:.3f if eval.INF_STACK >= 0.0 else '-'}</td>",
                f"<td>{eval.clashscore:.3f}</td>",
                f"<td>{eval.pvalue:.2e}</td>",
                f"<td>{eval.mcq:.2f}</td>",
                f"<td>{eval.tm:.4f}</td>",
                f"<td>{eval.best_sol_ndx + 1}</td>",
                f"<td><a href='/show/index.html?id={eval.file_name()}'>-></a></td></tr>"
            ]
            table.append("".join(row))
        
        template = """---
layout: post
title: "Puzzle {problem} table"
date: 2000-01-01 00:00:01 +0000
comments: false
categories: 
- table
---

<script src="{{ root_url }}/javascripts/sorttable.js"></script>
<script>
    window.onload = function() {{
        (document.getElementsByTagName( 'th' )[1]).click();
    }};
</script>
<br/>
<div align="center">
Please click the tags to sort data.<br/>
<table class="sortable" border=1>
  <tr>
    <th>Lab</th>
    <th>Num</th>
    <th>RMSD</th>
    <th>DI all</th>
    <th>INF all</th>
    <th>INF wc</th>
    <th>INF nwc</th>
    <th>INF stacking</th>
    <th>Clash Score</th>
    <th>P-value</th>
    <th>mcq</th>
    <th>TM-score</th>
    <th>best sol.</th>
    <th>Detail</th>
  </tr>
  {table}
</table>
-: No true positive interaction in the predicted structure.
</div>
"""
        content = template.format(
            problem=evals[0].problem if evals else "",
            table="\n".join(table)
        )
        FileHandler.write_file(fname, content)
    except Exception as e:
        print(f"FATAL: Error saving markdown file '{fname}': {str(e)}")
