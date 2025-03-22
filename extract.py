#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PDB文件残基提取工具

此脚本用于从PDB文件中提取指定的残基。支持多链和多残基范围的提取。

用法:
    python extract.py <input model> <residue list> <output model>

参数:
    input model: 输入的PDB文件
    residue list: 要提取的残基列表,格式为'chain:res_id:count,...,chain:res_id:count'
    output model: 输出的PDB文件名

示例:
    python extract.py INPUT.pdb A:1:10,A:21:5,B:2:9 OUTPUT.pdb
"""

import sys
from typing import List, Tuple, Optional
from dataclasses import dataclass
from Bio.PDB import PDBParser, PDBIO, Select

@dataclass
class ResidueRange:
    """残基范围数据类"""
    chain: str
    res_id: int
    count: int
    
    @classmethod
    def from_string(cls, data: str) -> 'ResidueRange':
        """从字符串创建残基范围
        
        Args:
            data: 格式为"chain:res_id:count"的字符串
            
        Returns:
            ResidueRange实例
            
        Raises:
            ValueError: 如果数据格式不正确
        """
        try:
            chain, res_id, count = data.split(":")
            return cls(chain, int(res_id), int(count))
        except ValueError as e:
            raise ValueError(f"Invalid residue range format: {data}") from e

class ResidueSelector(Select):
    """残基选择器类"""
    
    def __init__(self):
        """初始化选择器"""
        self.res_list: List[str] = []
        
    def config(self, res_list: List[ResidueRange]) -> None:
        """配置要选择的残基
        
        Args:
            res_list: 残基范围列表
        """
        self.res_list = []
        for res_range in res_list:
            for i in range(res_range.count):
                self.res_list.append(f"{res_range.chain}|{res_range.res_id + i}")
                
    def accept_residue(self, residue: Any) -> bool:
        """判断是否接受该残基
        
        Args:
            residue: 残基对象
            
        Returns:
            是否接受该残基
        """
        key = f"{residue.parent.id.strip()}|{residue.get_id()[1]}"
        return key in self.res_list

def write_pdb(struct: Any, select_class: ResidueSelector, file: str) -> None:
    """写入PDB文件
    
    Args:
        struct: 结构对象
        select_class: 残基选择器
        file: 输出文件名
    """
    io = PDBIO()
    io.set_structure(struct)
    io.save(file, select_class)

def parse_res_list(s: str) -> List[ResidueRange]:
    """解析残基列表字符串
    
    Args:
        s: 残基列表字符串
        
    Returns:
        残基范围列表
        
    Raises:
        ValueError: 如果数据格式不正确
    """
    res_list = []
    for piece in s.split(","):
        try:
            res_list.append(ResidueRange.from_string(piece))
        except ValueError as e:
            msgs.fatal(f"Invalid residue data: {piece}")
    return res_list

def extract_pdb(input_file: str, res_list_str: str, output_file: str) -> None:
    """从PDB文件中提取残基
    
    Args:
        input_file: 输入PDB文件
        res_list_str: 残基列表字符串
        output_file: 输出PDB文件
    """
    try:
        # 解析输入结构
        parser = PDBParser()
        structure = parser.get_structure("SI", input_file)
        
        # 准备选择器
        selector = ResidueSelector()
        res_list = parse_res_list(res_list_str)
        selector.config(res_list)
        
        # 写入输出文件
        write_pdb(structure, selector, output_file)
        
    except Exception as e:
        msgs.fatal(f"Error processing PDB file: {str(e)}")

def print_usage() -> None:
    """打印使用说明"""
    print("\n" + "- " * 40)
    print("extract.py - extracts the selected residues from the input PDB file")
    print("- " * 40 + "\n")
    print("Usage:")
    print("$ python extract.py <input model> <residue list> <output model>\n")
    print("<input model> - Initial model in PDB format.")
    print("<residue list> - Residues to extract")
    print("<output model> - Name of the new file\n")
    print("Residue lists should be in the following format:\n'chain:res_id:count,...,chain:res_id:count'\n")
    print("Examples:")
    print("$ python extract.py INPUT.pdb A:1:10,A:21:5,B:2:9 OUTPUT.pdb")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print_usage()
        sys.exit(1)
        
    extract_pdb(sys.argv[1], sys.argv[2], sys.argv[3])
	

