import sys
from enum import Enum
from typing import Optional

class LogLevel(Enum):
    """日志级别枚举"""
    DEBUG = "DEBUG"
    INFO = "INFO" 
    WARNING = "WARNING"
    ERROR = "ERROR"
    FATAL = "FATAL"

class MessageHandler:
    """消息处理类"""
    
    def __init__(self):
        """初始化消息处理器"""
        self._min_level = LogLevel.INFO
        self._max_prefix_length = 10
        
    def set_min_level(self, level: LogLevel) -> None:
        """设置最小日志级别
        
        Args:
            level: 最小日志级别
        """
        self._min_level = level
        
    def show(self, 
             prefix: str, 
             message: str, 
             new_line: bool = True, 
             back: bool = False,
             level: Optional[LogLevel] = None) -> None:
        """显示消息
        
        Args:
            prefix: 消息前缀
            message: 消息内容
            new_line: 是否换行
            back: 是否回退
            level: 日志级别,如果为None则从prefix推断
        """
        # 确定日志级别
        if level is None:
            try:
                level = LogLevel(prefix.upper())
            except ValueError:
                level = LogLevel.INFO
                
        # 检查是否达到最小日志级别
        if level.value < self._min_level.value:
            return
            
        # 格式化消息
        new_line_txt = "\n" if new_line else ""
        back_txt = "\b" * 50 if back else ""
        prefix = prefix[:self._max_prefix_length]
        
        # 输出消息
        sys.stderr.write(
            f"{back_txt}{prefix:<{self._max_prefix_length}} >> {message}{new_line_txt}"
        )
        
        # 如果是回退模式,立即刷新
        if back:
            sys.stderr.flush()
            
        # 如果是致命错误,终止程序
        if level == LogLevel.FATAL:
            sys.stderr.write("Abrupt termination!\n")
            sys.exit(1)
            
    def debug(self, message: str, new_line: bool = True) -> None:
        """显示调试消息"""
        self.show("DEBUG", message, new_line, level=LogLevel.DEBUG)
        
    def info(self, message: str, new_line: bool = True) -> None:
        """显示信息消息"""
        self.show("INFO", message, new_line, level=LogLevel.INFO)
        
    def warning(self, message: str, new_line: bool = True) -> None:
        """显示警告消息"""
        self.show("WARNING", message, new_line, level=LogLevel.WARNING)
        
    def error(self, message: str, new_line: bool = True) -> None:
        """显示错误消息"""
        self.show("ERROR", message, new_line, level=LogLevel.ERROR)
        
    def fatal(self, message: str, new_line: bool = True) -> None:
        """显示致命错误消息并终止程序"""
        self.show("FATAL", message, new_line, level=LogLevel.FATAL)

# 创建全局消息处理器实例
msgs = MessageHandler()

# 为了保持向后兼容,保留原来的show函数
def show(prefix: str, message: str, new_line: bool = True, back: bool = False) -> None:
    """向后兼容的show函数
    
    Args:
        prefix: 消息前缀
        message: 消息内容
        new_line: 是否换行
        back: 是否回退
    """
    msgs.show(prefix, message, new_line, back)