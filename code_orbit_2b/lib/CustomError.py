# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: CustomError.py
@date: 6/8/23 21:58
@desc: 
"""

import traceback
import sys

class CustomError(Exception):
    def __init__(self, ErrorInfo):
        super().__init__(self)  # 初始化父类
        self.errorinfo = ErrorInfo

    def __str__(self):
        return self.errorinfo


if __name__ == '__main__':
    msg = "XYZ"
    try:
        raise CustomError(msg)
    except CustomError as e:
        traceback.print_exc()
        sys.exit(-1)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
