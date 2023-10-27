# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: test_notification.py
@date: 9/22/23 21:43
@desc: 
"""


# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


import requests
import os

#获取文件名（含后缀）
name=os.path.basename(__file__)
path = os.getcwd() # 获取当前工作目录路径
mesg = path+'/'+name
requests.post("https://ntfy.sh/wallen1732_notfy_python",
            data=mesg.encode(encoding='utf-8'))




