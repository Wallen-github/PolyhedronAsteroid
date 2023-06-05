# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: test.py
@date: 6/5/23 17:16
@desc: 
"""


def set_inf(name,age):
  if not 0 < age < 120:
    raise ValueError('超出范围')
  else:
    print('%s is %s years old' % (name,age))
def set_inf2(name,age):
  assert 0 < age < 120,'超出范围'
  print('%s is %s years old' % (name,age))
if __name__ == '__main__':
  try:
   set_inf('bob',200)
  except ValueError as e:
    print('无效值:',e)
  set_inf2('bob',200)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
