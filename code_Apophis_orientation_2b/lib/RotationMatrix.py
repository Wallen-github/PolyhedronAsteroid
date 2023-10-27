# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: RotationMatrix.py
@date: 6/8/23 21:04
@desc: 
"""


import numpy as np

def RotationMatrix(theta,axis):
    if axis == 'X':
        RotationX = np.array([[1, 0, 0],
                              [0, np.cos(theta), np.sin(theta)],
                              [0, -np.sin(theta), np.cos(theta)]])
        return RotationX
    elif axis == 'Y':
        RotationY = np.array([[np.cos(theta), 0, -np.sin(theta)],
                              [0, 1, 0],
                              [np.sin(theta), 0, np.cos(theta)]])
        return RotationY
    elif axis == 'Z':
        RotationZ = np.array([[np.cos(theta), np.sin(theta), 0],
                              [-np.sin(theta), np.cos(theta), 0],
                              [0, 0, 1]])
        return RotationZ
    else:
        print('error in func RotationMatrix: the "axis" should be one of X,Y,Z.')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    r1 = np.array([0,1,0])

    RM1 = RotationMatrix(90 * np.pi / 180, 'Z')
    RM2 = RotationMatrix(90 * np.pi / 180, 'X')
    RM3 = RotationMatrix(90 * np.pi / 180, 'Z')

    r2 = np.dot(RM1, r1)
    print('r2 = ', r2)

    r3 = np.dot(RM2, r2)
    print('r3 = ', r3)

    r4 = np.dot(RM3, r2)
    print('r4 = ', r4)

    DCM = np.dot(RM3, np.dot(RM2, RM1))
    r5 = np.dot(DCM, r1)
    print('r5 = ', r5)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
