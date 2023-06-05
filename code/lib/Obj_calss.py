# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: Obj_calss.py
@date: 6/5/23 16:02
@desc: 
"""
import numpy as np


# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


class Obj_class:
    """
    This is a OBJ file class
    """
    def __init__(self, mat_name,mat_value, vertices, faces):
        self.mat_name = mat_name
        self.mat_value = mat_value
        self.vertices = vertices
        self.faces = faces

    def displayObjInfo(self):
        """
        This function can display the OBJ class information to terminal
        Returns:

        """
        print("-----------OBJ file information display-------------")
        print("material name: ",self.mat_name)
        print("material properties: ", self.mat_value.keys())
        print("vertices size: ", self.vertices.shape)
        print("faces size: ", self.faces.shape)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    mat_name = "Material000"
    mat_value = {'kd': [0.121569, 0.466667, 0.705882]}
    vertices1 = np.array([[0,50,10],[10,50,0],[0,60,0],[0,50,0]])
    faces1 = np.array([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
    vertices2 = np.array([[0, 0, 10], [10, 10, 0], [-10, 10, 0], [0, -10, 0], [0, 0, 0]])
    faces2 = np.array([[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 5], [2, 4, 5]])

    OBJ = []
    OBJ.append(Obj_class(mat_name,mat_value, vertices1, faces1))
    OBJ.append(Obj_class(mat_name, mat_value, vertices2, faces2))
    OBJ[0].displayObjInfo()
    OBJ[1].displayObjInfo()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
