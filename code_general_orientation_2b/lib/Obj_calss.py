# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: Obj_calss.py
@date: 6/5/23 16:02
@desc: 
"""
import numpy as np
import sys
import traceback
from CustomError import CustomError
from RotationMatrix import RotationMatrix


# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


class Obj_class:
    """
    This is a OBJ file class
    """
    def __init__(self, mat_name = None,mat_value = None, vertices = None, faces = None):
        self.mat_name = mat_name
        self.mat_value = mat_value
        self.vertices = vertices
        self.faces = faces
        self.verNum = np.size(vertices, 0)       # Number of vertices
        self.facNum = np.size(faces, 0)      # Number of faces

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

    def objTranslation(self, dS, direction):
        """
        This function translates the polyhedron object in 'direction'
        :param dS: Translation length
        :param direction: Translation direction
        :return: translated vertices
        """
        if direction=='X':
            self.vertices = self.vertices + np.array([dS, 0, 0])
        elif direction == 'Y':
            self.vertices = self.vertices + np.array([0, dS, 0])
        elif direction == 'Z':
            self.vertices = self.vertices + np.array([0, 0, dS])
        else:
            self.showError("The translation direction must be:\n 'X', 'Y', or 'Z'")

    def objRotation(self, rad, axis):
        """
        This function rotates the polyhedron object
        :param rad: Rotation radian
        :param axis: Rotation axis
        :return: rotated vertices
        """
        ver = np.zeros([self.verNum,3])
        if axis == 'X' or 'Y' or 'Z':
            DCM = RotationMatrix(rad, axis)
        else:
            self.showError("The rotation axis must be:\n 'X', 'Y', or 'Z'")

        for i in range(self.verNum):
            ver[i,:] = DCM.dot(self.vertices[i,:])
        self.vertices = ver

    def showError(self, msg):
        """
        this function prints an error message and stops the execution the program.
        :param msg: the message to be shown
        :return: exit program
        """
        try:
            raise CustomError(msg)
        except CustomError as e:
            traceback.print_exc()
            sys.exit(-1)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    mat_name = "Material000"
    mat_value = {'Kd': [0.121569, 0.466667, 0.705882]}
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
