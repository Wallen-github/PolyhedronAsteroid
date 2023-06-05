# coding=utf-8
"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: readObj.py
@date: 6/5/23 13:43
@desc: 
"""
import sys
import numpy as np
from Obj_calss import Obj_class
# sys.path.append("../data")

def readObj(Objfile):

    try:
        file = open(Objfile, "r")
    except IOError:
        print('Error: The object file " ' + Objfile + ' " can not be found.')
    else:
        print('Reading Object file : '+ Objfile)
        OBJ = Objreader(file)
        file.close()

    return OBJ

def Objreader(file):

    ObjNum = 0
    nvs = 0
    OBJ_list = []
    while 1:
        vertices = []
        faces = []

        line = file.readline()
        line = line.strip('\n')

        if not line:
            print('Error: This OBJ file is empty.')
            break
        strs = line.split(" ")

        if strs[0] == "mtllib":
            matfilepath = strs[1]
            try:
                matfile = open(matfilepath, "r")
            except IOError:
                raise IOError('Error: The object file " ' + matfilepath + ' " can not be found.')
            else:
                print('Reading Object file : ' + matfilepath)
                mat = Matreader(matfilepath)
                matfile.close()

        if strs[0] == "usemtl":
            OBJ_list[ObjNum] = Obj_class(mat[0]["mat_name"],mat[0]["mat_value"], vertices, faces)
            nvs = nvs + verNum
            verNum = 0
        if strs[0] == "v":
            vertices.append(np.array([float(strs[1]), float(strs[2]), float(strs[3])]))
            OBJ_list[ObjNum].vertices = vertices
            verNum += 1
        if strs[0] == "f":
            faces.append(np.array([int(strs[1])-nvs, int(strs[2])-nvs, int(strs[3])-nvs]))
            OBJ_list[ObjNum].faces = faces

    return OBJ_list

def Matreader(matfile):
    """
    This function is used to read the material file (.mtl)
    :param  matfile: the opened material file (.mtl)
    :return: Mat: the material list, each element is a dictionary for one object.
                    Two keys exist in each dictionary,
                        'mat_name': name of material, str
                        'kd': diffuse color, array(1,3)
                    for example: the second object
                        Mat[1]["mat_name"] = "Material000"
                        Mat[1]["kd"] = [1,2,3]
    """

    Mat = []
    while 1:
        line = matfile.readline()
        line = line.strip('\n')

        if not line:
            print('Error: This mtllib file is empty.')
            break

        strs = line.split(" ")
        if strs[0] == "newmtl":
            mat_name = strs[1]
            Mat.append({"mat_name": mat_name})
        if strs[0] == "kd":
            Mat[-1]["kd"] = [float(strs[1]),float(strs[2]),float(strs[3])]

    return Mat


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    Objfile = "../data/Apophis_N3_v200.obj"
    readObj(Objfile)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
