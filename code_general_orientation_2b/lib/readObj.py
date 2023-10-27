# coding=utf-8
"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: readObj.py
@date: 6/5/23 13:43
@desc: 
"""
import numpy as np
from Obj_calss import Obj_class

def readObj(Objpath,Objfile):
    """
    This function is a API of Objreader.
    :param Objpath: files path, str
    :param Objfile: files name, str
    :return: OBJ --- polyhedron objects list, list
                OBJ[i] --- The i-th polyhedron object.
                OBJ[i].mat_name --- material name, str.
                OBJ[i].mat_values["Kd"] --- diffuse color, array([1,3]).
                OBJ[i].vertices --- Polyhedron vertices, array([N,3]).
                OBJ[i].faces --- Polyhedron faces, array([N,3]).
    """
    try:
        file = open(Objpath + Objfile, "r")
    except IOError:
        raise IOError('Error: The object file " ' + Objfile + ' " can not be found.')
    else:
        print('Reading Object file : '+ Objfile)
        OBJ = Objreader(Objpath,file)
        file.close()

    return OBJ

def Objreader(Objpath,file):
    """
    The function read and rearrange the Obj file from V-HACD algothrim.
    :param Objpath: files path, str
    :param file: file object, class
    :return: OBJ_list --- polyhedron objects list, list
                OBJ_list[i] --- The i-th polyhedron object.
                OBJ_list[i].mat_name --- material name, str.
                OBJ_list[i].mat_values["Kd"] --- diffuse color, array([1,3]).
                OBJ_list[i].vertices --- Polyhedron vertices, array([N,3]).
                OBJ_list[i].faces --- Polyhedron faces, array([N,3]).
    """
    ObjNum = 0
    nvs = 0
    verNum = 0
    OBJ_list = []
    while 1:
        vertices = []
        faces = []

        line = file.readline()
        line = line.strip('\n')

        if not line:
            print('Finish read OBJ file.')
            break
        strs = line.split(" ")

        if strs[0] == "mtllib":
            matfilepath = Objpath+strs[1]
            try:
                matfile = open(matfilepath, "r")
            except IOError:
                raise IOError('Error: The object file " ' + matfilepath + ' " can not be found.')
            else:
                print('Reading Object file : ' + matfilepath)
                mat = Matreader(matfile)
                matfile.close()

        if strs[0] == "usemtl":
            ObjNum = ObjNum+1
            OBJ_list.append(Obj_class(mat[0]["mat_name"],mat[0]["Kd"], vertices, faces))
            nvs = nvs + verNum
            verNum = 0
            facNum = 0
        if strs[0] == "v":
            vertices=[float(strs[1]), float(strs[2]), float(strs[3])]
            OBJ_list[ObjNum-1].vertices.append(vertices)
            verNum += 1
        if strs[0] == "f":
            faces=[int(strs[1])-nvs, int(strs[2])-nvs, int(strs[3])-nvs]
            OBJ_list[ObjNum-1].faces.append(faces)
            facNum += 1

    for i in range(ObjNum):
        OBJ_list[i].vertices = np.array(OBJ_list[i].vertices)
        OBJ_list[i].faces = np.array(OBJ_list[i].faces)

    return OBJ_list

def Matreader(matfile):
    """
    This function is used to read the material file (.mtl)
    :param  matfile: the opened material file (.mtl)
    :return: Mat --- the material list, each element is a dictionary for one object.
                    Two keys exist in each dictionary,
                        'mat_name' --- name of material, str
                        'kd': diffuse color, array(1,3)
                    for example --- the second object
                        Mat[1]["mat_name"] = "Material000"
                        Mat[1]["kd"] = [1,2,3]
    """

    Mat = []
    while 1:
        line = matfile.readline()
        line = line.strip('\n')

        if not line:
            print('Finish read mtllib file.')
            break

        strs = line.split(" ")
        if strs[0] == "newmtl":
            mat_name = strs[1]
            Mat.append({"mat_name": mat_name})
        if strs[0] == "Kd":
            Mat[-1]["Kd"] = np.array([float(strs[1]),float(strs[2]),float(strs[3])])

    return Mat


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    Objpath = "../data/"
    Objfile = "Apophis_N3_v200.obj"
    OBJ = readObj(Objpath,Objfile)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
