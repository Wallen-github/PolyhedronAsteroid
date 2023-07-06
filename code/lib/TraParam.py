# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: TraParam.py
@date: 6/8/23 16:33
@desc: 
"""
import sys
sys.path.append("..")
import numpy as np
from pylmgc90 import chipy

def Get_EnergyMomentum(nbR3,GG=1):

    momentum = np.empty([1,4], dtype=float)
    mass = np.empty([nbR3], dtype=float)
    coor = np.empty([nbR3, 6], dtype=float)
    vel = np.empty([nbR3, 6], dtype=float)
    vel_cm = np.zeros([3], dtype=float)
    Energy = 0
    Knetic = 0
    potent = 0

    for i in range(nbR3):
        coor[i,:] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        mass[i] = chipy.RBDY3_GetMass(i + 1)
        vel[i,:] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)

    Msum = np.sum(mass)
    for i in range(nbR3 - 1):
        for j in range(i + 1, nbR3):
            vij = vel[j, 0:3] - vel[i, 0:3]
            rij = coor[j, 0:3] - coor[i, 0:3]
            Knetic += mass[i] * mass[j] * np.dot(vij, vij) / (2. * Msum)
            momentum[0,0:3] += mass[i] * mass[j] * (np.cross(rij, vij)) / Msum
            potent += - GG * mass[i] * mass[j] / np.linalg.norm(rij) #gravpot[i]  #

    for i in range(nbR3):
        vel_cm += mass[i] * vel[i, 0:3] / Msum
        inertia = chipy.RBDY3_GetGlobInertia(i + 1)
        omega_i = vel[i, 3:6]
        Knetic += 0.5 * np.dot(omega_i,np.dot(inertia,omega_i))
        momentum[0,0:3] += np.dot(inertia,omega_i)
    Energy = Knetic + potent
    momentum[0,3] = np.dot(momentum[0,0:3],momentum[0,0:3])
    return Energy, Knetic, momentum
