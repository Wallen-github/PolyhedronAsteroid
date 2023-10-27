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
import math
from pylmgc90 import chipy

def Set_Unit():
    GMA = 2.650         # m^3/s^2, Apophis gravitaional parameter
    GG = 6.67430e-11    # G, N*m^2/kg^2 = m^3/kg/s^2
    RA = 163            # m
    UnitM = GMA/GG
    UnitL = RA
    UnitT = np.sqrt(UnitL**3/GG/UnitM)

    return UnitM, UnitL, UnitT

def Get_TotalMomentOfInertia_com(nbR3):
    inertia_global = np.zeros([3, 3], 'd')
    coor = np.empty([nbR3, 6], dtype=float)
    for i in range(nbR3):
        inertia_global += chipy.RBDY3_GetGlobInertia(i+1)
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        density = chipy.RBDY3_GetBodyDensity(i+1)
        mass = chipy.RBDY3_GetMass(i+1)
        volume = mass/density
        for j in range(0, 3):
            d = np.array(coor[i, 0:3])
            d[i] = 0.
            inertia_global[j, j] += volume * np.dot(d, d)

        # contribution de la distance a l'axe aux termes extra-diagonaux
        d = np.array(coor[i, 0:3])
        inertia_global[0, 1] -= volume * d[0] * d[1]
        inertia_global[1, 0] -= volume * d[0] * d[1]
        inertia_global[0, 2] -= volume * d[0] * d[2]
        inertia_global[2, 0] -= volume * d[0] * d[2]
        inertia_global[1, 2] -= volume * d[1] * d[2]
        inertia_global[2, 1] -= volume * d[1] * d[2]

    nb = 0
    # pour chaque ligne
    for i in range(0, 3):
        # pour chaque terme au-dessus de la diagonale
        for j in range(i + 1, 3):
            # si le terme extra-diagonal courant est negligeable devant le terme
            # diagonal courant
            if math.fabs(inertia_global[i, j]) < 1.e-14 * math.fabs(inertia_global[i, i]):
                # on incremente le nombre de termes au-dessus de la diagonale
                # annulles
                nb += 1
                # on annule le terme extra-diagonal au-dessus de la diagonale
                inertia_global[i, j] = 0.e0
                # et le terme extra-diagonal au-dessous de la diagonale (matrice
                # symetrique)
                inertia_global[j, i] = 0.e0

    # diagonalisation de la matrice

    # si la matrice "nettoyee" est deja diagonale
    if nb == 3:
        # on obtient immediatement les valeurs propres (valeurs diagonales)
        I_diag = np.array([inertia_global[0, 0], inertia_global[1, 1], inertia_global[2, 2]])
        # et la matrice de passage (matrice identite)
        P = np.eye(3, 3)
    # sinon,
    else:
        # on appelle la routine de calcul des valeurs propres et vecteurs propres
        # disponible dans scipy
        I_diag, P = np.linalg.eigh(inertia_global)
        # on s'assure que le repere est direct en recalculant la troisieme direction
        # comme le produit vectoriel des deux premieres
        P[:, 2] = np.cross(P[:, 0], P[:, 1])

    return inertia_global, I_diag, P

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


