# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: SubFunc.py
@date: 7/18/23 13:07
@desc: 
"""

# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import sys
sys.path.append("..")
import numpy as np
import math
from pylmgc90 import chipy

def Get_TotalMomentOfInertia_gen_v2(poly_list):
    nbR3 =  len(poly_list)
    Ic = np.zeros([3, 3], dtype=np.float32)
    MR2 = np.zeros([3, 3], dtype=np.float32)
    I_global = np.zeros([3, 3], dtype=np.float32)
    coor = np.zeros([nbR3, 3], dtype=np.float32)
    mass = np.zeros(nbR3, dtype=np.float32)

    for i in range(nbR3):
        poly0 = poly_list[i]
        mass[i] = poly0.contactors[0].volume
        Ic = poly0.contactors[0].I
        coor[i, :] = poly0.nodes[1].coor
        d = np.array(coor[i, 0:3])
        MR2 = np.array([[d[2]**2 + d[1]**2, -d[0]*d[1], -d[0]*d[2]],
                        [-d[0]*d[1], d[2]**2 + d[0]**2, -d[1]*d[2]],
                        [-d[0]*d[2], -d[1]*d[2], d[1]**2 + d[0]**2]])
        I_global += MR2 + Ic

    I_diag, P = np.linalg.eigh(I_global)
    P[:, 2] = np.cross(P[:, 0], P[:, 1])
    return I_global, I_diag, P
def Get_TotalMomentOfInertia_gen(poly_list):
    '''

    :param poly_list:
    :return:
    '''

    inertia_total = np.zeros([3, 3], 'd')
    for i in range(len(poly_list)):
        poly0 = poly_list[i]
        inertia_total += poly0.contactors[0].I
        for j in range(0, 3):
            d = np.array(poly0.nodes[1].coor)
            # d = np.array(poly0.contactors[0].shift)
            d[j] = 0.
            inertia_total[j, j] += poly0.contactors[0].volume * np.dot(d, d)

        # contribution de la distance a l'axe aux termes extra-diagonaux
        d = poly0.nodes[1].coor
        inertia_total[0, 1] -= poly0.contactors[0].volume * d[0] * d[1]
        inertia_total[1, 0] -= poly0.contactors[0].volume * d[0] * d[1]
        inertia_total[0, 2] -= poly0.contactors[0].volume * d[0] * d[2]
        inertia_total[2, 0] -= poly0.contactors[0].volume * d[0] * d[2]
        inertia_total[1, 2] -= poly0.contactors[0].volume * d[1] * d[2]
        inertia_total[2, 1] -= poly0.contactors[0].volume * d[1] * d[2]

    # diagonalisation de la matrice d'inertie

    # nettoyage de la matrice :

    # on initialise le nombre de termes au-dessus de la diagonale annulles
    nb = 0
    # pour chaque ligne
    for i in range(0, 3):
        # pour chaque terme au-dessus de la diagonale
        for j in range(i + 1, 3):
            # si le terme extra-diagonal courant est negligeable devant le terme
            # diagonal courant
            if math.fabs(inertia_total[i, j]) < 1.e-14 * math.fabs(inertia_total[i, i]):
                # on incremente le nombre de termes au-dessus de la diagonale
                # annulles
                nb += 1
                # on annule le terme extra-diagonal au-dessus de la diagonale
                inertia_total[i, j] = 0.e0
                # et le terme extra-diagonal au-dessous de la diagonale (matrice
                # symetrique)
                inertia_total[j, i] = 0.e0

    # diagonalisation de la matrice

    # si la matrice "nettoyee" est deja diagonale
    if nb == 3:
        # on obtient immediatement les valeurs propres (valeurs diagonales)
        I_diag = np.array([inertia_total[0, 0], inertia_total[1, 1], inertia_total[2, 2]])
        # et la matrice de passage (matrice identite)
        P = np.eye(3, 3)
    # sinon,
    else:
        # on appelle la routine de calcul des valeurs propres et vecteurs propres
        # disponible dans scipy
        I_diag, P = np.linalg.eigh(inertia_total)
        # on s'assure que le repere est direct en recalculant la troisieme direction
        # comme le produit vectoriel des deux premieres
        P[:, 2] = np.cross(P[:, 0], P[:, 1])

    return inertia_total, I_diag, P

def DCM2EA_313(DCM):

    theta1 = math.atan2(DCM[2,0],-DCM[2,1])
    if (DCM[2,2] - 1) < 1e-5:
        DCM[2,2] = 1
    theta2 = math.acos(DCM[2,2])
    theta3 = math.atan2(DCM[0,2],DCM[1,2])

    return theta1,theta2,theta3

def Set_Unit():
    GMA = 2.650         # m^3/s^2, Apophis gravitaional parameter
    GG = 6.67430e-11    # G, N*m^2/kg^2 = m^3/kg/s^2
    RA = 168            # m
    UnitM = GMA/GG
    UnitL = RA
    UnitT = np.sqrt(UnitL**3/GG/UnitM)
    UnitG = 1.

    return UnitM, UnitL, UnitT, UnitG

def Get_CenterMass(coor,mass):
    '''
    This function computes the center of mass
    Args:
        coor (numpy array): positions in the inertial frame
        mass (numpy array): mass

    Returns:
        CM center of mass
        coor_cm coordinate origin at center of mass
    '''
    CM = np.zeros([3])
    coor_cm = np.zeros([len(mass),3])
    for i in range(len(mass)):
        CM += mass[i]*coor[i,:]
    CM = CM/np.sum(mass)

    for i in range(len(mass)):
        coor_cm[i,:] = coor[i,:] - CM

    return CM,coor_cm

def Detect_State(nbR3, max_mass_index, relative_distance, max_radius):

    # Get hill radiu, minimum is 0.49606892, which is not work
    # UnitM, UnitL, UnitT, UnitG = Set_Unit()
    # massE = 5.974227245203837E24 / UnitM
    # mass = chipy.RBDY3_GetMass(max_mass_index+1)
    # R_hill = np.linalg.norm(posE)*(mass/3/massE)**(1/3)
    Hill_radiu = 5
    loss_num = 0
    nomove_num = 0
    shifting_num = 0
    circular_num = 0
    info = 1

    for i in range(nbR3):
        if i != max_mass_index:
            dis = max_radius[max_mass_index] + max_radius[i]/2
            if relative_distance[i]>Hill_radiu:
                loss_num += 1
                if loss_num >= nbR3-1:
                    info = 0
            elif relative_distance[i] < 1e-4:
                nomove_num += 1
            elif relative_distance[i] > 1e-4 and relative_distance[i] < dis:
                shifting_num += 1
            elif relative_distance[i] > dis and relative_distance[i] < Hill_radiu:
                circular_num += 1

    return info, loss_num, nomove_num, shifting_num, circular_num

def Get_RelativeDistance(nbR3, mass, coor):
    '''
    This function compute the relative distance between the max-mass body and other bodies.
    :param nbR3: number of bodies
    :param mass: mass list of bodies
    :param coor: position list of bodies
    :return: R2, list of relative distance
    '''

    R2 = np.zeros(nbR3, dtype=np.float32)
    # Find the max-mass object
    max_mass_index = np.argmax(mass)

    for i in range(nbR3):
        R2[i] = np.linalg.norm(coor[i, 0:3] - coor[max_mass_index, 0:3])

    return R2

def Get_RelativeAngel(nbR3, mass, coor, NB_pre):
    """
    This function aims to get the angles between body positions and the body-fixed vector of max-mass body.
    The body-fixed vector is the maximum inertial principal direction
    :param nbR3: number or bodies
    :param mass: mass list of bodies
    :param coor: position list of bodies
    :param NB_pre: The previous NB matrix
    :return: theta_list, angle list
    """

    # Find the max-mass object
    max_mass_index = np.argmax(mass)

    # choose the max eigenvalue direction as a body-fixed vector
    I_global = chipy.RBDY3_GetGlobInertia(int(max_mass_index+1))
    I_diag, NB = np.linalg.eigh(I_global)
    NB[:, 2] = np.cross(NB[:, 0], NB[:, 1])
    # detect if the eigenvector is opposite direction
    for i in range(3):
        dis = np.linalg.norm(NB[:,i] - NB_pre[:,i])
        # if the distance betwee previous eigenvector and current eigenvector larger than 1.1,
        # the current eigenvalue is point to opposite direction, then I correct it here.
        # (because the initial NB_pre is [0,0,0] and the dis ~ 1, 1.1 would be a suitable value.)
        if dis > np.sqrt(2):
            NB[:,i] = -NB[:,i]
    NB_pre = NB
    max_eig_index = np.argmax(I_diag)
    vector_bodyfix = NB[:,max_eig_index]

    theta_list = []
    for i in range(nbR3):
        if i != max_mass_index:
            relativePos = np.dot(NB.T,coor[i,0:3] - coor[max_mass_index,0:3])
            R2 = coor[i,0:3] - coor[max_mass_index,0:3]
            costheta = np.dot(vector_bodyfix,R2)/np.linalg.norm(vector_bodyfix)/np.linalg.norm(R2)
            theta = math.acos(costheta)
            theta_list.append(theta)

    return np.array(theta_list), NB_pre, relativePos

def Get_TotalMomentOfInertia_v2(nbR3):
    Ic = np.zeros([3, 3], dtype=np.float32)
    MR2 = np.zeros([3, 3], dtype=np.float32)
    I_global = np.zeros([3, 3], dtype=np.float32)
    coor = np.zeros([nbR3, 6], dtype=np.float32)
    mass = np.zeros(nbR3, dtype=np.float32)

    for i in range(nbR3):
        mass[i] = chipy.RBDY3_GetMass(i + 1)
        Ic = chipy.RBDY3_GetBodyInertia(i+1)
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        d = np.array(coor[i, 0:3])
        MR2 = np.array([[d[2]**2 + d[1]**2, -d[0]*d[1], -d[0]*d[2]],
                        [-d[0]*d[1], d[2]**2 + d[0]**2, -d[1]*d[2]],
                        [-d[0]*d[2], -d[1]*d[2], d[1]**2 + d[0]**2]])
        I_global += MR2 + Ic

    I_diag, P = np.linalg.eigh(I_global)
    P[:, 2] = np.cross(P[:, 0], P[:, 1])
    return I_global, I_diag, P
def Get_TotalMomentOfInertia_com(nbR3):
    inertia_global = np.zeros([3, 3], dtype=np.float32)
    coor = np.zeros([nbR3, 6], dtype=np.float32)
    for i in range(nbR3):
        inertia_global += chipy.RBDY3_GetGlobInertia(i+1)
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        volume = chipy.RBDY3_GetVolume(i+1)
        for j in range(0, 3):
            d = np.array(coor[i, 0:3])
            d[j] = 0.
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

    Momentum = np.zeros([4], dtype=np.float32)
    mass = np.zeros([nbR3], dtype=np.float32)
    coor = np.zeros([nbR3, 6], dtype=np.float32)
    vel = np.zeros([nbR3, 6], dtype=np.float32)
    vel_cm = np.zeros([3], dtype=np.float32)
    Energy = 0
    Kinetic = 0
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
            Kinetic += mass[i] * mass[j] * np.dot(vij, vij) / (2. * Msum)
            Momentum[0:3] += mass[i] * mass[j] * (np.cross(rij, vij)) / Msum
            potent += - GG * mass[i] * mass[j] / np.linalg.norm(rij) #gravpot[i]  #

    for i in range(nbR3):
        vel_cm += mass[i] * vel[i, 0:3] / Msum
        inertia = chipy.RBDY3_GetGlobInertia(i + 1)
        omega_i = vel[i, 3:6]
        Kinetic += 0.5 * np.dot(omega_i,np.dot(inertia,omega_i))
        Momentum[0:3] += np.dot(inertia,omega_i)
    Energy = Kinetic + potent
    Momentum[3] = np.dot(Momentum[0:3],Momentum[0:3])
    return Energy, Kinetic, Momentum

def Get_EffectiveInertiaSpin(Kinetic, Momentum):

    I_d = Momentum[3]/2/Kinetic
    w_e = 2*Kinetic/np.sqrt(Momentum[3])

    return I_d, w_e



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print('Hi PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
