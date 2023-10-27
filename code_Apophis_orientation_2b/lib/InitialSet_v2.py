# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: InitialSet_v2.py
@date: 7/16/23 23:01
@desc: 
"""


# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import math

def Get_TotalMomentOfInertia(poly_list):
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
            d[i] = 0.
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
    theta2 = math.acos(DCM[2,2])
    theta3 = math.atan2(DCM[0,2],DCM[1,2])

    return theta1,theta2,theta3



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
