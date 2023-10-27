# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: InitialSet.py
@date: 7/10/23 14:22
@desc: 
"""

from pylmgc90.pre import *
from lib.RotationMatrix import *


class ApophisPram:
    # Apophis Parameters
    GravParam = 2.605  # m^3/s^2
    GG = 6.6742e-11  # G, N*m^2/kg^2 = m^3/kg/s^2
    Density = 2000  # kg/m^3
    longitude = 247 # deg
    latitude = -59 # deg

    # Normalized unit
    UnitM = GravParam / GG
    UnitL = (UnitM / Density) ** (1 / 3)
    UnitT = np.sqrt(UnitL ** 3 / GG / UnitM)


def InitialPolySet(rho,dim,OBJ):

    mat = material(name='STONE', materialType='RIGID', density=rho)
    mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

    poly_list,COM = normalizePoly(OBJ,rho, mat, mod)
    poly_list_new = GetPoly(poly_list, COM)

    return poly_list_new, mat, mod

def normalizePoly(OBJ,rho,mat,mod):

    mass_total = 0
    volume_total = 0
    for i in range(len(OBJ)):
        vertices = OBJ[i].vertices
        faces = OBJ[i].faces
        poly = rigidPolyhedron(model=mod, material=mat, color='BLEUx',
                               generation_type='full', vertices=vertices,
                               faces=faces, tol=0., number=None, seed=None,
                               xr=1., yr=1., zr=1.)
        mass_total += poly.contactors[0].volume * rho
        volume_total += poly.contactors[0].volume
        del poly,vertices,faces

    # Creat poly list and find COM
    Lscale = (volume_total) ** (1 / 3)
    mass_total = 0
    volume_total = 0
    COM = np.zeros([3])
    poly_list = []
    for i in range(len(OBJ)):
        # vertices = OBJ[i].vertices/Lscale*ApophisPram.UnitL
        vertices = OBJ[i].vertices / Lscale
        faces = OBJ[i].faces
        poly = rigidPolyhedron(model=mod, material=mat, color='BLEUx',
                               generation_type='full', vertices=vertices,
                               faces=faces, tol=0., number=None, seed=None,
                               xr=1., yr=1., zr=1.)
        mass_total += poly.contactors[0].volume * rho
        volume_total += poly.contactors[0].volume
        COM += poly.contactors[0].volume * rho * poly.nodes[1].coor
        poly_list.append(poly)
        del poly,vertices,faces

    # Center of Mass
    COM = COM / mass_total

    return poly_list,COM

def GetPoly(poly_list, COM):

    # Euler angle
    phi = np.deg2rad(0)
    theta = np.deg2rad(15.208605645340988)
    psi = np.deg2rad(90)
    # longitude and latitude
    alpha = np.deg2rad(ApophisPram.longitude)
    beta = np.deg2rad(ApophisPram.latitude)
    # Rotation matrix
    RM1 = RotationMatrix(phi, 'Z')
    RM2 = RotationMatrix(theta, 'X')
    RM3 = RotationMatrix(psi, 'Z')
    RMb = RotationMatrix(np.pi / 2 - beta, 'Y')
    RMa = RotationMatrix(alpha, 'Z')
    BN = RM3 @ RM2 @ RM1 @ RMb @ RMa
    # Spin rate
    # omegaB_BN = [0.202078195115783E-4, 0, 0.542630931734589E-4]  # rad/sec
    # omegaB_BN = [i *ApophisPram.UnitT for i in omegaB_BN]
    omegaB_BN = [0,0,2*np.pi]
    omegaN_BN = BN.T @ omegaB_BN
    vel = [0, 0, 0, omegaN_BN[0], omegaN_BN[1], omegaN_BN[2]]
    # vel = [i * 10000000 for i in vel]

    poly_list_new = []
    for i in range(len(poly_list)):
        poly = poly_list[i]
        # Translate to COM
        poly.translate(-COM[0], -COM[1], -COM[2])
        # Rotate to initial attitude, 3-1-3 rotation
        poly.rotate(description='axis', alpha=alpha, axis=[0., 0., 1.])
        poly.rotate(description='axis', alpha=np.pi / 2 - beta, axis=[0., 1., 0.])
        poly.rotate(description='axis', alpha=phi, axis=[0., 0., 1.])
        poly.rotate(description='axis', alpha=theta, axis=[1., 0., 0.])
        poly.rotate(description='axis', alpha=psi, axis=[0., 0., 1.])
        # Set initial spin rate
        poly.imposeInitValue(component=[1, 2, 3, 4, 5, 6], value=vel)
        poly_list_new.append(poly)

    return poly_list_new



# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    print(ApophisPram.GravParam)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
