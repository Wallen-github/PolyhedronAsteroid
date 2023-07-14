# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: gen_sample.py
@date: 1/10/23 13:38
@desc: 
"""

import sys
from lib.RotationMatrix import *
sys.path.append('./lib')
import numpy as np
from pylmgc90.pre import *
from lib.readObj import readObj


if not os.path.isdir('./DATBOX'):
    os.mkdir('./DATBOX')
if 'norand' in sys.argv:
    seed = list(range(13))
else:
    seed = None

# Set dimension
dim = 3
# Set uniform density
rho =1

Omega = 3e-5

# Creat empty class
bodies = avatars()
mat    = materials()
svs    = see_tables()
tacts  = tact_behavs()
stone = material(name='STONE', materialType='RIGID', density=rho)
mat.addMaterial(stone)
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

# Read OBJ data
Objpath = "./data/"
Objfile = "Apophis_N2_v2000.obj"
OBJ = readObj(Objpath,Objfile)

volume_total = 0
for i in range(len(OBJ)):
    vertices = OBJ[i].vertices
    faces = OBJ[i].faces
    poly = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                           generation_type='full', vertices=vertices,
                           faces=faces, tol=0., number=None, seed=None,
                           xr=1., yr=1., zr=1.)
    volume_total += poly.contactors[0].volume
    del poly,vertices,faces

# Creat poly list and find COM
Lscale = (volume_total) ** (1 / 3)

# Creat polyhedron body and compute center of mass
COM = np.zeros([3])
mass_total = 0
poly_list = []
for i in range(len(OBJ)):
    vertices = OBJ[i].vertices/Lscale + np.array([-i*1,0,0])
    faces = OBJ[i].faces
    poly = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', vertices=vertices,
                       faces=faces, tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)
    COM += poly.contactors[0].volume * rho * poly.nodes[1].coor
    mass_total += poly.contactors[0].volume * rho
    poly_list.append(poly)
    del poly

# Center of Mass
COM = COM/mass_total
# Euler angle
phi = np.deg2rad(0)
theta = np.deg2rad(15.208605645340988)
psi = np.deg2rad(90)
# longitude and latitude
alpha = np.deg2rad(247)
beta = np.deg2rad(-59)
# Rotation matrix
RM1 = RotationMatrix(phi, 'Z')
RM2 = RotationMatrix(theta, 'X')
RM3 = RotationMatrix(psi, 'Z')
RMb = RotationMatrix(np.pi/2 - beta, 'Y')
RMa = RotationMatrix(alpha, 'Z')
BN = RM3@RM2@RM1@RMb@RMa
# Spin rate
# omegaB_BN = [0.202078195115783E-4,0,0.542630931734589E-4] # rad/sec
omegaB_BN = [0, 0, 2*np.pi/10]
# omegaN_BN = BN.T@omegaB_BN
omega = [0,2*np.pi/10,0]

for i in range(len(poly_list)):
    poly = poly_list[i]
    # Translate to COM
    poly.translate(-COM[0], -COM[1], -COM[2])
    # Rotate to initial attitude, 3-1-3 rotation
    # poly.rotate(description='axis', alpha=alpha, axis=[0., 0., 1.])
    # poly.rotate(description='axis', alpha=np.pi / 2 - beta, axis=[0., 1., 0.])
    # poly.rotate(description='axis', alpha=phi, axis=[0., 0., 1.])
    # poly.rotate(description='axis', alpha=theta, axis=[1., 0., 0.])
    # poly.rotate(description='axis', alpha=psi, axis=[0., 0., 1.])
    # Set initial spin rate
    r = poly.nodes[1].coor.copy()
    velocity = np.cross(omega,r)
    vel = list(np.concatenate((velocity, omega)))
    poly.imposeInitValue(component=[1, 2, 3, 4, 5, 6], value=vel)
    poly.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
    bodies.addAvatar(poly)

LawSPSPx = tact_behav(name='rst01', law='RST_CLB', rstn=1.0, rstt=0.0, fric=0.)
tacts   += LawSPSPx
svSPSPx = see_table(CorpsCandidat='RBDY3', candidat='POLYR',colorCandidat='BLEUx', behav=LawSPSPx,
                     CorpsAntagoniste='RBDY3', antagoniste='POLYR', colorAntagoniste='BLEUx', alert=0.)
svs += svSPSPx

writeBodies(bodies, chemin='DATBOX/')
writeBulkBehav(mat, chemin='DATBOX/', dim=dim , gravy=[0.,0.,0.])
writeTactBehav(tacts, svs, chemin='DATBOX/')
writeDrvDof(bodies, chemin='DATBOX/')
writeDofIni(bodies, chemin='DATBOX/')
writeVlocRlocIni(chemin='DATBOX/')

try:
  visuAvatars(bodies,with_axis=True)
except:
  pass

