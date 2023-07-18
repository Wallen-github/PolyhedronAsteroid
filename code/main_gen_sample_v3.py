# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: gen_sample.py
@date: 1/10/23 13:38
@desc: This version provides several avatar configuration to test how to set the spin rate
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

# Read OBJ data
Objpath = "./data/"
Objfile = "Apophis_N2_v2000.obj"
OBJ = readObj(Objpath,Objfile)

dim = 3
rho = 10.
Omega = 3e-5
bodies = avatars()
mat    = materials()
svs    = see_tables()
tacts  = tact_behavs()
stone = material(name='STONE', materialType='RIGID', density=rho)
mat.addMaterial(stone)
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

vertices1 = OBJ[0].vertices + np.array([0.8,0,0])
faces1 = OBJ[0].faces
vertices2 = OBJ[1].vertices
faces2 = OBJ[1].faces
poly01 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', vertices=vertices1,
                       faces=faces1, tol=0., number=None, seed=None)
poly02 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', vertices=vertices2,
                       faces=faces2, tol=0., number=None, seed=None)

# vertices1 = OBJ[0].vertices
# faces1 = OBJ[0].faces
# vertices2 = OBJ[1].vertices
# faces2 = OBJ[1].faces

# vertices1 = np.array([[0.5,-0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,-0.5],[0.5,-0.5,-0.5],
#                       [0.,0.,2.],[0.,2.,0.],[0.,0.,-2.],[0.,-2.,0.],[5.,0.6,0.7]]) + np.array([5,0,0])
# vertices2 = np.array([[0.,-2.,2.],[0.,2.,2.],[0.,2.,-2.],[0.,-2.,-2.],
#                       [-1.,-1.,1.],[-1.,1.,1.],[-1.,1.,-1.],[-1.,-1.,-1.],
#                       [0,0,-10]])
# poly01 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
#                        generation_type='vertices', vertices=vertices1)
# poly02 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
#                        generation_type='vertices', vertices=vertices2)

# vertices1 = np.array([[0.5,-0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,-0.5],[0.5,-0.5,-0.5],
#                       [0.,0.,2.],[0.,2.,0.],[0.,0.,-2.],[0.,-2.,0.]]) + np.array([5,0,0])
# vertices2 = np.array([[0.,-2.,2.],[0.,2.,2.],[0.,2.,-2.],[0.,-2.,-2.],
#                       [-1.,-1.,1.],[-1.,1.,1.],[-1.,1.,-1.],[-1.,-1.,-1.]])
# poly01 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
#                        generation_type='vertices', vertices=vertices1)
# poly02 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
#                        generation_type='vertices', vertices=vertices2)

volumetotal = poly01.contactors[0].volume+poly02.contactors[0].volume
GG = 6.6742e-11
UnitL = volumetotal ** (1 / 3)

vertices1 = vertices1/UnitL
vertices2 = vertices2/UnitL
# poly1 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
#                        generation_type='full', vertices=vertices1,faces=faces1, tol=0., number=None, seed=None)
# poly2 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
#                        generation_type='full',  vertices=vertices2,faces=faces2, tol=0., number=None, seed=None)

poly1 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='vertices', vertices=vertices1)
poly2 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='vertices',  vertices=vertices2)
COM = (poly1.contactors[0].volume*rho*poly1.nodes[1].coor
       +poly2.contactors[0].volume*rho*poly2.nodes[1].coor)/(
        poly1.contactors[0].volume*rho+poly2.contactors[0].volume*rho)
# poly1.translate(-COM[0],-COM[1],-COM[2])
# poly2.translate(-COM[0],-COM[1],-COM[2])
# poly1.translate(0.2/UnitL,0,0.)
# poly2.translate(-0.5,0.,0.)

# poly1.rotate(description='axis',alpha=np.pi/4.,axis=[0.,0.,1.])
# poly2.rotate(description='axis',alpha=np.pi/4.,axis=[0.,0.,1.])
# poly2.rotate(description='axis',alpha=np.pi/2.,axis=[1.,0.,0.])
# poly1.rotate(description='axis',alpha=np.pi/4.,axis=[0.,0.,1.])
# poly2.rotate(description='axis',alpha=np.pi/4.,axis=[0.,0.,1.])

omegaN = np.array([0.,0.,2.*np.pi/15])
omegaB1 = np.dot(poly1.bulks[0].axis.T,omegaN)
omegaB2 = np.dot(poly2.bulks[0].axis.T,omegaN)
r1 = poly1.nodes[1].coor
r2 = poly2.nodes[1].coor
vel1 = np.cross(omegaN,r1)
vel2 = np.cross(omegaN,r2)
velo1 = [vel1[0],vel1[1],vel1[2],omegaB1[0],omegaB1[1],omegaB1[2]]
velo2 = [vel2[0],vel2[1],vel2[2],omegaB2[0],omegaB2[1],omegaB2[2]]
poly1.imposeInitValue(component=[1,2,3,4,5,6], value=velo1)
poly2.imposeInitValue(component=[1,2,3,4,5,6], value=velo2)
# poly1.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
# poly2.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')

bodies.addAvatar(poly1)
bodies.addAvatar(poly2)

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

