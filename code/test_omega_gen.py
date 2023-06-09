# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: gen_sample.py
@date: 1/10/23 13:38
@desc:
"""

import sys

import numpy as np

sys.path.append('./lib')
from pylmgc90.pre import *
from lib.Obj_calss import Obj_class

if not os.path.isdir('./DATBOX'):
    os.mkdir('./DATBOX')
if 'norand' in sys.argv:
    seed = list(range(13))
else:
    seed = None


dim = 3
rho = 5164
Omega = 3e-5
bodies = avatars()
mat    = materials()
svs    = see_tables()
tacts  = tact_behavs()
stone = material(name='STONE', materialType='RIGID', density=rho)
mat.addMaterial(stone)
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

displacement = np.array([1,0,0])

vertices1 = np.array([[2,1,-1],[2,-1,-1],[2,-1,1],[2,1,1],
                      [0,1,-1],[0,-1,-1],[0,-1,1],[0,1,1]])
faces1 = np.array([[1,2,3],[1,4,3],[1,4,8],[1,5,8],[1,2,6],[1,5,6],
                   [7,8,4],[7,3,4],[7,3,2],[7,6,2],[7,8,5],[7,6,5]])
vertices2 = np.array([[0,2,-2],[-3,2,-2],[-3,2,2],[0,2,2],
                      [0,-2,-2],[-3,-2,-2],[-3,-2,2],[0,-2,2]])
faces2 = np.array([[1,2,3],[1,4,3],[1,4,8],[1,5,8],[1,2,6],[1,5,6],
                   [7,8,4],[7,3,4],[7,3,2],[7,2,6],[7,8,5],[7,6,5]])

OBJ1 = Obj_class(vertices = vertices1, faces = faces1)
OBJ1.objTranslation(dS = 1,direction='X')
ver = OBJ1.objRotation(rad = 2*np.pi, axis='X')
OBJ2 = Obj_class(vertices = vertices2, faces = faces2)
OBJ2.objTranslation(dS = -1,direction='X')
# OBJ2.objRotation(rad = np.pi/2, axis='X')

poly1 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', nb_vertices=4, vertices=OBJ1.vertices,
                       faces=OBJ1.faces, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)
poly2 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', nb_vertices=4, vertices=OBJ2.vertices,
                       faces=OBJ2.faces, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)

vel1 = [0,0,0,0,0,2*np.pi]
vel2 = [0,0,0,0,0,2*np.pi]
poly1.imposeInitValue(component=[1,2,3,4,5,6], value=vel1)
poly2.imposeInitValue(component=[1,2,3,4,5,6], value=vel2)
poly1.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
poly2.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')

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
  visuAvatars(bodies)
except:
  pass

