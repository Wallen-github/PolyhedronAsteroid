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
from pylmgc90.pre import *
from lib.readObj import readObj

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

Objpath = "./data/"
Objfile = "Apophis_N4_v500.obj"
OBJ = readObj(Objpath,Objfile)
for i in range(len(OBJ)):
    vertices = OBJ[i].vertices
    faces = OBJ[i].faces
    poly = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', vertices=vertices,
                       faces=faces, tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)
    poly.imposeInitValue(component=[i+1], value=10)
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
  visuAvatars(bodies)
except:
  pass

