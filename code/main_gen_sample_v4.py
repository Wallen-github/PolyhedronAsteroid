# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: main_gen_sample_v4.py
@date: 7/13/23 16:56
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

# Read OBJ data
Objpath = "./data/"
Objfile = "Apophis_N2_v2000.obj"
OBJ = readObj(Objpath,Objfile)

# Set dimension and density
dim = 3
rho = 1.
bodies = avatars()
mat    = materials()
svs    = see_tables()
tacts  = tact_behavs()
stone = material(name='STONE', materialType='RIGID', density=rho)
mat.addMaterial(stone)
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

# Read vertices and faces
volume_total = 0
COM = np.zeros([3])
for i in range(len(OBJ)):
    vertices = OBJ[i].vertices
    faces = OBJ[i].faces
    poly0 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                             generation_type='full', vertices=vertices,faces=faces)
    volume_total += poly0.contactors[0].volume
    COM += poly0.contactors[0].volume * poly0.nodes[1].coor

UnitL = volume_total ** (1 / 3)
COM = COM/volume_total/UnitL

# Set scaled polyhedron
poly_list = []
omegaN = np.array([0.,0.,2.*np.pi/10])
for i in range(len(OBJ)):
    vertices = OBJ[i].vertices/UnitL
    faces = OBJ[i].faces
    poly = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                            generation_type='full', vertices=vertices, faces=faces)
    poly.translate(-COM[0], -COM[1], -COM[2])

    omegaB = np.dot(poly.bulks[0].axis.T, omegaN)
    r = poly.nodes[1].coor
    vel = np.cross(omegaN, r)
    velo = [vel[0], vel[1], vel[2], omegaB[0], omegaB[1], omegaB[2]]
    poly.imposeInitValue(component=[1, 2, 3, 4, 5, 6], value=velo)
    poly_list.append(poly)
    bodies.addAvatar(poly)


LawSPSPx = tact_behav(name='rst01', law='RST_CLB', rstn=1.0, rstt=0.0, fric=0.5)
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