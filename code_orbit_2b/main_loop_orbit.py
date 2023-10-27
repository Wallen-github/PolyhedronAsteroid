# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: main_loop_orbit.py
@date: 9/21/23 14:51
@desc: 
"""


from __future__ import print_function
import sys
sys.path.append('./lib')
import requests
import os
import numpy as np
from pylmgc90.pre import *
from lib.readObj import readObj
from lib.RotationMatrix import *
from lib.ExternalTorques import *
import time
import pandas as pd
from lib.SubFunc import *
from pylmgc90 import chipy
if not os.path.isdir('./DATBOX'):
    os.mkdir('./DATBOX')
if 'norand' in sys.argv:
    seed = list(range(13))
else:
    seed = None
sys.path.append("..")


def gen_sample(Objpath, Objfile, alpha, beta, fric, rstn=1.0, rstt=0.0):
    # Read OBJ data
    OBJ = readObj(Objpath, Objfile)

    # change surface
    OBJ[0].vertices[:,0] = OBJ[0].vertices[:,0] + 0.05
    OBJ[1].vertices[:,0] = OBJ[1].vertices[:,0] - 0.05
    contact_surface = 0.00005 #0.00005 #
    tri_length = np.sqrt(4*contact_surface/np.sqrt(3))
    tri_vertex = np.array([[0.1, 0., np.sqrt(3)/4*tri_length],
                           [0.1, tri_length/2, -np.sqrt(3)/4*tri_length],
                           [0.1, -tri_length/2, -np.sqrt(3)/4*tri_length]])
    OBJ[0].vertices = np.vstack((OBJ[0].vertices, tri_vertex))
    OBJ[1].vertices = np.vstack((OBJ[1].vertices, tri_vertex))

    # Set dimension and density
    dim = 3
    rho = 1.
    bodies = avatars()
    mat = materials()
    svs = see_tables()
    tacts = tact_behavs()
    stone = material(name='STONE', materialType='RIGID', density=rho)
    mat.addMaterial(stone)
    mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

    # Set vertices and faces
    volume_total = 0
    COM = np.zeros([3])
    poly0_list = []
    faces_list = []
    for i in range(len(OBJ)):
        vertices = OBJ[i].vertices
        faces = OBJ[i].faces
        # if i==0:
        #     vertices += np.array([0.5,0,0])
        poly0 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                                generation_type='vertices', vertices=vertices)
        volume_total += poly0.contactors[0].volume
        COM += poly0.contactors[0].volume * poly0.nodes[1].coor
        poly0_list.append(poly0)

    # Get the scale length to normalized total mass to unity
    UnitL = volume_total ** (1 / 3)
    COM = COM / volume_total / UnitL

    # Get Total moment of inertia and DCM
    I_global, I_diag, DCM = Get_TotalMomentOfInertia_gen(poly0_list)
    # I_global, I_diag, DCM = Get_TotalMomentOfInertia_gen_v2(poly0_list)
    theta1, theta2, theta3 = DCM2EA_313(DCM)
    # theta1 = np.pi/4
    # theta2 = 0
    # theta3 = 0

    # Set scaled polyhedron
    poly_list = []
    omegaB = np.array([0.027030895893665, 0., 0.072584774502740])
    # omegaB = np.array([0.,0.,0.02*np.pi])
    max_radius = np.zeros(len(OBJ), dtype=np.float32)

    for i in range(len(OBJ)):
        vertices = OBJ[i].vertices / UnitL
        faces = poly0_list[i].contactors[0].connectivity
        poly = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                               generation_type='full', vertices=vertices, faces=faces)
        # Translate and Rotate to Principal axis frame
        poly.translate(-COM[0], -COM[1], -COM[2])
        poly.rotate(description='Euler', phi=theta1, theta=theta2, psi=theta3)
        poly.rotate(description='axis', alpha=np.pi, axis=[0., 0., 1.])

        # User defined orientation
        poly.rotate(description='axis', alpha=alpha, axis=[0., 1., 0.])
        poly.rotate(description='axis', alpha=beta, axis=[0., 0., 1.])
        RM1 = RotationMatrix(alpha, 'Y')
        RM2 = RotationMatrix(beta, 'Z')
        RM = np.dot(RM2,RM1)
        omegaU = np.dot(RM, omegaB)

        omegaBi = np.dot(poly.bulks[0].axis.T, omegaU)
        r = poly.nodes[1].coor
        vel = np.cross(omegaU, r)
        velo = [vel[0], vel[1], vel[2], omegaBi[0], omegaBi[1], omegaBi[2]]
        poly.imposeInitValue(component=[1, 2, 3, 4, 5, 6], value=velo)
        poly_list.append(poly)
        bodies.addAvatar(poly)
        max_radius[i] = max(np.linalg.norm(poly.contactors[0].vertices - poly.contactors[0].shift, axis=1))

    LawSPSPx = tact_behav(name='rst01', law='RST_CLB', rstn=rstn, rstt=rstt, fric=fric)
    tacts += LawSPSPx
    svSPSPx = see_table(CorpsCandidat='RBDY3', candidat='POLYR', colorCandidat='BLEUx', behav=LawSPSPx,
                        CorpsAntagoniste='RBDY3', antagoniste='POLYR', colorAntagoniste='BLEUx', alert=0.)
    svs += svSPSPx

    writeBodies(bodies, chemin='DATBOX/')
    writeBulkBehav(mat, chemin='DATBOX/', dim=dim, gravy=[0., 0., 0.])
    writeTactBehav(tacts, svs, chemin='DATBOX/')
    writeDrvDof(bodies, chemin='DATBOX/')
    writeDofIni(bodies, chemin='DATBOX/')
    writeVlocRlocIni(chemin='DATBOX/')

    return bodies,max_radius

def computation(pos_p, vel_inf, max_radius):
    chipy.Initialize()

    chipy.checkDirectories()

    chipy.utilities_DisableLogMes()

    chipy.SetDimension(3)

    UnitM, UnitL, UnitT, UnitG = Set_Unit()
    total_time = 3 * 86400 / UnitT

    dt = 1.e-3
    theta = 0.5
    nb_steps = int(total_time / dt)

    echo = 0
    freq_display = 400
    ref_radius = 5.
    freq_write = 400

    # freq_update = 600

    chipy.nlgs_3D_DiagonalResolution()
    chipy.RBDY3_NewRotationScheme()

    chipy.PRPRx_UseCpCundallDetection(300)
    chipy.PRPRx_LowSizeArrayPolyr(20)

    # type = 'Stored_Delassus_Loops         '
    stype = 'Exchange_Local_Global         '
    norm = 'QM/16'
    tol = 0.1e-3
    relax = 1.0
    gs_it1 = 10
    gs_it2 = 200

    chipy.utilities_logMes('INIT TIME STEPPING')
    chipy.TimeEvolution_SetTimeStep(dt)
    chipy.Integrator_InitTheta(theta)

    chipy.utilities_logMes('READ BODIES')
    chipy.ReadBodies()
    chipy.LoadTactors()

    chipy.utilities_logMes('READ INI DOF')
    chipy.ReadIniDof()

    chipy.utilities_logMes('READ BEHAVIOURS')
    chipy.ReadBehaviours()
    chipy.LoadBehaviours()

    chipy.utilities_logMes('READ INI Vloc Rloc')
    chipy.ReadIniVlocRloc()

    chipy.utilities_logMes('READ DRIVEN DOF')
    chipy.ReadDrivenDof()

    chipy.WriteBodies()
    #
    chipy.WriteBehaviours()

    chipy.WriteDrivenDof()

    chipy.ComputeMass()

    # chipy.RBDY3_FatalDamping()

    nbR3 = chipy.RBDY3_GetNbRBDY3()

    mass = np.zeros(nbR3, dtype=np.float32)
    for i in range(nbR3):
        mass[i] = chipy.RBDY3_GetMass(i + 1)

    # Get the initial position of Earth
    massA = np.sum(mass)
    # PosVecE0 = InitialPosVel_Earth(massA, total_time)
    PosVecE0 = InitialPosVel_Earth_v2(massA, vel_inf, pos_p, total_time)

    coor = np.zeros([nbR3, 6], dtype=np.float32)
    vbeg = np.zeros([nbR3, 6], dtype=np.float32)
    fext = np.zeros([nbR3, 6], dtype=np.float32)
    fint = np.zeros([nbR3, 6], dtype=np.float32)

    Record_info = np.zeros([nb_steps, 4], dtype=np.float32)

    chipy.OpenDisplayFiles()
    chipy.WriteDisplayFiles(1)

    for k in range(1, nb_steps + 1):
        if k % int(nb_steps / 10) == 0: print(k, '/', (nb_steps + 1))

        chipy.IncrementStep()

        chipy.ComputeFext()

        # Get the bodies position and velocity
        for i in range(0, nbR3, 1):
            coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
            vbeg[i, :] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)
            fint[i, :] = chipy.RBDY3_GetBodyVector('Fint_', i + 1)

        # Get Energy, momentum and other tracking parameters
        Energy, Kinetic, Momentum = Get_EnergyMomentum(nbR3, GG=1)
        if k == 1:
            I_d = 0
            w_e = 0
        else:
            I_d, w_e = Get_EffectiveInertiaSpin(Kinetic, Momentum)
        I_global, I_body, DCM_NB = Get_TotalMomentOfInertia_com(nbR3)

        # Get Earth acceleration and position
        posE, PosVec = EarthPos(massA, k, dt, PosVecE0)
        PosVecE0 = PosVec
        # Get N-body force
        GG = UnitG
        massE = 5.974227245203837E24 / UnitM  # kg, Earth mass
        for i in range(0, nbR3, 1):
            F_nbody = np.zeros([6], dtype=np.float32)
            for j in range(0, nbR3, 1):
                if j != i:
                    r_norm = np.linalg.norm(coor[j, 0:3] - coor[i, 0:3])
                    F_nbody[0:3] += GG * mass[j] * (coor[j, 0:3] - coor[i, 0:3]) / r_norm ** 3
            r0_norm = np.linalg.norm(posE - coor[i, 0:3])
            posE_norm = np.linalg.norm(posE)
            F_nbody[0:3] += GG * massE * (posE - coor[i, 0:3]) / r0_norm ** 3 - GG * massE * posE / posE_norm ** 3
            fext[i, 0:3] = mass[i] * F_nbody[0:3]

        for i in range(0, nbR3, 1):
            chipy.RBDY3_PutBodyVector('Fext_', i + 1, fext[i, :])

        max_mass_index = np.argmax(mass)
        if k == 1:
            dis0 = Get_RelativeDistance(nbR3, mass, coor)
        relative_distance = abs(Get_RelativeDistance(nbR3, mass, coor) - dis0)
        info, loss_num, nomove_num, shifting_num, circular_num = Detect_State(nbR3, max_mass_index, relative_distance,
                                                                              max_radius)
        Record_info[k - 1, :] = np.array([info, np.sum(relative_distance), w_e, np.sqrt(Momentum[3])])
        if info == 0:
            print('System breakup')
            Record_info[- 1, :] = Record_info[k - 1, :]
            break

        chipy.ComputeBulk()
        chipy.ComputeFreeVelocity()

        chipy.SelectProxTactors()
        chipy.RecupRloc()

        chipy.ExSolver(stype, norm, tol, relax, gs_it1, gs_it2)

        chipy.StockRloc()

        chipy.ComputeDof()
        chipy.UpdateStep()

        # chipy.WriteDisplayFiles(freq_display, ref_radius)

        chipy.WriteOutDof(freq_write)
        chipy.WriteOutVlocRloc(freq_write)

        chipy.overall_CleanWriteOutFlags()

    chipy.WriteLastDof()
    chipy.CloseDisplayFiles()

    chipy.Finalize()

    Record_info = np.delete(Record_info, 0, axis=0)

    return Record_info

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    time_start = time.time()

    Objpath = "./data/"
    Objfile = "Apophis_N2_v2000.obj"
    bodies, max_radius = gen_sample(Objpath, Objfile, alpha=0., beta=0., fric=0.5, rstn=0.55, rstt=0.55)

    RowNum = 10
    ColNum = 10
    periapsis_max = 6.0
    periapsis_min = 1.1
    periapsis_step = (periapsis_max - periapsis_min) / (RowNum - 1)
    vel_max = 20E3
    vel_min = 0E3
    vel_step = (vel_max - vel_min) / (ColNum - 1)
    Rec_info = np.zeros([RowNum * ColNum, 6], dtype=np.float32)

    time_start = time.time()
    for i in range(RowNum):
        for j in range(ColNum):
            print('-------------------')
            print('i = ', i, 'j = ', j)
            print('-------------------')
            pos_p = (periapsis_min + periapsis_step * i) * 6378.14E3  # m # 38012 km for Apophis
            vel_inf = vel_min + vel_step * j  # 5.946E3 # m/s

            info = computation(pos_p, vel_inf, max_radius)

            Rec_info[j + i * ColNum, 0] = (periapsis_min + periapsis_step * i) # RE
            Rec_info[j + i * ColNum, 1] = vel_inf/1E3  # km
            Rec_info[j + i * ColNum, 2] = info[-1, 1]  # ds - ds0
            Rec_info[j + i * ColNum, 3] = info[0, 2] / info[-1, 2]  # we0/we
            Rec_info[j + i * ColNum, 4] = info[-1, 3] / info[0, 3]  # H/H0
            Rec_info[j + i * ColNum, 5] = info[-1, 0]  # breakup info
            print('(perigee, vel_inf, ds - ds0, we0/we, H/H0, info)')
            print('Rec_info = ', Rec_info[j + i * ColNum, :])

    time_end = time.time()
    print('!-------------------------------------!')
    print('     Computation Cost: ', time_end - time_start, ' second')

    # Write Data
    out_info = pd.DataFrame(Rec_info)
    # Row index: 1. alpha, 2. beta, 3. ds - ds0, 4. we0/we, 5. H/H0, 6. info
    out_info.to_csv('./data/ds_orbit_2b_v2.txt', sep=',', index=False, header=False)

    # 获取文件名（含后缀）
    name = os.path.basename(__file__)
    path = os.getcwd()  # 获取当前工作目录路径
    mesg = 'Finish: ds_orbit_2b_v2'
    requests.post("https://ntfy.sh/wallen1732_notfy_python",
                  data=mesg.encode(encoding='utf-8'))

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
