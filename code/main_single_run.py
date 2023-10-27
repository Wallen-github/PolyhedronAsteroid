# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: main_single_run.py
@date: 9/19/23 11:19
@desc: 
"""
from __future__ import print_function
import os, sys
import traceback
import requests
from lib.CustomError import CustomError
sys.path.append('./lib')
import numpy as np
from pylmgc90.pre import *
from lib.readObj import readObj
from lib.SubFunc import *
from lib.RotationMatrix import *
if not os.path.isdir('./DATBOX'):
    os.mkdir('./DATBOX')
if 'norand' in sys.argv:
    seed = list(range(13))
else:
    seed = None
import math as mt
from lib.ExternalTorques import *
from lib.PlotFunc import *
import time

sys.path.append("..")
from lib.SubFunc import *
from pylmgc90 import chipy
# from pykdgrav import *


def gen_sample(Objpath,Objfile, alpha=0., beta=0., ctacSurf=0.005, shapechange=False, fric=0., rstn=1.0, rstt=0.0):
    # Read OBJ data
    OBJ = readObj(Objpath, Objfile)

    # change surface
    if shapechange==True:
        OBJ[0].vertices[:,0] = OBJ[0].vertices[:,0] + 0.05
        OBJ[1].vertices[:,0] = OBJ[1].vertices[:,0] - 0.05
        contact_surface = ctacSurf # 0.005 #0.004330127018922
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
    # I_global, I_diag, DCM = Get_TotalMomentOfInertia_gen(poly0_list)
    I_global, I_diag, DCM = Get_TotalMomentOfInertia_gen_v2(poly0_list)
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
        poly.translate(-COM[0], -COM[1], -COM[2])
        poly.rotate(description='Euler', phi=theta1, theta=theta2, psi=theta3)
        poly.rotate(description='axis', alpha=np.pi, axis=[0., 0., 1.])

        # User defined orientation
        poly.rotate(description='axis', alpha=alpha, axis=[0., 1., 0.])
        poly.rotate(description='axis', alpha=beta, axis=[0., 0., 1.])
        RM1 = RotationMatrix(alpha, 'Y')
        RM2 = RotationMatrix(beta, 'Z')
        RM = np.dot(RM2, RM1)
        omegaU = np.dot(RM, omegaB)
        omegaBi = np.dot(poly.bulks[0].axis.T, omegaU)

        # omegaBi = np.dot(poly.bulks[0].axis.T, omegaB)
        r = poly.nodes[1].coor
        vel = np.cross(omegaB, r)
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

    print('max vertex radius:', max_radius)
    return bodies,max_radius

def computation(max_radius, pos_p = 0., vel_inf= 0., trajType='general'):
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
    if trajType == 'apophis':
        PosVecE0 = InitialPosVel_Earth(massA, total_time)
    elif trajType == 'general':
        if pos_p == 0:
            msg = 'Plase input the perigee value'
            try:
                raise CustomError(msg)
            except CustomError as e:
                traceback.print_exc()
                sys.exit(-1)
        # pos_p = 1.9 * 6378.14E3 #m # 38012 km for Apophis
        # vel_inf = 5E3 # 5.946E3 # m/s
        PosVecE0 = InitialPosVel_Earth_v2(massA, vel_inf, pos_p, total_time)


    coor = np.zeros([nbR3, 6], dtype=np.float32)
    vbeg = np.zeros([nbR3, 6], dtype=np.float32)
    fext = np.zeros([nbR3, 6], dtype=np.float32)
    fint = np.zeros([nbR3, 6], dtype=np.float32)

    # varibles for plot
    Plot_T = np.zeros([nb_steps, 1], dtype=np.float32)
    Plot_Momentum = np.zeros([nb_steps, 4], dtype=np.float32)
    Plot_Kinetic = np.zeros([nb_steps, 1], dtype=np.float32)
    Plot_Energy = np.zeros([nb_steps, 1], dtype=np.float32)
    Plot_Effective = np.zeros([nb_steps, 2], dtype=np.float32)
    Plot_Distance = np.zeros([nb_steps, nbR3], dtype=np.float32)
    Plot_Vbeg = np.zeros([nb_steps, nbR3 * 6], dtype=np.float32)
    Plot_inertia = np.zeros([nb_steps, 3], dtype=np.float32)
    Plot_EA = np.zeros([nb_steps, 3], dtype=np.float32)
    Plot_PosVecE = np.zeros([nb_steps, 6], dtype=np.float32)
    Record_info = np.zeros([nb_steps, 4], dtype=np.float32)

    chipy.OpenDisplayFiles()
    chipy.WriteDisplayFiles(1)

    time_start = time.time()

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
        # I_global, I_body, DCM_NB = Get_TotalMomentOfInertia_v2(nbR3)
        theta1, theta2, theta3 = DCM2EA_313(DCM_NB)

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
        Record_info[k - 1, :] = np.array([info, I_d, w_e, np.sqrt(Momentum[3])])
        if info == 0:
            print('System breakup')
            Record_info[- 1, :] = Record_info[k - 1, :]
            break
        Plot_Distance[k - 1, :] = relative_distance
        Plot_T[k - 1, 0] = dt * (k - 1) * UnitT / 3600
        Plot_Energy[k - 1, 0] = Energy
        Plot_Kinetic[k - 1, 0] = Kinetic
        Plot_Momentum[k - 1, :] = Momentum
        Plot_inertia[k - 1, :] = I_body
        Plot_Effective[k - 1, :] = [np.linalg.det(I_global), w_e]
        Plot_EA[k - 1, :] = [theta1, theta2, theta3]
        Plot_PosVecE[k - 1, :] = np.append(PosVec[0:3] * UnitL / 1E3, PosVec[3:6] * UnitL / 1E3 / UnitT)
        for i in range(0, nbR3, 1):
            Plot_Vbeg[k - 1, i * 6:i * 6 + 6] = vbeg[i, :]

        chipy.ComputeBulk()
        chipy.ComputeFreeVelocity()

        chipy.SelectProxTactors()
        chipy.RecupRloc()

        chipy.ExSolver(stype, norm, tol, relax, gs_it1, gs_it2)

        chipy.StockRloc()

        chipy.ComputeDof()
        chipy.UpdateStep()

        chipy.WriteDisplayFiles(freq_display, ref_radius)

        chipy.WriteOutDof(freq_write)
        chipy.WriteOutVlocRloc(freq_write)

        chipy.overall_CleanWriteOutFlags()

    chipy.WriteLastDof()
    chipy.CloseDisplayFiles()

    chipy.Finalize()

    time_end = time.time()
    print('!-------------------------------------!')
    print('     Computation Cost: ', time_end - time_start, ' second')
    print('Final state: ', Record_info[-1, :], '(loss_num, nomove_num, shifting_num, circular_num)')

    Plot_T = np.delete(Plot_T, 0, axis=0)
    Plot_Momentum = np.delete(Plot_Momentum, 0, axis=0)
    Plot_Energy = np.delete(Plot_Energy, 0, axis=0)
    Plot_Effective = np.delete(Plot_Effective, 0, axis=0)
    Plot_Kinetic = np.delete(Plot_Kinetic, 0, axis=0)
    Plot_Distance = np.delete(Plot_Distance, 0, axis=0)
    Plot_Vbeg = np.delete(Plot_Vbeg, 0, axis=0)
    Plot_EA = np.delete(Plot_EA, 0, axis=0)
    Plot_inertia = np.delete(Plot_inertia, 0, axis=0)
    Record_info = np.delete(Record_info, 0, axis=0)

    return Plot_T, Plot_Effective, Plot_Momentum, Plot_Distance,Plot_inertia


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    Objpath = "./data/"
    Objfile = "Apophis_N2_v2000.obj"
    alpha = np.deg2rad(0)
    beta = np.deg2rad(0)
    ctacSurf = 0.005
    bodies, max_radius = gen_sample(Objpath,Objfile, alpha, beta, ctacSurf, shapechange=False, fric=0.5, rstn=0.55, rstt = 0.55)

    vel_max = 20E3
    vel_min = 0E3
    vel_step = (vel_max - vel_min) / (10 - 1)
    vel_inf = 3E3 # vel_min + vel_step * 2  # 5.946E3 # m/s
    pos_p = 1.5 * 6378.14E3  # m # 38012 km for Apophis
    Plot_T, Plot_Effective, Plot_Momentum, Plot_Distance,Plot_inertia = computation(max_radius,pos_p,vel_inf, trajType='general')

    plot_momentum(Plot_T, Plot_Momentum)
    plot_Effective(Plot_T, Plot_Effective)
    plot_distance(Plot_T, Plot_Distance)
    plot_inertia(Plot_T, Plot_inertia)

    # Obtain the file name including suffix
    name = os.path.basename(__file__)
    path = os.getcwd()  # Obtain the current path
    mesg = 'Finish: ds_orbit_2b_test'
    requests.post("https://ntfy.sh/wallen1732_notfy_python",
                  data=mesg.encode(encoding='utf-8'))

    # try:
    #     visuAvatars(bodies, with_axis=True)
    # except:
    #     pass

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
