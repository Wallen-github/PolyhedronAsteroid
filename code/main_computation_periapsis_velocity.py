# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: computation_periapsis_velocity.py
@date: 9/10/23 00:03
@desc: 
"""
from __future__ import print_function
import os, sys
import numpy as np
from lib.ExternalTorques import *
import time
import pandas as pd

sys.path.append("..")
from lib.SubFunc import *
from pylmgc90 import chipy

def Initial_computation():
    chipy.Initialize()

    chipy.checkDirectories()

    chipy.utilities_DisableLogMes()
    # timer gravity
    # chipy.timer_InitializeTimers()
    timer_id = chipy.timer_GetNewTimer('gravity computation')
    return timer_id

def Compute_periapsis_velocity(pos_p, vel_inf, max_radius, timer_id):

    chipy.SetDimension(3)

    UnitM, UnitL, UnitT, UnitG = Set_Unit()
    total_time = 3 * 86400 / UnitT

    dt = 1.e-3
    theta = 0.5
    nb_steps = int(total_time / dt)

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
    PosVecE0 = InitialPosVel_Earth_v2(massA, vel_inf, pos_p, total_time)

    coor = np.zeros([nbR3, 6], dtype=np.float32)
    vbeg = np.zeros([nbR3, 6], dtype=np.float32)
    vbeg0 = np.zeros([nbR3, 6], dtype=np.float32)
    fext = np.zeros([nbR3, 6], dtype=np.float32)
    relative_distance = np.zeros(nbR3, dtype=np.float32)
    Record_info = np.zeros([nb_steps, 5], dtype=np.float32)

    chipy.OpenDisplayFiles()
    chipy.WriteDisplayFiles(1)

    # time_start = time.time()

    for k in range(1, nb_steps + 1):
        if k % int(nb_steps / 10) == 0: print(k, '/', (nb_steps + 1))

        chipy.IncrementStep()

        chipy.ComputeFext()

        if k==2:
            for i in range(0, nbR3, 1):
                vbeg0[i, :] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)

        # Get the bodies position and velocity
        for i in range(0, nbR3, 1):
            coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
            vbeg[i, :] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)

        chipy.timer_StartTimer(timer_id)

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

        chipy.timer_StopTimer(timer_id)

        for i in range(0, nbR3, 1):
            chipy.RBDY3_PutBodyVector('Fext_', i + 1, fext[i, :])

        max_mass_index = np.argmax(mass)
        if k == 1:
            dis0 = Get_RelativeDistance(nbR3, mass, coor)
        relative_distance = abs(Get_RelativeDistance(nbR3, mass, coor) - dis0)
        info, loss_num, nomove_num, shifting_num, circular_num = Detect_State(nbR3, max_mass_index, relative_distance,
                                                                              max_radius)
        Record_info[k - 1, :] = np.array([info, loss_num, nomove_num, shifting_num, circular_num])

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

        if info == 0:
            print('System breakup')
            break

    chipy.WriteLastDof()
    chipy.CloseDisplayFiles()
    # chipy.timer_ClearAll()

    chipy.Finalize()

    # time_end = time.time()
    # print('!-------------------------------------!')
    # print('     Computation Cost: ', time_end - time_start, ' second')

    return Record_info, np.linalg.norm(vbeg[max_mass_index,3:6])/np.linalg.norm(vbeg0[max_mass_index,3:6])



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # max_radius = np.array([0.6077247, 0.7198549])
    max_radius = np.array([0.656534, 0.73158497, 0.48802143])
    # pos_p = 1.1 * 6378.14E3  # m # 38012 km for Apophis
    # vel_inf = 8E3  # 5.946E3 # m/s

    timer_id  = Initial_computation()

    RowNum = 10
    ColNum = 10
    periapsis_max = 2.5
    periapsis_min = 1.1
    periapsis_step = (periapsis_max - periapsis_min) / (RowNum-1)
    vel_max = 30E3
    vel_min = 5E3
    vel_step = (vel_max - vel_min) / (ColNum-1)
    Rec_info = np.zeros([RowNum, ColNum], dtype=np.float32)
    Rec_period = np.zeros([RowNum, ColNum], dtype=np.float32)
    Rec_pos = np.zeros([RowNum, ColNum], dtype=np.float32)
    Rec_vel = np.zeros([RowNum, ColNum], dtype=np.float32)

    time_start = time.time()
    for i in range(RowNum):
        for j in range(ColNum):
            print('-------------------')
            print('i = ', i , 'j = ', j)
            print('-------------------')
            pos_p = (periapsis_min + periapsis_step*i) * 6378.14E3  # m # 38012 km for Apophis
            vel_inf = vel_min + vel_step*j  # 5.946E3 # m/s
            # [info, loss_num, nomove_num, shifting_num, circular_num]
            info,wbeg = Compute_periapsis_velocity(pos_p, vel_inf, max_radius, timer_id)

            print('info = ', info[-1,:], '(info, loss_num, nomove_num, shifting_num, circular_num)')

            Rec_pos[i,j] = (periapsis_min + periapsis_step*i) # R_E
            Rec_vel[i,j] = vel_inf/1E3 # km/s
            Rec_period[i,j] = 1/wbeg
            if info[-1, 0] ==0:
                Rec_info[i,j] = -1 # breakup
            elif info[-1, 2] != 0:
                Rec_info[i,j] = 0 # no move
            elif info[-1,3] != 0:
                Rec_info[i,j] = 1 # shifting
            elif info[-1, 4] != 0:
                Rec_info[i,j] = 2 # circular
            elif info[-1, 1] !=0:
                Rec_info[i,j] = 3 # loss some bodies

    time_end = time.time()
    print('!-------------------------------------!')
    print('     Computation Cost: ', time_end - time_start, ' second')

    # Write Data
    out_info = pd.DataFrame(Rec_info)
    out_period = pd.DataFrame(Rec_period)
    out_pos = pd.DataFrame(Rec_pos)
    out_vel = pd.DataFrame(Rec_vel)

    out_info.to_csv('./data/pervel3b_info.txt', sep=',', index=False, header=None)
    out_period.to_csv('./data/pervel3b_period.txt', sep=',', index=False, header=None)
    out_pos.to_csv('./data/pervel3b_periapsis.txt', sep=',', index=False, header=None)
    out_vel.to_csv('./data/pervel3b_velinf.txt', sep=',', index=False, header=None)
