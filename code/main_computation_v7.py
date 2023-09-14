# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: main_computation_v7.py
@date: 9/10/23 00:10
@desc: 
"""
from __future__ import print_function
import os, sys
import numpy as np
import math as mt
from lib.ExternalTorques import *
from lib.PlotFunc import *
import time

sys.path.append("..")
from lib.SubFunc import *
from pylmgc90 import chipy
# from pykdgrav import *

chipy.Initialize()

chipy.checkDirectories()

chipy.utilities_DisableLogMes()
# timer gravity
timer_id = chipy.timer_GetNewTimer('gravity computation')

chipy.SetDimension(3)

UnitM, UnitL, UnitT, UnitG  = Set_Unit()
total_time = 3*86400/UnitT

dt = 1.e-3
theta = 0.5
nb_steps = int(total_time/dt)

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
pos_p = 1.1 * 6378.14E3 #m # 38012 km for Apophis
vel_inf = 7.9E3 # 5.946E3 # m/s
PosVecE0 = InitialPosVel_Earth_v2(massA, vel_inf, pos_p, total_time)

coor = np.zeros([nbR3, 6], dtype=np.float32)
vbeg = np.zeros([nbR3, 6], dtype=np.float32)
fext = np.zeros([nbR3, 6], dtype=np.float32)
fext_pykd = np.zeros([nbR3, 6], dtype=np.float32)
NB_pre = np.zeros([3, 3], dtype=np.float32)

# varibles for plot
Plot_T = np.zeros([nb_steps,1],dtype=np.float32)
Plot_Momentum = np.zeros([nb_steps,4],dtype=np.float32)
Plot_Kinetic = np.zeros([nb_steps,1],dtype=np.float32)
Plot_Energy = np.zeros([nb_steps,1],dtype=np.float32)
Plot_Angle = np.zeros([nb_steps,nbR3],dtype=np.float32)
Plot_Distance = np.zeros([nb_steps,1],dtype=np.float32)
Plot_vector = np.zeros([nb_steps,3],dtype=np.float32)
Plot_Vbeg = np.zeros([nb_steps, nbR3*6], dtype=np.float32)
Plot_inertia = np.zeros([nb_steps, 3], dtype=np.float32)
Plot_EA = np.zeros([nb_steps, 3], dtype=np.float32)
Plot_PosVecE = np.zeros([nb_steps, 6], dtype=np.float32)

chipy.OpenDisplayFiles()
chipy.WriteDisplayFiles(1)

time_start = time.time()

for k in range(1, nb_steps + 1):
    if k % int(nb_steps/10) == 0: print(k, '/', (nb_steps + 1))

    chipy.IncrementStep()

    chipy.ComputeFext()

    # Get the bodies position and velocity
    for i in range(0, nbR3, 1):
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        vbeg[i, :] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)

    # Get Energy, momentum and other tracking parameters
    Energy, Kinetic, momentum = Get_EnergyMomentum(nbR3, GG=1)
    # I_global, I_body, DCM_NB = Get_TotalMomentOfInertia_com(nbR3)
    I_global, I_body, DCM_NB = Get_Get_TotalMomentOfInertia_v2(nbR3)
    theta1, theta2, theta3 = DCM2EA_313(DCM_NB)

    chipy.timer_StartTimer(timer_id)

    # Get Earth acceleration and position
    posE, PosVec = EarthPos(massA, k, dt, PosVecE0)
    PosVecE0 = PosVec
    CMP, coor_cm = Get_CenterMass(coor[:, 0:3], mass)

    # GG = UnitG
    # massE = 5.974227245203837E24 / UnitM  # kg, Earth mass
    # fext_pykd[:, 0:3] = Accel(coor[:, 0:3], mass, G=UnitG)
    # for i in range(0, nbR3, 1):
    #     r0_norm = np.linalg.norm(posE - coor[i, 0:3])
    #     posE_norm = np.linalg.norm(posE)
    #     fext_pykd[i,0:3] += GG*massE*(posE-coor[i,0:3])/r0_norm**3 - GG*massE*posE/posE_norm**3
    #     fext_pykd[i, :] = fext_pykd[i, :] * mass[i]

    #Get N-body force
    GG = UnitG
    massE = 5.974227245203837E24 / UnitM  # kg, Earth mass
    for i in range(0, nbR3, 1):
        F_nbody = np.zeros([6], dtype=np.float32)
        for j in range(0, nbR3, 1):
            if j != i:
                r_norm = np.linalg.norm(coor[j,0:3]-coor[i,0:3])
                F_nbody[0:3] += GG*mass[j]*(coor[j,0:3]-coor[i,0:3])/r_norm**3
        r0_norm = np.linalg.norm(posE-coor[i,0:3])
        posE_norm = np.linalg.norm(posE)
        # F_nbody[0:3] += GG*massE*(posE-coor[i,0:3])/r0_norm**3 - GG*massE*posE/posE_norm**3
        fext[i,0:3] = mass[i] * F_nbody[0:3]

    chipy.timer_StopTimer(timer_id)

    for i in range(0, nbR3, 1):
        chipy.RBDY3_PutBodyVector('Fext_', i + 1, fext[i, :])

    # Record data for plot
    theta_list, NB_pre, vector_bodyfix = Get_RelativeAngel(nbR3, mass, coor, NB_pre)

    Plot_T[k-1,0] = dt*(k-1)*UnitT/3600
    Plot_Energy[k-1,0] = Energy
    Plot_Kinetic[k-1,0] = Kinetic
    Plot_Angle[k-1,:] = np.rad2deg(theta_list)
    Plot_Distance[k-1, :] = np.linalg.norm(coor[0,0:3] - coor[1,0:3])
    Plot_vector[k - 1, :] = vector_bodyfix
    Plot_Momentum[k-1,:] = momentum
    Plot_inertia[k-1,:] = I_body
    Plot_EA[k-1,:] = [theta1, theta2, theta3]
    Plot_PosVecE[k - 1, :] = np.append(PosVec[0:3]*UnitL/1E3,PosVec[3:6]*UnitL/1E3/UnitT)
    for i in range(0, nbR3, 1):
        Plot_Vbeg[k-1,i*6:i*6+6] = vbeg[i, :]

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
print('     Computation Cost: ',time_end - time_start, ' second')

# Plot
Plot_T = np.delete(Plot_T, 0, axis=0)
Plot_Momentum = np.delete(Plot_Momentum, 0, axis=0)
Plot_Energy = np.delete(Plot_Energy, 0, axis=0)
Plot_Kinetic = np.delete(Plot_Kinetic, 0, axis=0)
Plot_Angle = np.delete(Plot_Angle, 0, axis=0)
Plot_Distance = np.delete(Plot_Distance, 0, axis=0)
Plot_vector = np.delete(Plot_vector, 0, axis=0)
Plot_Vbeg = np.delete(Plot_Vbeg, 0, axis=0)
Plot_EA = np.delete(Plot_EA, 0, axis=0)
Plot_inertia = np.delete(Plot_inertia, 0, axis=0)

plot_momentum(Plot_T,Plot_Momentum)
# plot_kinetic(Plot_T,Plot_Kinetic)
# plot_energy(Plot_T,Plot_Energy)
# plot_omega(Plot_T, Plot_Vbeg, nbR3, 2, 1)
# plot_velocity(Plot_T, Plot_Vbeg, nbR3, 2, 1)
# plot_eulerangle(Plot_T, Plot_EA)
plot_inertia(Plot_T, Plot_inertia)
plot_angel(Plot_T, Plot_Angle, Plot_vector)
plot_distance(Plot_T, Plot_Distance)
# plot_earthtraj(Plot_PosVecE)
