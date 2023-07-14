# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: main_computation_v4.py
@date: 7/13/23 16:55
@desc: 
"""
from __future__ import print_function
import os, sys
import numpy as np
from lib.PlotFunc import *

sys.path.append("..")
from lib.TraParam import Get_EnergyMomentum
from pylmgc90 import chipy
from pykdgrav import *

chipy.Initialize()

chipy.checkDirectories()

chipy.utilities_DisableLogMes()
# timer gravity
timer_id = chipy.timer_GetNewTimer('gravity computation')

chipy.SetDimension(3)
UnitT = 122405.34764413263
time = 2*np.pi/UnitT
dt = 1.e-3
theta = 0.5
nb_steps = 100000

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

mass = np.zeros(nbR3)
for i in range(nbR3):
    mass[i] = chipy.RBDY3_GetMass(i + 1)

coor = np.empty([nbR3, 6], dtype=float)
p_coor = np.empty([nbR3, 3], dtype=float)
vbeg = np.empty([nbR3, 6], dtype=float)
fext = np.empty([nbR3, 6], dtype=float)

# varibles for plot
Plot_T = np.empty([nb_steps,1],dtype=float)
Plot_Momentum = np.empty([nb_steps,4],dtype=float)
Plot_Kinetic = np.empty([nb_steps,1],dtype=float)
Plot_Energy = np.empty([nb_steps,1],dtype=float)
Plot_Vbeg = np.empty([nb_steps, nbR3*6], dtype=float)

chipy.OpenDisplayFiles()
chipy.WriteDisplayFiles(1)

for k in range(1, nb_steps + 1):
    print(k, '/', (nb_steps + 1))
    # Globinertia = chipy.RBDY3_GetGlobInertia(2)
    # Bodyinertia = chipy.RBDY3_GetBodyInertia(2)
    # print('body 2, Global inertia = ', Globinertia, 'Body inertia = ', Bodyinertia)
    #
    chipy.IncrementStep()

    chipy.ComputeFext()

    fext = np.empty([nbR3, 6], dtype=float)

    for i in range(0, nbR3, 1):
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        vbeg[i, :] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)
        p_coor[i, 0:3] = coor[i, 0:3]
        fext[i, :] = 0.

    chipy.timer_StartTimer(timer_id)
    # fext[:, 0:3] = Accel(p_coor, mass, G=6.6742e-11)
    fext[:, 0:3] = Accel(p_coor, mass, G=1.)
    for i in range(0, nbR3, 1):
        fext[i, :] = fext[i, :] * mass[i]
    chipy.timer_StopTimer(timer_id)

    for i in range(0, nbR3, 1):
        chipy.RBDY3_PutBodyVector('Fext_', i + 1, fext[i, :])

    Energy, Kinetic, momentum = Get_EnergyMomentum(nbR3, GG=1)
    # print('Momentum = ', momentum)

    # Record data for plot
    Plot_T[k-1,0] = dt*(k-1)
    Plot_Energy[k-1,0] = Energy
    Plot_Kinetic[k-1,0] = Kinetic
    Plot_Momentum[k-1,:] = momentum
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

# Plot
Plot_T = np.delete(Plot_T, 0, axis=0)
Plot_Momentum = np.delete(Plot_Momentum, 0, axis=0)
Plot_Energy = np.delete(Plot_Energy, 0, axis=0)
Plot_Kinetic = np.delete(Plot_Kinetic, 0, axis=0)
Plot_Vbeg = np.delete(Plot_Vbeg, 0, axis=0)

plot_momentum(Plot_T,Plot_Momentum)
plot_kinetic(Plot_T,Plot_Kinetic)
plot_energy(Plot_T,Plot_Energy)
plot_omega(Plot_T, Plot_Vbeg, nbR3, 2, 1)
plot_velocity(Plot_T, Plot_Vbeg, nbR3, 2, 1)