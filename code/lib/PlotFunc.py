# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: PlotFunc.py
@date: 7/6/23 13:28
@desc: 
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_momentum(Plot_T,Plot_Momentum):
    # plt.figure()
    # plt.plot(Plot_T,Plot_Momentum[:,3])
    # plt.xlabel('time')
    # plt.ylabel('Momentum Square')
    # plt.show()

    plt.figure()
    plt.plot(Plot_T, (Plot_Momentum[:, 3] - Plot_Momentum[0, 3])/Plot_Momentum[0, 3])
    plt.xlabel('time')
    plt.ylabel('(H-H0)/H0')
    plt.show()

    plt.figure()
    plt.plot(Plot_T, Plot_Momentum[:, 0:3])
    plt.xlabel('time')
    plt.ylabel('Momentum')
    plt.legend(['$H_1$','$H_2$','$H_3$'])
    plt.show()

def plot_kinetic(Plot_T,Plot_Kinetic):
    plt.figure()
    plt.plot(Plot_T, Plot_Kinetic)
    plt.xlabel('time')
    plt.ylabel('Kinetic')
    plt.show()

def plot_energy(Plot_T,Plot_Energy):
    plt.figure()
    plt.plot(Plot_T, Plot_Energy)
    plt.xlabel('time')
    plt.ylabel('Energy')
    plt.show()

def plot_omega(Plot_T, Plot_Vbeg, nbR3, row, col):
    plt.figure()
    for i in range(0, nbR3, 1):
        plt.subplot(row,col,i+1)
        plt.plot(Plot_T, Plot_Vbeg[:, i*6+3:i*6+6])
        plt.xlabel('time')
        plt.ylabel('Spin Rate')
        plt.title('body #'+str(i))
        plt.legend(['$\omega_1$','$\omega_2$','$\omega_3$'])
    plt.show()

def plot_eulerangle(Plot_T, Plot_EA):
    plt.figure()
    plt.plot(Plot_T, Plot_EA[:, 0:3])
    plt.xlabel('time')
    plt.ylabel('Euler Angle (rad)')
    plt.legend(['$\phi$','$\theta$','$\psi$'])
    plt.show()

def plot_inertia(Plot_T, Plot_inertia):
    plt.figure()
    plt.plot(Plot_T, Plot_inertia[:, 0:3])
    plt.xlabel('time')
    plt.ylabel('I_{body}')
    plt.legend(['$I_1$','$I_2$','$I_3$'])
    plt.show()

def plot_velocity(Plot_T, Plot_Vbeg, nbR3, row, col):
    plt.figure()
    for i in range(0, nbR3, 1):
        plt.subplot(row,col, i + 1)
        plt.plot(Plot_T, Plot_Vbeg[:, i * 6:i * 6 + 3])
        plt.xlabel('time')
        plt.ylabel('Velocity')
        plt.title('body #' + str(i))
        plt.legend(['$v_1$', '$v_2$', '$v_3$'])
    plt.show()

def plot_earthtraj(Plot_PosVecE):

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(Plot_PosVecE[:,0],Plot_PosVecE[:,1],Plot_PosVecE[:,2])
    ax.scatter3D(0,0,0,color='red')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_zlabel('z (km)')
    ax.set_title('Earth trajectory')
    plt.show()

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    Plot_T = np.linspace(0,5,5)
    Plot_Momentum = np.linspace(0,20,20)
    Plot_Momentum.resize(5,4)
    plot_momentum(Plot_T, Plot_Momentum)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
