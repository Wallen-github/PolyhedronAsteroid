# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ExternalTorques.py
@date: 7/18/23 16:10
@desc: 
"""
import math

# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
from lib.SubFunc import *
from pylmgc90 import chipy
from scipy.integrate import solve_ivp


def EarthAccel(posE, coor):
    AccelE = np.zeros([6])
    # Get unit
    UnitM, UnitL, UnitT, GG = Set_Unit()

    # Integrate Earth trajectory
    massE = 5.974227245203837E24 / UnitM  # kg, Earth mass

    # Compute the Earth acceleration
    r0i = posE - coor[0:3]
    AccelE[0:3] = GG * massE / (np.linalg.norm(r0i)) ** 3 * r0i

    return AccelE

def EarthPos(massA, k, dt, PosVecE0):

    # Get unit
    UnitM, UnitL, UnitT, GG = Set_Unit()

    # Integrate Earth trajectory
    massE = 5.974227245203837E24 / UnitM  # kg, Earth mass
    tol = 1E-13
    timespan = [(k - 1) * dt, k * dt]
    mu = GG*(massA+massE)
    PosVecSol = solve_ivp(fun=FlybyOrbit, t_span=timespan, y0=PosVecE0, args=(mu,), method='RK45',
                          rtol=tol, atol=tol)
    PosVec = PosVecSol.y[:, -1].copy()
    posE = PosVec[0:3]
    # vecE = PosVec[3:6]

    return posE, PosVec

def InitialPosVel_Earth(massA, time):
    """
    This function generate a real Apophis initial position before close approach
    :param massA: Asteroid mass
    :param time: Days before close approach
    :return: initial position and velocity
    """

    # Get unit
    UnitM, UnitL, UnitT, GG = Set_Unit()
    massE = 5.974227245203837E24 / UnitM  # kg, Earth mass

    # Compute the Apophis's initial position w.r.t Earth center
    # Closest Approach at Apr 13 2029 21:46 (38012 km radius)
    PosVec = np.array([-1.918693981897831E+04,3.225869496868482E+04,6.007972503646030E+03,
                         6.332252780541534E+00,3.405149229749669E+00,1.844426571255006E+00])
    PosVecCA = np.zeros([6])
    PosVecCA[0:3] = PosVec[0:3]*1E3/UnitL
    PosVecCA[3:6] = PosVec[3:6]*1E3*UnitT/UnitL

    timespan = [0, -time/2]
    tol = 1E-13
    mu = GG * (massA + massE)
    PosVecSol = solve_ivp(fun=FlybyOrbit, t_span=timespan, y0=PosVecCA, args=(mu,), method='RK45',
                          rtol=tol, atol=tol)

    # Get the Earth's initial position w.r.t Apophis center
    PosVecE0 = - PosVecSol.y[:, -1].copy()

    return PosVecE0

def InitialPosVel_Earth_v2(massA, vel_inf, pos_p, time):
    """
    This function generate a planar hyperbolic flyby orbit, which are defined by the encounter velocity at
    infinity, vel_inf, and the perigee distance, pos_p.
    :param massA: Asteroid mass
    :param vel_inf: the encounter velocity at infinity
    :param pos_p: the perigee distance
    :param time: Days before close approach
    :return: initial position and velocity
    """

    # For Apophis, we have pos_p = 38012 km, ecc = 4.229, axi = 1.152E4 km, vel_inf = 5.882 km/s

    # Get unit
    UnitM, UnitL, UnitT, GG = Set_Unit()
    massE = 5.974227245203837E24 / UnitM  # kg, Earth mass
    muE = GG*massE

    axi = muE/vel_inf**2
    ecc = pos_p/axi + 1
    f = 0
    inc = 0
    LOme = 0
    COme = 0
    p = axi*(ecc**2-1)
    r = p/(1+ecc*np.cos(f))
    # vel_p = np.sqrt(muE*(2.0/pos_p+1./axi))
    vel_p = np.sqrt(muE/p)*(1+ecc)

    Phat = np.array([np.cos(COme) * np.cos(LOme) - np.sin(LOme)*np.sin(LOme)*np.cos(inc),
                     np.sin(COme) * np.cos(LOme) + np.cos(LOme)*np.sin(LOme)*np.cos(inc),
                     np.sin(LOme)*np.sin(inc)])
    Qhat = np.array([-np.cos(COme) * np.sin(LOme) - np.sin(LOme)*np.cos(LOme)*np.cos(inc),
                     -np.sin(COme) * np.sin(LOme) + np.cos(LOme)*np.cos(LOme)*np.cos(inc),
                     np.cos(LOme)*np.sin(inc)])

    PosVecCA = np.zeros([6])
    PosVecCA[0:3] = r*np.cos(f)*Phat + r*np.sin(f)*Qhat
    PosVecCA[3:6] = np.sqrt(muE/p)*(-np.sin(f)*Phat + (np.cos(f) + ecc)*Qhat)

    PosVecCA[0:3] = PosVecCA[0:3] / UnitL
    PosVecCA[3:6] = PosVecCA[3:6] *UnitT/UnitL

    # PosVecCA[0:3] = np.array([pos_p, 0., 0.])/UnitL
    # PosVecCA[3:6] = np.array([0., -vel_p, 0.])*UnitT/UnitL

    timespan = [0, -time / 2]
    tol = 1E-13
    mu = GG * (massA + massE)
    PosVecSol = solve_ivp(fun=FlybyOrbit, t_span=timespan, y0=PosVecCA, args=(mu,), method='RK45',
                          rtol=tol, atol=tol)

    # Get the Earth's initial position w.r.t Apophis center
    PosVecE0 = - PosVecSol.y[:, -1].copy()


    return PosVecE0

def InitialPosVel_Asteroid(massA, time):

    # Get unit
    UnitM, UnitL, UnitT, GG = Set_Unit()
    massE = 5.974227245203837E24 / UnitM  # kg, Earth mass

    # Compute the Apophis's initial position w.r.t Earth center
    # Closest Approach at Apr 13 2029 21:46 (38012 km radius)
    PosVec = np.array([-1.918693981897831E+04,3.225869496868482E+04,6.007972503646030E+03,
                         6.332252780541534E+00,3.405149229749669E+00,1.844426571255006E+00])
    PosVecCA = np.zeros([6])
    PosVecCA[0:3] = PosVec[0:3]*1E3/UnitL
    PosVecCA[4:6] = PosVec[4:6]*1E3*UnitT/UnitL

    timespan = [0, -time/2]
    tol = 1E-13
    mu = GG * (massA + massE)
    PosVecSol = solve_ivp(fun=FlybyOrbit, t_span=timespan, y0=PosVecCA, args=(mu,), method='RK45',
                          rtol=tol, atol=tol)

    # Get the Earth's initial position w.r.t Apophis center
    PosVecA0 = PosVecSol.y[:, -1].copy()

    return PosVecA0

def FlybyOrbit(t,PosVec,mu):
    '''
    This function provide a EOM of two-body problem
    Args:
        t: time
        PosVec: position and velocity
        mu: mass, parameter, G * M

    Returns:
        Accel: acceleration
    '''

    Accel = np.zeros(6)

    # Zero order term
    Accel[0:3] = PosVec[3:6]
    Accel[3:6] = - mu * PosVec[0:3] / np.linalg.norm(PosVec[0:3])**3

    return Accel


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print('PyCharm')

    UnitM, UnitL, UnitT, GG = Set_Unit()
    massA = 3.970453830364233e+10/UnitM
    time = 4*86400/UnitT
    PosVecE0 = InitialPosVel_Earth(massA, time)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
