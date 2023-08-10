# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ExternalTorques.py
@date: 7/18/23 16:10
@desc: 
"""


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
