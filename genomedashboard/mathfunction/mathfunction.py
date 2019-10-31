# -*- coding: utf-8 -*-

"""math functions module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
from convert import convert as cv
import numpy as np

def straight_twisted_line(Rise, Twist, step_number):
    """
    Given Rise, Twist, and number of steps, return a list of HPs.
    """
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    for i in range(step_number):
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,Rise,0.0,0.0,Twist)))
    return hps

def circular_DNA(Rise, V1, V2, step_number, step_size=1.0):
    """
    Given
    """
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    V1=V1/step_size
    Rise = Rise*step_size
    Twist = 360.0*V2/V1
    for i in range(step_number):
        s = i
        Tilt = np.sin(Twist*s*np.pi/180.0)*360.0/V1
        Roll = np.cos(Twist*s*np.pi/180.0)*360.0/V1
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,Rise,Tilt,Roll,Twist)))
    return hps

def helix_torsion(Rise, Twist, V1, V2, step_number, step_size=1.0):
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    V1=V1/step_size
    Twist=Twist*step_size
    Rise=Rise*step_size
    V2=V2*step_size
    for i in range(step_number):
        s = i
        Tilt = np.sin((Twist+V2)*s*np.pi/180.0)*360.0/V1
        Roll = np.cos((Twist+V2)*s*np.pi/180.0)*360.0/V1
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,Rise,Tilt,Roll,Twist)))
    return hps

def helix_shear(Rise, Twist, V1, V2, step_number, step_size=1.0):
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    V1=V1/step_size
    Twist=Twist*step_size
    Rise=Rise*step_size
    V2=V2*step_size
    for i in range(step_number):
        s = i
        Tilt = np.sin(Twist*s*np.pi/180.0)*360.0/V1
        Roll = np.cos(Twist*s*np.pi/180.0)*360.0/V1
        Shift = V2*np.sin(Twist*s*np.pi/180)
        Slide = V2*np.cos(Twist*s*np.pi/180)
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(Shift,Slide,Rise,Tilt,Roll,Twist)))
    return hps

def circular_DNA_RD(Radius, V1, V2, step_number, step_size=1.0):
    V1 = V1/step_size
    Twist = 2*np.pi*V2/V1
    rd=[]
    for s in range(step_number):
        r = np.array([-Radius*np.cos(s*2*np.pi/V1),0,Radius*np.sin(s*2*np.pi/V1)])+np.array([Radius,0,0])
        d = np.dot(np.array([[np.cos(Twist*s),-np.sin(Twist*s),0],[np.sin(Twist*s),np.cos(Twist*s),0],[0,0,1]]),np.array([[np.cos(s*2*np.pi/V1),0,-np.sin(s*2*np.pi/V1)],[0,1,0],[np.sin(s*2*np.pi/V1),0,np.cos(s*2*np.pi/V1)]]))
        rd.append(ds.RD(r,d))
    return rd
