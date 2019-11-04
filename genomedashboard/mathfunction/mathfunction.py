# -*- coding: utf-8 -*-

"""math functions module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
from convert import convert as cv
import numpy as np


def module_name():
    print("Module: genomedashboard.io.io.py.")

def straight_twisted_line(Rise, Twist, step_number):
    """
    Given Rise, Twist, and number of steps, return a list of HPs.
    """
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    for i in range(step_number):
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,Rise,0.0,0.0,Twist)))
    return hps

def circular_DNA(Rise, V1, V2, step_number, step_size=1.0, phase=0):
    """
    Given
    """
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    V1=V1/step_size
    Rise = Rise*step_size
    Twist = 360.0*V2/V1
    for s in range(step_number):
        s = s+phase
        Tilt = np.sin(Twist*s*np.pi/180.0)*360.0/V1
        Roll = np.cos(Twist*s*np.pi/180.0)*360.0/V1
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,Rise,Tilt,Roll,Twist)))
    return hps

def helix_torsion(Rise, Twist, V1, V2, step_number, step_size=1.0, phase=0):
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    V1=V1/step_size
    Twist=Twist*step_size
    Rise=Rise*step_size
    V2=V2*step_size
    for s in range(step_number):
        s = s+phase
        Tilt = np.sin((Twist+V2)*s*np.pi/180.0)*360.0/V1
        Roll = np.cos((Twist+V2)*s*np.pi/180.0)*360.0/V1
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,Rise,Tilt,Roll,Twist)))
    return hps

def helix_shear(Rise, Twist, V1, V2, step_number, step_size=1.0, phase=0):
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    V1=V1/step_size
    Twist=Twist*step_size
    Rise=Rise*step_size
    V2=V2*step_size
    for s in range(step_number):
        s = s+phase
        Tilt = np.sin(Twist*s*np.pi/180.0)*360.0/V1
        Roll = np.cos(Twist*s*np.pi/180.0)*360.0/V1
        Shift = V2*np.sin(Twist*s*np.pi/180)
        Slide = V2*np.cos(Twist*s*np.pi/180)
        hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(Shift,Slide,Rise,Tilt,Roll,Twist)))
    return hps

def circular_DNA_RD(Rise, V1, V2, step_number, step_size=1.0):
    Rise = Rise*step_size
    V1 = V1/step_size
    Radius = (Rise/2)/np.sin(np.pi/V1)
    Twist = 2*np.pi*V2/V1
    rd=[]
    for s in range(step_number):
        r = np.array([-Radius*np.cos(s*2*np.pi/V1),0,Radius*np.sin(s*2*np.pi/V1)])+np.array([Radius,0,0])
        d = np.dot(np.array([[np.cos(Twist*s),-np.sin(Twist*s),0],[np.sin(Twist*s),np.cos(Twist*s),0],[0,0,1]]).T,np.array([[np.cos(s*2*np.pi/V1),0,-np.sin(s*2*np.pi/V1)],[0,1,0],[np.sin(s*2*np.pi/V1),0,np.cos(s*2*np.pi/V1)]]))
        rd.append(ds.RD(r,d))
    return rd

def helix_torsion_RD(Rise, Twist, V1, V2, step_number, step_size):
    V1 = V1/step_size
    V2 = V2*step_size
    Rise = Rise*step_size
    Twist = Twist*step_size*np.pi/180
    Radius = np.sqrt(Rise**2/(4*(np.sin(np.pi/V1))**2+(np.tan(V2*np.pi/180))**2))
    h = Radius*np.tan(V2*np.pi/180)
    rd=[]
    for s in range(step_number):
        r = np.array([-Radius*np.cos(s*2*np.pi/V1),h*s,Radius*np.sin(s*2*np.pi/V1)])+np.array([Radius,0,0])
        tor = V2*np.pi/180
        phi = np.arctan(h/(np.sqrt(Rise**2-h**2)*np.cos(np.pi/V1)))
        X = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
        Y = np.array([[np.cos(s*2*np.pi/V1),0,-np.sin(s*2*np.pi/V1)],[0,1,0],[np.sin(s*2*np.pi/V1),0,np.cos(s*2*np.pi/V1)]])
        Z = np.array([[np.cos((Twist-tor)*s),-np.sin((Twist-tor)*s),0],[np.sin((Twist-tor)*s),np.cos((Twist-tor)*s),0],[0,0,1]])
        d = np.dot(Z.T,np.dot(X,Y))
        rd.append(ds.RD(r,d))
    return rd

def helix_shear_RD(Rise,Twist,V1,V2,step_number,step_size):
    V1 = V1/step_size
    V2 = V2*step_size
    Rise = Rise*step_size
    Twist = Twist*step_size*np.pi/180
    Radius = (np.sqrt(Rise**2-V2**2)/2)/np.sin(np.pi/V1)
    rd=[]
    for s in range(step_number):
        r = np.array([-Radius*np.cos(s*2*np.pi/V1),s*V2,Radius*np.sin(s*2*np.pi/V1)])+np.array([Radius,0,0])
        phi = np.arctan(V2/(np.sqrt(Rise**2-V2**2)*np.cos(np.pi/V1)))
        X = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
        Y = np.array([[np.cos(s*2*np.pi/V1),0,-np.sin(s*2*np.pi/V1)],[0,1,0],[np.sin(s*2*np.pi/V1),0,np.cos(s*2*np.pi/V1)]])
        Z = np.array([[np.cos(Twist*s),-np.sin(Twist*s),0],[np.sin(Twist*s),np.cos(Twist*s),0],[0,0,1]])
        d = np.dot(Z.T,np.dot(X,Y))
        rd.append(ds.RD(r,d))
    return rd
