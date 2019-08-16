# -*- coding: utf-8 -*-

"""convert (data format conversions) module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
import numpy as np
import scipy.linalg as la

def module_name():
    print("Module: genomedashboard.convert.convert.py.")

def Ry(theta):
    ry=np.array([[np.cos(theta),0.,np.sin(theta)],[0.,1.,0.],[-np.sin(theta),0.,np.cos(theta)]])
    return ry

def Rz(theta):
    rz=np.array([[np.cos(theta),-np.sin(theta),0.],[np.sin(theta),np.cos(theta),0.],[0.,0.,1.]])
    return rz

def rodr2(k):
    """rodr2 is used in HP2RD, CURVES solution"""
    phi = np.linalg.norm(k)
    if(np.isclose(phi,0)):
        #print "zero rotation what kind of DNA is this?"
        return np.eye(3)
    k =  k/phi
    I = np.eye(3)
    X = np.array([[0,-k[2][0],k[1][0]],[k[2][0],0,-k[0][0]],[-k[1][0],k[0][0],0]])
    Q =  np.dot(np.cos(phi),I) + np.dot((1 - np.cos(phi)),np.dot(k,np.transpose(k)) )  +  np.dot(np.sin(phi),X)
    return Q

def HP2RD(shi, sli, ris, til, rol, twi, hptype='3DNA'):
    """
        Convert HP to RD, which require inputs of 6 inter HPs, and outputs R and D.
        Options of 3DNA and CURVES are provided in use of different types of HPs.
    """
    pi=np.pi/180
    til=til*pi
    rol=rol*pi
    twi=twi*pi
    r = np.zeros((3,1))
    D=np.array([[shi],[sli],[ris]])
    if hptype=='3DNA':
        bend=np.sqrt(til*til+rol*rol)
        if (bend == 0):
            phi = 0.0
        elif (rol == 0):
            phi = np.pi/2
        elif (rol < 0):
            phi = np.arctan(til/rol) + np.pi
        else:
            phi = np.arctan(til/rol)
        RZ1=Rz(twi*0.5-phi)
        Tmst = np.dot(RZ1,np.dot(Ry(bend*0.5),Rz(phi)))
        Ti = np.dot(RZ1,np.dot(Ry(bend),Rz(twi*0.5+phi)))
        r = r + np.dot(Tmst,D)
    elif hptype=='CURVES':
        u = np.array([[til],[rol],[twi]])
        Ti = rodr2(u)
        H = np.real(la.sqrtm(Ti))
        r = r + np.dot(H,D)
    else:
        print ('Please provide a valid type, "3DNA" or "CURVES"')
    return r,Ti


