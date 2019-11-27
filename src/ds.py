# -*- coding: utf-8 -*-

"""ds (data structure) module."""

import numpy as np

def module_name():
    print("Module: genomedashboard.ds.ds.py.")


class HP_intra(object):
    """Data Structure for intra Helical parameters, including Shear, Stretch, Stagger, Buckle, Prop-Tw and Opening."""

    def __init__(self, she, str, sta, buc, pro, ope):
        self.she = she
        self.str = str
        self.sta = sta
        self.buc = buc
        self.pro = pro
        self.ope = ope


class HP_inter(object):
    """Data Structure for inter Helical parameters, including Shift, Slide, Rise, Tilt, Roll, Twist."""

    def __init__(self, shi, sli, ris, til, rol, twi):
        self.shi = shi
        self.sli = sli
        self.ris = ris
        self.til = til
        self.rol = rol
        self.twi = twi


class HP(object):
    """Data Structure for Helical Parameters, including intra and inter"""

    def __init__(self, HP_intra, HP_inter, hptype='3DNA'):
        self.HP_intra = HP_intra
        self.HP_inter = HP_inter
        self.hptype = hptype


class RD(object):
    """RD includes: r, a set of coordinates of points on the centerline of a space curve, and d, the direction frame of these points."""

    def __init__(self, r, d):
        self.r = r
        self.d = d


class Mask_1D(object):
    """Data Structure/Place Holder for any one dimensional data"""
    dimension = 1

    def __init__(self, values, des=None):
        self.values = values
        self.des = des


class Mask_2D(object):
    """Data Structure/Place Holder for any two dimensional data"""
    dimension = 2

    def __init__(self, values, des=None):
        self.values = values
        self.des = des


class Mask_3D(object):
    """
    Data Structure/Place Holder for any three dimensional data.
    A 3D Mask does not necessary to be a space curve, but an entry and an exit(could be the same) for a 3D Mask embedding on space curve is required
    In this case, r, d are inputs as in Mask_3D
    """
    dimension = 3

    def __init__(self, values, RD_entry=RD(np.zeros(3),np.eye(3)), RD_exit=RD(np.zeros(3),np.eye(3)), skip=0, des=None):
        self.values = values
        self.RD_entry = RD_entry
        self.RD_exit = RD_exit
        self.skip = skip
        self.des = des


class SC(object):
    """
    Data Structure for a Space Curve
    A Space Curve contains all the information including Helical Parameters(HP), RD, mass(m), moment of inertia(I), Energy const(K)
    Masks are input here as kwargs, which allows user input any numbers and any types of masks.
    """

    def __init__(self, HP=None, RD=None, m=None, I=None, K=None, **kwargs):
        self.HP = HP
        self.RD = RD
        self.m = m
        self.I = I
        self.K = K
        self.__dict__.update(kwargs)

class SEQ(object):
    """
    Data Structure for DNA sequence,
    GDash model is sequence specified model, sequence has its own data structure as a special 1D Mask.
    """
    def __init__(self, seq):
        self.seq=seq

    def tolist(self):
        if type(self.seq)==list:
            return self.seq
        else:
            seqlist = list(self.seq)
            seqlist = ['A-T' if x=='A' else x for x in seqlist]
            seqlist = ['T-A' if x=='T' else x for x in seqlist]
            seqlist = ['C-G' if x=='C' else x for x in seqlist]
            seqlist = ['G-C' if x=='G' else x for x in seqlist]
            return seqlist

    def tostring(self):
        if type(self.seq)==str:
            return self.seq
        else:
            seqstr=''.join(x[0] for x in self.seq)
            return seqstr

    def tostep(self):
        sequence=self.tostring()
        step=[]
        for i in range(len(sequence)-1):
            step.append(sequence[i]+'-'+sequence[i+1])
        return step

class PDB_std(object):
    """
    Data structure for standard pdb format.
    ATOM and CONECT are the only two components that take consider in current version.
    PDB_std is a container that hold two lines of PDB information(atom+connection)
    """
    def __init__(self, atom=None, serial=None, name=None, altLoc=None, resName=None, chainID=None, resSeq=None, iCode=None, x=None, y=None, z=None, occupancy=None, tempFactor=None, segID=None,element=None, charge=None, CONECT=None):
        self.atom=atom
        self.serial=serial
        self.name=name
        self.altLoc=altLoc
        self.resName=resName
        self.chainID=chainID
        self.resSeq=resSeq
        self.iCode=iCode
        self.x=x
        self.y=y
        self.z=z
        self.occupancy=occupancy
        self.tempFactor=tempFactor
        self.segID=segID
        self.element=element
        self.charge=charge
        self.CONECT=CONECT
