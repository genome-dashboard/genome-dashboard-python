# -*- coding: utf-8 -*-

"""ds (data structure) module."""


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
        self.des = des


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

    def __init__(self, HP_intra, HP_inter):
        self.HP_intra = HP_intra
        self.HP_inter = HP_inter


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
    """Data Structure/Place Holder for any three dimensional data"""
    """A 3D Mask does not necessary to be a space curve, but an entry and an exit(could be the same) for a 3D Mask embedding on space curve is required"""
    """In this case, r, d are inputs as in Mask_3D"""
    dimension = 3

    def __init__(self, values, r, d, des=None):
        self.values = values
        self.r = r
        self.d = d
        self.des = des


class SC(object):
    """Data Structure for a Space Curve"""
    """A Space Curve contains all the information including Helical Parameters(HP), RD, mass(m), moment of inertia(I), Energy const(K)"""
    """Masks are input here as kwargs, which allows user input any numbers and any types of masks."""

    def __init__(self, HP=None, RD=None, m=None, I=None, K=None, **kwargs):
        self.HP = HP
        self.RD = RD
        self.m = m
        self.I = I
        self.K = K
        self.__dict__.update(kwargs)
