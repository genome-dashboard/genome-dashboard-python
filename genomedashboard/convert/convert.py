# -*- coding: utf-8 -*-

"""convert (data format conversions) module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
import numpy as np
import scipy.linalg as la
import copy

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

def skal(v1,v2):
    skal=0.0
    for k in range(3):
        skal=skal+v1[k]*v2[k]
    return skal

def vnormal(v):
    n = np.linalg.norm(v)
    if (np.isclose(n,0)):
        return v
    else:
        return v/n

def vprod(v1,v2):
    return np.array([v1[1]*v2[2]-v2[1]*v1[2],-v1[0]*v2[2]+v2[0]*v1[2],v1[0]*v2[1]-v2[0]*v1[1]])

def vrot(v,vn,fir):
    vout=np.zeros(3)
    vpom=vprod(vn,v)
    for k in range(3):
        vout[k]=v[k]*np.cos(fir)+vn[k]*skal(vn,v)*(1.0-np.cos(fir))+vpom[k]*np.sin(fir)
    return vout

def HP2RD(HP, hptype='3DNA'):
    """
    Convert HP to RD, which require inputs of 6 inter HPs, and outputs R and D.
    Options of 3DNA and CURVES are provided in use of different types of HPs.
    """
    pi=np.pi/180
    til=HP.HP_inter.til*pi
    rol=HP.HP_inter.rol*pi
    twi=HP.HP_inter.twi*pi
    r = np.zeros((3,1))
    D=np.array([[HP.HP_inter.shi],[HP.HP_inter.sli],[HP.HP_inter.ris]])
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
    rd = ds.RD(r.T.reshape(3),Ti.T)
    return rd

def RD_stack(rd1, rd2, value=None):
    """
    Given two steps of RD, previous one is global RD, and latter one is local one.
    Calculate the global RD for the latter one.
    If a 3D Mask is given, then the global value of the Mask can be calculated.
    """
    R = np.dot(rd2.r,rd1.d) + rd1.r
    D = np.dot(rd2.d,rd1.d)
    rd = ds.RD(R,D)
    if value is None:
        return rd
    else:
        value = np.dot(value,rd1.d) + rd1.r
        return rd, value

def RD_loc2g(rdlist):
    """
    Given a list of local RD, generate the global RD used in SC(Space Curve).
    """
    rd = [rdlist[0]]
    for i in range(1,len(rdlist)):
        rd.append(RD_stack(rd[i-1],rdlist[i]))
    return rd

def RD2HP(rd1,rd2,hptype='3DNA'):
    """
    Given two global rd steps, calculate HP.
    Has two options, 3DNA and CURVES
    """
    if hptype=='3DNA':
        x1 = vnormal(rd1.d[0])
        y1 = vnormal(rd1.d[1])
        z1 = vnormal(rd1.d[2])
        x2 = vnormal(rd2.d[0])
        y2 = vnormal(rd2.d[1])
        z2 = vnormal(rd2.d[2])
        pgama=skal(z1,z2)
        if (pgama>0.99999):
            pgama=1.0
        if (pgama<-0.99999):
            pgama=-1.0
        gama=np.arccos(pgama)
        if (z1[0]==z2[0] and z1[1]==z2[1] and z1[2]==z2[2]):
            rt=np.zeros(3)
            x1p=copy.copy(x1)
            x2p=copy.copy(x2)
            y1p=copy.copy(y1)
            y2p=copy.copy(y2)
        else:
            rt=vnormal(vprod(z1,z2))
            x1p=vrot(x1,rt,gama/2)
            y1p=vrot(y1,rt,gama/2)
            x2p=vrot(x2,rt,-gama/2)
            y2p=vrot(y2,rt,-gama/2)
        xm=vnormal(x1p+x2p)
        ym=vnormal(y1p+y2p)
        zm=vnormal(z1+z2)
        pomega=skal(y1p,y2p)
        if (pomega>0.99999):
            pomega=1.0
        if (pomega<-0.99999):
            pomega=-1.0
        omega=np.arccos(pomega)
        ypom=vprod(y1p,y2p)
        if (skal(ypom,zm)<0.0):
            omega=-omega
        pfi=skal(rt,ym)
        if (pfi>0.99999):
            pfi=1.0
        if (pfi<-0.99999):
            pfi=-1.0
        fi=np.arccos(pfi)
        rpom=vprod(rt,ym)
        if (skal(rpom,zm)<0.0):
            fi=-fi
        rol=gama*np.cos(fi)*180.0/np.pi
        til=gama*np.sin(fi)*180.0/np.pi
        disp=rd2.r-rd1.r
        shi=skal(disp,xm)
        sli=skal(disp,ym)
        ris=skal(disp,zm)
        twi=omega*180.0/np.pi
        interHP = ds.HP_inter(shi,sli,ris,til,rol,twi)
    elif hptype=='CURVES':
        G1=rd1.d
        G2=rd2.d
        dq=rd2.r-rd1.r
        L = np.dot(G2,np.linalg.inv(G1))
        Q = L - L.T
        phi = np.arccos((np.trace(L)-1.0)/2.0)
        theta = (1.0/(1.0+np.trace(L)))*np.array([[Q[2][1]],[Q[0][2]],[Q[1][0]]])
        scale = np.tan(phi/2.0)
        if scale==0:
            trt = np.zeros((3,1))
        else:
            trt = -(phi*180/np.pi)*theta/scale
        H = np.real(la.sqrtm(L))
        ssr = np.dot(dq.T,np.dot(G1.T,H.T))
        interHP = ds.HP_inter(ssr[0],ssr[1],ssr[2],trt[0][0],trt[1][0],trt[2][0])
    else:
        print('Please provide a valid type, "3DNA" or "CURVES"')
    return ds.HP(None,interHP,hptype=hptype)

def HP_T(HP,T):
    T=np.sqrt(T/298.0)
    HP.HP_inter.shi = np.random.normal(HP.HP_inter.shi, T*0.76)
    HP.HP_inter.sli = np.random.normal(HP.HP_inter.sli, T*0.68)
    HP.HP_inter.ris = np.random.normal(HP.HP_inter.ris, T*0.37)
    HP.HP_inter.til = np.random.normal(HP.HP_inter.til, T*4.6)
    HP.HP_inter.rol = np.random.normal(HP.HP_inter.rol, T*7.2)
    HP.HP_inter.twi = np.random.normal(HP.HP_inter.twi, T*7.3)
    return HP

def SEQ2HP(seq,HP_dic,occ=[],nuc_type=[],T=0):
    """
        Given sequence, dictionary of HP(e.g. {A-A: [HP], oct: [HP...HP]}), occupancy(e.g. [1, 500 , 789]), nucleosome type at each occupancy(e.g. ['oct','tet','hex']), and temperature.
        Return the HPs associate with the given sequence.
    """
    seqstep=seq.tostep()
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    j=0
    while j<len(seqstep):
        if j in occ:
            hps.extend(HP_dic[nuc_type[occ.index(j)]][1:])
            j=j+len(HP_dic[nuc_type[occ.index(j)]])-1
        else:
            hps.append(HP_T(HP_dic[seqstep[j]][0],T))
            j=j+1
    return hps

def SEQ2RD(seq, RD):
    """
    Given sequence and RD,
    Return RDs associate with the given sequence.
    """
    

def elastic_energy(seq,HP_free,HP_nuc,K):
    """
    Given the HP for free DNA, and the HP for nucleosome DNA, and K.
    Calculate the elastic energy.
    """
    seqstep=seq.tostep()
    Enum = len(HP_free)-len(HP_nuc)+1
    E = np.zeros(Enum)
    p = np.array([[i.HP_inter.shi,i.HP_inter.sli,i.HP_inter.ris,i.HP_inter.til,i.HP_inter.rol,i.HP_inter.twi] for i in HP_nuc[1:]])
    d = np.array([[i.HP_inter.shi,i.HP_inter.sli,i.HP_inter.ris,i.HP_inter.til,i.HP_inter.rol,i.HP_inter.twi] for i in HP_free[1:]])
    for i in range(0,Enum):
        pd = p-d[i:i+len(p)]
        x = [0.5*np.dot(np.dot(j,K[seqstep[i+idx]]),j.T) for idx,j in enumerate(pd)]
        E[i] = sum(x)
    return E

def E2Occ(seq_length, nuc_nbp, E, occu, lk, phase=0):
    """
    Given the length of sequence, number of base pairs of the nucleosome, occupancy percentage,
    linker length, and phase.
    Calculate the occupancy according to the Energy.
    Energy could be replaced by any 1D Mask, and the occupancy of local lowest Energy(1D informatics) in then calculated
    """
    denominator = nuc_nbp + lk
    Etmp=copy.copy(E)
    if occu==1:
        numnuc = int((seq_length-phase+lk)/denominator)
        occ = np.zeros(numnuc)
        for i in range(numnuc):
            occ[i] = phase + i*denominator
    else:
        numnuc = int(seq_length*occu/denominator)
        occ = np.zeros(numnuc)
        i=0
        while i<numnuc:
            tmpmin=np.argmin(Etmp)
            if i>0:
                k=0
                for x in occ:
                    if np.abs(tmpmin-x)<denominator:
                        k+=1
                        break
                if k==0:
                    occ[i]=tmpmin
                    Etmp[tmpmin]=np.inf
                    i+=1
                else:
                    Etmp[tmpmin]=np.inf
            else:
                occ[i]=tmpmin
                Etmp[tmpmin]=np.inf
                i+=1
    return [int(x+1) for x in occ]


##########Two Angle Model###########

def deflection(xyz1,xyz2,xyz3):
    """
    Calculate deflection angle(alpha)
    """
    v1=xyz1-xyz2;
    v2=xyz3-xyz2;
    alpha = np.arccos(np.dot(v1,v2)/(np.sqrt(v1[0]**2+v1[1]**2+v1[2]**2)*(np.sqrt(v2[0]**2+v2[1]**2+v2[2]**2))))*180/np.pi
    return alpha

def dihedral(xyz1,xyz2,xyz3,xyz4):
    """
    Calculate dihedral angle(beta)
    """
    q1=xyz2-xyz1
    q2=xyz3-xyz2
    q3=xyz4-xyz3
    q1_x_q2 = np.cross(q1,q2)
    q2_x_q3 = np.cross(q2,q3)
    n1 = q1_x_q2/np.sqrt(np.dot(q1_x_q2,q1_x_q2))
    n2 = q2_x_q3/np.sqrt(np.dot(q2_x_q3,q2_x_q3))
    u1 = n2
    u3 = q2/(np.sqrt(np.dot(q2,q2)))
    u2 = np.cross(u3,u1)
    cos_phi = np.dot(n1,u1)
    sin_phi = np.dot(n1,u2)
    beta = -math.atan2(sin_phi,cos_phi)*180/np.pi
    return beta

def two_angle_plot(Mask_3D):
    """
    Given 3D Mask and plot the two angle plot.
    """

