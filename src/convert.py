# -*- coding: utf-8 -*-

"""convert (data format conversions) module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
import numpy as np
import scipy.linalg as la
from scipy.integrate import odeint
import matplotlib.pyplot as plt

"""
NOTE: There are good reasons for why scientific libraries don't usually have matplotlib as a direct dependency:
it is very platform dependent and it is usually not essential for what these libraries are trying to accomplish.
In our case we can try and handle the OS X sidecases with a try/except block (see below).
"""
# try:
#    import matplotlib.pyplot as plt
# except RuntimeError as e:
#     # Handle OSX exceptions fopr matplotlib.
#    if 'Python is not installed as a framework.' in e.message:
#      # warnings.warn(.. some warning about disabled plotting...)
#      import matplotlib
#      # matplotlib.use('PS')
#      matplotlib.use('TkAgg')
#      import matplotlib.pyplot as plt
"""
Solution is to get users to include an rc config file in the pip library after installaiton:
Create a file $VENV/../site-packages/matplotlib/matplotlibrc there and add the following code: `backend: TkAgg`
"""

import math
import copy


def module_name():
    print("Module: genomedashboard.convert.convert.py.")

##########################
####Preprocessor macro####
##########################

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

def odeSC(y,s,hp_list):
    """
    The set up for the ode function of Space Curve. dr/ds=D*Gamma and ddi/ds=(D*Omega)xdi.
    Gamma = [Shift, Slide, Rise], Omega=[Tilt, Roll, Twist]
    """
    pi=np.pi/180
    hp=hp_list[int(s)]
    til=hp.HP_inter.til*pi
    rol=hp.HP_inter.rol*pi
    twi=hp.HP_inter.twi*pi
    Dmat=y.reshape(12,1)[3:12].reshape(3,3)
    gamma = np.dot(Dmat.T,np.array([[hp.HP_inter.shi],[hp.HP_inter.sli],[hp.HP_inter.ris]]))
    omega = np.dot(Dmat.T,np.array([[til],[rol],[twi]]))
    #omega = np.dot(Dmat.T,np.sin(np.array([[til],[rol],[twi]])))
    dydt = np.zeros((4,3))
    dydt[0] = gamma.reshape(1,3)
    dydt[1:4] = np.cross(omega.reshape(1,3),Dmat)
    #dydt[1:4] = 2*np.sin(np.cross(omega.reshape(1,3),Dmat)/2)
    return dydt.reshape(12,)

def odeSC_d(y,s,hp_list):
    """
    The set up for the ode function of Space Curve. ddi/ds=(D*Omega)xdi.
    Omega=[Tilt, Roll, Twist]
    """
    if s>len(hp_list):
        s = len(hp_list)-1
    pi=np.pi/180
    hp=hp_list[int(s)]
    til=hp.HP_inter.til*pi
    rol=hp.HP_inter.rol*pi
    twi=hp.HP_inter.twi*pi
    Dmat=y.reshape(9,1).reshape(3,3)
    #gamma = np.dot(Dmat.T,np.array([[hp.HP_inter.shi],[hp.HP_inter.sli],[hp.HP_inter.ris]]))
    omega = np.dot(Dmat.T,np.array([[til],[rol],[twi]]))
    #omega = np.dot(Dmat.T,np.sin(np.array([[til],[rol],[twi]])))
    dydt = np.cross(omega.reshape(1,3),Dmat)
    #dydt[1:4] = 2*np.sin(np.cross(omega.reshape(1,3),Dmat)/2)
    return dydt.reshape(9,)

def odeSC_r(y,s,hp_list,d):
    """
    The set up for the ode function of Space Curve. dr/ds=D*Gamma
    Gamma = [Shift, Slide, Rise]
    """
    pi=np.pi/180
    if s>len(hp_list):
        s=len(hp_list)-1
    hp=hp_list[int(s)]
    til=hp.HP_inter.til*pi
    rol=hp.HP_inter.rol*pi
    twi=hp.HP_inter.twi*pi
    #Dmat=np.dot(np.real(la.sqrtm(np.dot(d[int(s+1)],d[int(s)].T))),d[int(s)])
    Dmat=d[int(s)]
    gamma = np.dot(Dmat.T,np.array([[hp.HP_inter.shi],[hp.HP_inter.sli],[hp.HP_inter.ris]]))
    #omega = np.dot(Dmat.T,np.array([[til],[rol],[twi]]))
    #omega = np.dot(Dmat.T,np.sin(np.array([[til],[rol],[twi]])))
    dydt = gamma.reshape(1,3)
    #dydt[1:4] = 2*np.sin(np.cross(omega.reshape(1,3),Dmat)/2)
    return dydt.reshape(3,)


##########################
######Basic Functions#####
##########################

def HP2RD(HP, hptype='3DNA'):
    """
    Convert HP to RD, which require inputs of HPs(HP_inter is required, HP_intra could be None), and outputs R and D.
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
        sys.exit(0)
    rd = ds.RD(r.T.reshape(3),Ti.T)
    return rd

def RD_stack(rd1, rd2):
    """
    Given two steps of RD, previous one is global RD, and latter one is local one.
    Calculate the global RD for the latter one.
    If a 3D Mask is given, then the global value of the Mask can be calculated.
    """
    R = np.dot(rd2.r,rd1.d) + rd1.r
    D = np.dot(rd2.d,rd1.d)
    rd = ds.RD(R,D)
    return rd

def docking_Mask_3D(rd, Mask_3D):
    """
    Docking 3D Mask values onto rd
    """
    entry = copy.copy(Mask_3D.RD_entry)
    exit = copy.copy(Mask_3D.RD_exit)
    value = copy.copy(Mask_3D.values)
    entrydT = entry.d.T
    value = np.dot(value - entry.r,entrydT)
    exit.r = np.dot(exit.r - entry.r,entrydT)
    exit.d = np.dot(exit.d,entrydT)
    value = np.dot(value,rd.d) + rd.r
    entry.r = rd.r
    entry.d = rd.d
    exit.r = np.dot(exit.r,rd.d) +rd.r
    exit.d = np.dot(exit.d,rd.d)
    Mask_new = ds.Mask_3D(value,RD_entry=entry,RD_exit=exit,skip=Mask_3D.skip,des=Mask_3D.des)
    return Mask_new

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
        if (pgama>1.0):
            pgama=1.0
        if (pgama<-1.0):
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
        if (pomega>1.0):
            pomega=1.0
        if (pomega<-1.0):
            pomega=-1.0
        omega=np.arccos(pomega)
        ypom=vprod(y1p,y2p)
        if (skal(ypom,zm)<0.0):
            omega=-omega
        pfi=skal(rt,ym)
        if (pfi>1.0):
            pfi=1.0
        if (pfi<-1.0):
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
        sys.exit(0)
    return ds.HP(None,interHP,hptype=hptype)

def HP_T(HP,T):
    """
    Temperature
    """
    HP_new=copy.deepcopy(HP)
    T=np.sqrt(T/298.0)
    HP_new.HP_inter.shi = np.random.normal(HP_new.HP_inter.shi, T*0.76)
    HP_new.HP_inter.sli = np.random.normal(HP_new.HP_inter.sli, T*0.68)
    HP_new.HP_inter.ris = np.random.normal(HP_new.HP_inter.ris, T*0.37)
    HP_new.HP_inter.til = np.random.normal(HP_new.HP_inter.til, T*4.6)
    HP_new.HP_inter.rol = np.random.normal(HP_new.HP_inter.rol, T*7.2)
    HP_new.HP_inter.twi = np.random.normal(HP_new.HP_inter.twi, T*7.3)
    return HP_new

def SEQ2HP(seq,HP_dic,occ=[],nuc_type=[],T=0):
    """
    Given sequence, dictionary of HP(e.g. {A-A: [HP], oct: [HP...HP],...}), occupancy(e.g. [1, 500 , 789]), nucleosome type at each occupancy(e.g. ['oct','tet','hex']), and temperature.
    Return the HPs associate with the given sequence.
    """
    seqstep=seq.tostep()
    hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
    j=0
    while j<len(seqstep):
        hps.append(HP_T(HP_dic[seqstep[j]][0],T))
        j=j+1
    for i in occ:
        hps[i:i+len(HP_dic[nuc_type[occ.index(i)]])-1] = HP_dic[nuc_type[occ.index(i)]][1:]
    return hps

def SEQ2RD(seq, RD_dic,occ=[],nuc_type=[]):
    """
    Given sequence, dictionary of local RD(e.g. {A-A: [RD], oct:[RD...RD],...}),
    occupancy(e.g. [1, 500 , 789]), nucleosome type at each occupancy(e.g. ['oct','tet','hex'])
    Return RDs associate with the given sequence.
    """
    seqstep=seq.tostep()
    rd=[ds.RD(np.zeros(3),np.eye(3))]
    j=0
    while j<len(seqstep):
        rd.append(RD_dic[seqstep[j]][-1])
        j=j+1
    for i in occ:
        rd[i:i+len(RD_dic[nuc_type[occ.index(i)]])-1] = RD_dic[nuc_type[occ.index(i)]][1:]
    return rd

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
    occ_out=[int(x+1) for x in occ]
    occ_out.sort()
    return occ_out

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

def alpha_beta_list(xyz_list):
    """
    Given a list of ordered xyz coord, return alpha and beta
    """
    if len(xyz_list)<4:
        print('At least four elements required.')
        sys.exit(0)
    else:
        alpha=[]
        beta=[]
        for i in range(len(xyz_list)-3):
            alpha.append(deflection(xyz_list[i+1],xyz_list[i+2],xyz_list[i+3])*np.pi/180)
            beta.append(dihedral(xyz_list[i],xyz_list[i+1],xyz_list[i+2],xyz_list[i+3])*np.pi/180)
    return alpha,beta

######Distance Matrix#####

def distance_matrix(xyz_list, cut=0):
    """
    Given a list of ordered xyz coord, return distance matrix.
    cut is the minimal value that treat as colliding.
    """
    if len(xyz_list)<2:
        print('At least two elements required.')
        sys.exit(0)
    else:
        distance_matrix = np.zeros((len(xyz_list),len(xyz_list)))
        for i in range(len(xyz_list)):
            for j in range(len(xyz_list)):
                v = xyz_list[i]-xyz_list[j]
                distance_matrix[i][j] = np.sqrt(v[0]**2+v[1]**2+v[2]**2)
                if distance_matrix[i][j]<cut:
                    distance_matrix[i][j]=0.0
    return distance_matrix


##########################
######Specical Usages#####
##########################

def RD_loc2g(rdlist):
    """
    Given a list of local RD, generate the global RD used in SC(Space Curve).
    """
    rd = [rdlist[0]]
    for i in range(1,len(rdlist)):
        rd.append(RD_stack(rd[i-1],rdlist[i]))
    return rd

def HP2SC(hp_list,hptype='3DNA'):
    """
    Given a list of HP, calculate global RD.
    Return the Space Curve with both HP and RD.
    """
    if hptype=='3DNA' or hptype=='CURVES':
        local_rd = [HP2RD(hp_list[i],hptype) for i in range(len(hp_list))]
        rd_list = RD_loc2g(local_rd)
    elif hptype=='MATH':
        new_list=hp_list[1:]
        new_list.append(hp_list[0])
        y0d = np.eye(3)
        t=[i for i in range(len(hp_list))]
        yd = odeint(odeSC_d,y0d.reshape(9,),t,args=(new_list,))
        d = [i.reshape(3,3) for i in yd]
        d.append(np.eye(3))
        d_half = [np.dot(np.real(la.sqrtm(np.dot(d[int(s+1)],d[int(s)].T))),d[int(s)]) for s in range(len(d)-1)]
        y0r = np.zeros((1,3))
        tr=[i for i in range(len(hp_list))]
        yr = odeint(odeSC_r,y0r.reshape(3,),tr,args=(new_list,d_half,))
        rd_list = [ds.RD(yr[i],d[i]) for i in range(len(yr))]
        '''
    elif hptype=='MATH_3DNA':
        new_list=hp_list[1:]
        new_list.append(hp_list[0])
        new_list_d = [copy.copy(x) for x in new_list]
        new_list_d_half=[]
        for i in range(len(new_list_d)):
            new_list_d[i].HP_inter.til = new_list_d[i].HP_inter.til/2
            new_list_d[i].HP_inter.rol = new_list_d[i].HP_inter.rol/2
            new_list_d[i].HP_inter.twi = new_list_d[i].HP_inter.twi/2
            new_list_d_half.append(new_list_d[i])
            new_list_d_half.append(new_list_d[i])
        y0d = np.eye(3)
        t=[i for i in range(len(new_list_d_half))]
        yd = odeint(odeSC_d,y0d.reshape(9,),t,args=(new_list_d_half,))
        d = [j.reshape(3,3) for i,j in enumerate(yd) if i%2==1]
        d_half = [j.reshape(3,3) for i,j in enumerate(yd) if i%2==0]
        y0r = np.zeros((1,3))
        tr=[i for i in range(len(hp_list))]
        yr = odeint(odeSC_r,y0r.reshape(3,),tr,args=(new_list,d_half,))
        rd_list = [ds.RD(yr[i],d[i]) for i in range(len(yr))]
        '''
    else:
        print('Please provide a valid type, "3DNA", "CURVES" or "MATH"')
        sys.exit(0)
    return ds.SC(HP=hp_list,RD=rd_list)

def RD2SC(rd_list,hptype='3DNA',step_size=1):
    """
    Given a list of global RD, calculate HP.
    step_size is the step size between two steps.
    Return the Space Curve with both HP and RD
    """
    if hptype=='3DNA' or hptype=='CURVES':
        hp_list = [RD2HP(rd_list[i],rd_list[i],hptype) if i==0 else RD2HP(rd_list[i-1],rd_list[i],hptype) for i in range(len(rd_list))]
    elif hptype=='MATH':
        s = np.linspace(0,(len(rd_list)-1)*step_size,len(rd_list))
        rx = np.array([x.r[0] for x in rd_list])
        ry = np.array([x.r[1] for x in rd_list])
        rz = np.array([x.r[2] for x in rd_list])
        tx = np.diff(rx)/np.diff(s)
        ty = np.diff(ry)/np.diff(s)
        tz = np.diff(rz)/np.diff(s)
        t = np.array([tx,ty,tz]).T
        ssr = np.array([np.dot(t[i],np.dot(np.real(la.sqrtm(np.dot(rd_list[i+1].d,rd_list[i].d.T))),rd_list[i].d).T) for i in range(len(s)-1)])
        d1x = np.array([x.d[0][0] for x in rd_list])
        d1y = np.array([x.d[0][1] for x in rd_list])
        d1z = np.array([x.d[0][2] for x in rd_list])
        d2x = np.array([x.d[1][0] for x in rd_list])
        d2y = np.array([x.d[1][1] for x in rd_list])
        d2z = np.array([x.d[1][2] for x in rd_list])
        d3x = np.array([x.d[2][0] for x in rd_list])
        d3y = np.array([x.d[2][1] for x in rd_list])
        d3z = np.array([x.d[2][2] for x in rd_list])
        td1x = np.diff(d1x)/np.diff(s)
        td1y = np.diff(d1y)/np.diff(s)
        td1z = np.diff(d1z)/np.diff(s)
        td2x = np.diff(d2x)/np.diff(s)
        td2y = np.diff(d2y)/np.diff(s)
        td2z = np.diff(d2z)/np.diff(s)
        td3x = np.diff(d3x)/np.diff(s)
        td3y = np.diff(d3y)/np.diff(s)
        td3z = np.diff(d3z)/np.diff(s)
        deltad1 = np.array([td1x,td1y,td1z]).T
        deltad2 = np.array([td2x,td2y,td2z]).T
        deltad3 = np.array([td3x,td3y,td3z]).T
        Twist = [2*np.arcsin(np.dot(np.dot(deltad1[i],np.dot(np.real(la.sqrtm(np.dot(rd_list[i+1].d,rd_list[i].d.T))),rd_list[i].d).T),np.array([0,1,0]))/2)*180/np.pi for i in range(len(s)-1)]
        Roll = [2*np.arcsin(np.dot(np.dot(deltad3[i],np.dot(np.real(la.sqrtm(np.dot(rd_list[i+1].d,rd_list[i].d.T))),rd_list[i].d).T),np.array([1,0,0]))/2)*180/np.pi for i in range(len(s)-1)]
        Tilt = [2*np.arcsin(np.dot(np.dot(deltad2[i],np.dot(np.real(la.sqrtm(np.dot(rd_list[i+1].d,rd_list[i].d.T))),rd_list[i].d).T),np.array([0,0,1]))/2)*180/np.pi for i in range(len(s)-1)]
        #Twist = [np.dot(np.dot(2*np.arcsin(deltad1[i]/2),np.dot(np.real(la.sqrtm(np.dot(rd_list[i+1].d,rd_list[i].d.T))),rd_list[i].d).T),np.array([0,1,0]))*180/np.pi for i in range(len(s)-1)]
        #Roll = [np.dot(np.dot(2*np.arcsin(deltad3[i]/2),np.dot(np.real(la.sqrtm(np.dot(rd_list[i+1].d,rd_list[i].d.T))),rd_list[i].d).T),np.array([1,0,0]))*180/np.pi for i in range(len(s)-1)]
        #Tilt = [np.dot(np.dot(2*np.arcsin(deltad2[i]/2),np.dot(np.real(la.sqrtm(np.dot(rd_list[i+1].d,rd_list[i].d.T))),rd_list[i].d).T),np.array([0,0,1]))*180/np.pi for i in range(len(s)-1)]
        hps=[ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(0.0,0.0,0.0,0.0,0.0,0.0))]
        for i in range(len(s)-1):
            hps.append(ds.HP(ds.HP_intra(0.0,0.0,0.0,0.0,0.0,0.0),ds.HP_inter(ssr[i][0],ssr[i][1],ssr[i][2],Tilt[i],Roll[i],Twist[i])))
        hp_list=hps
    else:
        print('Please provide a valid type, "3DNA", "CURVES" or "MATH"')
        sys.exit(0)
    return ds.SC(HP=hp_list,RD=rd_list)

def atom_extract(Mask_3D, atom):
    """
    Given (a list of) 3D Mask, extract the coord's of needed atom.
    """
    xyz_list=[]
    if type(Mask_3D) is list:
        for i in Mask_3D:
            for x,y in enumerate(i.des):
                if y==atom:
                    if i.values.size==3:
                        i.values = i.values.reshape(1,3)
                    xyz_list.append(i.values[x])
    else:
        for x,y in enumerate(Mask_3D.des):
            if y==atom:
                if Mask_3D.values.size==3:
                    Mask_3D.values=Mask_3D.values.reshape(1,3)
                xyz_list.append(Mask_3D.values[x])
    return xyz_list

def RD2Mask_3D(RD):
    """
    Given a list of RD(global), convert it into a 3D Mask.
    """
    entry=RD[0]
    exit=RD[-1]
    values=np.array([x.r for x in RD])
    return ds.Mask_3D(values,RD_entry=entry,RD_exit=exit,skip=len(RD),des=['CA']*len(RD))

def docking_Mask_pdb(rd,pdb_list,entry=ds.RD(np.zeros(3),np.eye(3)), exit=ds.RD(np.zeros(3),np.eye(3)),skip=0):
    """
    A usage of docking_Mask_3D, but pdb has its own format.
    This function will input a list of pdb, docking the x y z values to the desired location,
    then output the pdb format list with all other information unchanged.
    """
    pdb_list_new=copy.deepcopy(pdb_list)
    l = [[i.x,i.y,i.z] for i in pdb_list_new]
    values = np.array(l,dtype=float)
    values_new = docking_Mask_3D(rd,ds.Mask_3D(values,RD_entry=entry,RD_exit=exit,skip=skip)).values
    if values_new.size==3:
        values_new=values_new.reshape(1,3)
    for i,j in enumerate(values_new):
        pdb_list_new[i].x = j[0]
        pdb_list_new[i].y = j[1]
        pdb_list_new[i].z = j[2]
    return pdb_list_new

def Mask_3D_stack(Mask_3D_1,Mask_3D_2):
    """
    Given two 3D Masks, the first one is global positioned, the latter one is local.
    calculate the latter one to the global position.
    """
    rd=Mask_3D_1.RD_exit
    new_Mask = docking_Mask_3D(rd, Mask_3D_2)
    return new_Mask

def DNA_allatom_pdb_combine(pdb_list):
    """
    Given a list of DNA pdb list, return a combined DNA pdb format data.
    Works only for double stranded DNA.
    """
    ###Identify the name for chain A and B
    chainA_name = pdb_list[0][0].chainID
    chainB_name = pdb_list[0][-1].chainID
    ###Separate chain A and B
    chainA=[]
    chainB=[]
    for i in pdb_list:
        chainA.append([x for x in i if x.chainID==chainA_name])
        chainB.append([x for x in i if x.chainID==chainB_name])
    ###Initialize serial and resSeq
    ser=1
    res=1
    ###Holder for combined DNA_pdb
    combined_pdb=[]
    ###Start combine chainA, tmp is for store the number of connection between O3' and P(next base pair)
    tmp=None
    for i in chainA:
        for j in i:
            j.serial=ser
            j.resSeq=res
            j.CONECT=[x-j.CONECT[0]+ser for x in j.CONECT]
            if j.name.split()[0]=="O3'":
                tmp=copy.copy(ser)
            if j.name.split()[0]=="P" and tmp is not None:
                j.CONECT.insert(1,tmp)
                combined_pdb[tmp-1].CONECT.append(ser)
            combined_pdb.append(j)
            ser+=1
        res+=1
    ###Start combine chainB, it is reverse.
    chainB.reverse()
    tmp=None
    for i in chainB:
        for j in i:
            j.serial=ser
            j.resSeq=res
            j.CONECT=[x-j.CONECT[0]+ser for x in j.CONECT]
            if j.name.split()[0]=="O3'":
                tmp=copy.copy(ser)
            if j.name.split()[0]=="P" and tmp is not None:
                j.CONECT.insert(1,tmp)
                combined_pdb[tmp-1].CONECT.append(ser)
            combined_pdb.append(j)
            ser+=1
        res+=1
    return combined_pdb

#########Plotting##########

def two_angle_plot(alpha,beta,filename):
    """
    Given list of alpha and beta and plot the two angle plot.
    """
    fig, ax = plt.subplots(figsize=(3,3))
    cm = plt.cm.get_cmap('RdYlBu')
    z = [float(x)/len(alpha) for x in range(len(alpha))]
    sc = plt.scatter(alpha,beta,c=z,cmap=cm)
    cbar = fig.colorbar(sc, ticks=[0, 1, 15])
    plt.xlim(0,np.pi)
    ax.set_xticks([0,np.pi/2,np.pi])
    ax.set_xticklabels(['0', '$\pi$/2', '$\pi$'], fontsize=12)
    ax.set_yticks([-np.pi, 0, np.pi])
    ax.set_yticklabels(['-$\pi$','0', '$\pi$'], fontsize=12)
    plt.ylim(-np.pi,np.pi)
    plt.xlabel(r'$\alpha$', fontsize=12)
    plt.ylabel(r'$\beta$', fontsize=12)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    fig.savefig(filename)
    plt.close()

def distance_matrix_plot(distance_matrix,filename,cut=0):
    """
    Input distance matrix, plot a heatmap of it.
    """
    fig, ax = plt.subplots(figsize=(3,3))
    sc = plt.imshow(distance_matrix,cmap='gist_heat',interpolation='nearest')
    cbar = fig.colorbar(sc, ticks=[0])
    if cut>0:
        cbar.ax.set_yticklabels(['<'+str(cut)])
    fig.savefig(filename)
    plt.close()
