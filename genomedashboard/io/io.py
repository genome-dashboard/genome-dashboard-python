# -*- coding: utf-8 -*-

"""io (input/read, output/write) module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
from urllib.parse import urlparse
import numpy as np
import twobitreader
import pyBigWig

def module_name():
    print("Module: genomedashboard.io.io.py.")

####################
### Check valid url
####################

def _CHECKURL(url):
    """Check whether the url is valid"""
    try:
        result = urlparse(url)
        return all([result.scheme, result.netloc])
    except:
        return False

####################
### Check valid file
####################

def _CHECKFILE(fp):
    """Check whether the file is valid"""
    if os.path.isfile(fp) == True:
        return True
    else:
        return False


class READ(object):
    """A class to read data into data structure from different sources."""
    def __init__(self,fp):
        if _CHECKFILE(fp)==True:
            self.fp = fp
        elif _CHECKURL(fp)==True:
            self.fp = fp
        else:
            print ('No such file, please provide a valid file.')
            sys.exit(0)

    def hp(self, hptype = '3DNA'):
        f = open(self.fp,'r')
        content = [x.rstrip('\n') for x in f]
        f.close()
        data = [x.split()[1:] for x in content[3:]]
        seq = ds.SEQ([x.split()[0] for x in content[3:]])
        hps = []
        for i in data:
            hp_intra_tmp = ds.HP_intra(float(i[0]),float(i[1]),float(i[2]),float(i[3]),float(i[4]),float(i[5]))
            hp_inter_tmp = ds.HP_inter(float(i[6]),float(i[7]),float(i[8]),float(i[9]),float(i[10]),float(i[11]))
            hp_tmp = ds.HP(hp_intra_tmp, hp_inter_tmp, hptype)
            hps.append(hp_tmp)
        sc = ds.SC(HP=hps,SEQ=seq)
        return sc

    def rd(self):
        """
        Read rd from xyz file, specifically for xyz only contains 'CA' and 'H1,H2,H3' in base pair resolution
        Read into SC class, as the rd for the space curve.
        """
        f = open(self.fp,'r')
        content = [x.rstrip('\n') for x in f]
        f.close()
        data = [x.split()[1:] for x in content[2:]]
        rd = []
        for i in range(int(len(data)/4)):
            rdtmp = ds.RD(np.array(data[i*4],dtype=float),np.array(data[(i*4+1):(i*4+4)],dtype=float)-np.array(data[i*4],dtype=float))
            rd.append(rdtmp)
        sc = ds.SC(RD=rd)
        return sc

    def xyz(self):
        """
        Read general xyz file, which contains:
        the first row as base-pair numbers
        the second row as comments
        data from the third row, with the first column as atom names, the 2-4 column as xyz coordinates
        Data will read into a 3D Mask.
        """
        f = open(self.fp,'r')
        content = [x.rstrip('\n') for x in f]
        f.close()
        data = [x.split()[1:] for x in content[2:]]
        atoms = [x.split()[0] for x in content[2:]]
        return ds.Mask_3D(np.array(data,dtype=float),des=atoms)


    def sequence_txt(self):
        """
        Given seqin.txt, either with one column of sequence or one/several rows of sequence
        Return sequence(string) with uppercase 'ACGT's
        """
        with open(self.fp) as f:
            seq=''.join(line.replace('\n', '') for line in f)
        seq = ds.SEQ(seq.upper())
        return seq

    def sequence_2bit(self,chromatin,start,end):
        """
        Given 2bit file, chromatin, start position, end position.
        Return sequence with uppercase 'ACGT's
        url is not work, need 2bit file
        """
        tbf = twobitreader.TwoBitFile(str(self.fp))
        seq = tbf[str(chromatin)][int(start):int(end)]
        seq = ds.SEQ(seq.upper())
        return seq

    def K(self):
        """
        e.g. MD-B.dat ; 6x6 or 12x12 stiffness matrix
        """
        f =  open(self.fp,'r')
        content = [x.rstrip('\n') for x in f]
        f.close()
        k = {}
        for i in range(16):
            step = [x.split()[0] for x in content[i*7+1:i*7+2]]
            k[step[0]] = np.array([x.split()[1:] for x in content[i*7+1:i*7+7]],dtype='float')
        return k

    def bigwig(self,chrom,start,end):
        """
        Read BigWig files, given chromosome number, start location and end location, return values
        """
        bw = pyBigWig.open(self.fp)
        contents = bw.values(chrom,start,end)
        return contents

    def PDB(self):
        """
        Read standard PDB file, return a list of PDB_std data structure that contains all the information in ATOM and CONECT
        """
        f = open(self.fp,'r')
        content = [x.rstrip('\n') for x in f]
        f.close()
        atoms = [x for x in content if x[0:4]=='ATOM']
        conn = [x for x in content if x[0:6]=='CONECT']
        conn_list=[]
        for i in conn:
            tmp=[]
            for j in range(int((len(i)-11)/5)):
                if i[j*5+6:j*5+11]!='     ':
                    tmp.append(int(i[j*5+6:j*5+11]))
            conn_list.append(tmp)
        pdb_list = []
        for i in atoms:
            tmppdb = ds.PDB_std(atom=i[0:6],serial=i[6:11],name=i[12:16],altLoc=i[16],resName=i[17:20],chainID=i[21],resSeq=i[22:26],iCode=i[26],x=i[30:38],y=i[38:46],z=i[46:54],occupancy=i[54:60],tempFactor=i[60:66],segID=i[72:76],element=i[76:78],charge=i[78:80])
            for j in conn_list:
                if int(tmppdb.serial)==j[0]:
                    tmppdb.CONECT = j
            pdb_list.append(tmppdb)
        return pdb_list

    def chrom(self):
        """
        Read chrom.bin file that provides chromosome length.
        Output the format that will feed to pyBigWig headers.
        """
        f = open(self.fp,'r')
        content = [x.rstrip('\n') for x in f]
        f.close()
        data = [x.split() for x in content]
        chrom_header=[(i[0],int(i[1])) for i in data]
        return chrom_header

class WRITE(object):
    """a class to write data into different format"""
    def __init__(self,fp):
        self.fp = fp

    def hp(self,HP,seq=None):
        """Given HP data structure, write into HP files"""
        if seq == None:
            seq = ['L-L']*len(HP)
        else:
            seq = seq.tolist()
        f = open(self.fp,'w')
        f.write(str(len(HP)) + ' base pairs' +'\n')
        f.write('   0  ***local base-pair & step parameters***' +'\n')
        f.write('       Shear  Stretch  Stagger Buckle Prop-Tw Opening   Shift  Slide    Rise    Tilt    Roll   Twist' + '\n')
        for i,j in enumerate(HP):
            f.write(seq[i]+' ')
            if j.HP_intra == None:
                f.write(str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8))
            else:
                f.write(str("%0.2f" % j.HP_intra.she).rjust(8)+str("%0.2f" % j.HP_intra.str).rjust(8)+str("%0.2f" % j.HP_intra.sta).rjust(8)+str("%0.2f" % j.HP_intra.buc).rjust(8)+str("%0.2f" % j.HP_intra.pro).rjust(8)+str("%0.2f" % j.HP_intra.ope).rjust(8))
            if j.HP_inter == None:
                f.write(str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8)+str("%0.2f" % 0.0).rjust(8))
            else:
                f.write(str("%0.2f" % j.HP_inter.shi).rjust(8)+str("%0.2f" % j.HP_inter.sli).rjust(8)+str("%0.2f" % j.HP_inter.ris).rjust(8)+str("%0.2f" % j.HP_inter.til).rjust(8)+str("%0.2f" % j.HP_inter.rol).rjust(8)+str("%0.2f" % j.HP_inter.twi).rjust(8))
            f.write('\n')
        f.close()

    def rd(self,RD):
        """Given RD data structure, write into xyz files, this is for dna in base pair level"""
        f = open(self.fp,'w')
        f.write(str(len(RD)*4)+'\n')
        f.write('COMMENT: zli' + '\n')
        for i in RD:
            f.write('CA'.ljust(7)+str("%0.5f" % i.r[0]).rjust(20)+str("%0.5f" % i.r[1]).rjust(20)+str("%0.5f" % i.r[2]).rjust(20)+'\n')
            f.write('H1'.ljust(7)+str("%0.5f" % (i.r+i.d[0])[0]).rjust(20)+str("%0.5f" % (i.r+i.d[0])[1]).rjust(20)+str("%0.5f" % (i.r+i.d[0])[2]).rjust(20)+'\n')
            f.write('H2'.ljust(7)+str("%0.5f" % (i.r+i.d[1])[0]).rjust(20)+str("%0.5f" % (i.r+i.d[1])[1]).rjust(20)+str("%0.5f" % (i.r+i.d[1])[2]).rjust(20)+'\n')
            f.write('H3'.ljust(7)+str("%0.5f" % (i.r+i.d[2])[0]).rjust(20)+str("%0.5f" % (i.r+i.d[2])[1]).rjust(20)+str("%0.5f" % (i.r+i.d[2])[2]).rjust(20)+'\n')
        f.close()

    def xyz(self,Masks_3D):
        """
        xyz is for 3D Masks, not limited by DNA, but it does not specific direction frames.
        Given a list of Masks, the atom/element name should provided in Mask_3D.des,
        if there is no des, ATOM will be the name,
        des should be a list, e.g. ['O','CA',...]
        """
        content=[]
        for i in Masks_3D:
            if i.values.size==3:
                i.values=i.values.reshape(1,3)
            for x,j in enumerate(i.values):
                if i.des is None:
                    atom='ATOM'
                else:
                    atom=i.des[x]
                content.append(atom.ljust(7)+str("%0.5f" % j[0]).rjust(20)+str("%0.5f" % j[1]).rjust(20)+str("%0.5f" % j[2]).rjust(20)+'\n')
        f=open(self.fp,'w')
        f.write(str(len(content))+'\n')
        f.write('COMMENT: zli' + '\n')
        for i in content:
            f.write(i)
        f.close()

    def pdb(self,pdb_list):
        """
        Given a list of pdb(class PDB_std), write into a standard pdb file.
        Current version supports ATOM and CONECT
        """
        f=open(self.fp,'w')
        for i,j in enumerate(pdb_list):
            if i>0:
                if j.chainID != pdb_list[i-1].chainID:
                    f.write('TER'+'\n')
            f.write(j.atom.ljust(6)+str(j.serial).rjust(5)+' '+j.name.ljust(4)+j.altLoc+j.resName.rjust(3)+' '+j.chainID+str(j.resSeq).rjust(4)+j.iCode+' '.rjust(3)+str("%0.3f" % float(j.x)).rjust(8)+str("%0.3f" % float(j.y)).rjust(8)+str("%0.3f" % float(j.z)).rjust(8)+j.occupancy.rjust(6)+j.tempFactor.rjust(6)+' '.rjust(6)+j.segID.ljust(4)+j.element.ljust(2)+j.charge.ljust(2))
            f.write('\n')
        f.write('TER'+'\n')
        for k in pdb_list:
            if k.CONECT is not None:
                f.write('CONECT')
                for x in k.CONECT:
                    f.write(str(x).rjust(5))
                f.write('\n')
        f.write('END')
        f.close()

    def pdb2xyz(self,pdb_list):
        """
        Given a list of pdb(class PDB_std), write into xyz file.
        """
        f=open(self.fp,'w')
        f.write(str(len(pdb_list))+'\n')
        f.write('COMMENT: zli' + '\n')
        for i in pdb_list:
            f.write(i.name.strip()[0].ljust(7)+str("%0.5f" % float(i.x)).rjust(20)+str("%0.5f" % float(i.y)).rjust(20)+str("%0.5f" % float(i.z)).rjust(20)+'\n')
        f.close()

    def bigwig(self,chrom_len,chrom,start,end,values):
        """
        Write bigwig format file.
        chrom_len can be obtained by chrom.bin file from UCSC
        """
        bw = pyBigWig.open(self.fp,'w')
        bw.addHeader(chrom_len,maxZooms=0)
        bw.addEntries(np.array([chrom]*len(start)),start,ends=end,values=values)
        bw.close()

    def sequence_txt(self,SEQ):
        """
        Given a (class SEQ) sequence, output the txt sequence file.
        """
        f=open(self.fp,'w')
        for i in SEQ.tostring():
            f.write(i+'\n')
        f.close()
