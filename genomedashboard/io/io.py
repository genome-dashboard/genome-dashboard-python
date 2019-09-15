# -*- coding: utf-8 -*-

"""io (input/read, output/write) module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
from urllib.parse import urlparse
import numpy as np

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

    def sequence_txt(self):
        """
        Given seqin.txt, either with one column of sequence or one/several rows of sequence
        Return sequence(string) with uppercase 'ACGT's
        """
        with open(self.fp) as f:
            seq=''.join(line.replace('\n', '') for line in f)
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

