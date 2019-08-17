# -*- coding: utf-8 -*-

"""io (input/read, output/write) module."""

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from ds import ds as ds
from urllib.parse import urlparse

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








