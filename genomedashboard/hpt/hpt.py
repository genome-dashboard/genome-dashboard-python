# -*- coding: utf-8 -*-

"""hpt (helix parameter tools) module."""

def hello():
    print("Hello from genomedashboard.hpt.hpt.py.")


class HP(object):
    def __init__(self, nbp, hps, comments=None):
        self.nbp = nbp
        self.hps = hps
        if self.comments==None:
            self.comments = '   0  ***local base-pair & step parameters***'
