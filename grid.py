__author__ = "Ayman H Alzraiee"
import os, sys
import numpy as np

class Grid(object):
    def __init__(self, ws = os.getcwd(), lenX=1.0, lenY = 1.0,
                 lenZ = 1.0, ncols = None, nrows = None, nlays = None,
                 origin = None):
        self.ws = ws
        self.lenX = lenX
        self.lenY = lenY
        self.lenZ = lenZ
        self.ncols = ncols
        self.nrows = nrows
        self.nlays = nlays
        self.origin = origin
        self._grid_dict = None
        self.del_xyz = None

    @property
    def grid_dict(self):
        gtemp = {'lenX': 1.0, 'lenY': 1.0, 'lenZ': 1.0,
              'nlays': 1, 'ncols': 50, 'nrows': 50,
              'origin': [0, 0, 0]}
        for key in gtemp.keys():
            gtemp[key] = getattr(self,key)
        self._grid_dict = gtemp
        return self._grid_dict


    @grid_dict.setter
    def grid_dict(self, gr):
        for key in gr.keys():
            setattr(self, key, gr[key])
        self._grid_dict = gr


    def plot(self):
        pass