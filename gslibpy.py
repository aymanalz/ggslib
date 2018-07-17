import os,sys
import numpy as np
from dataset import dataset
from grid import grid


class Gsbase(object):
    def __init__(self, name=None, work_space = os.getcwd(), exe_folder = ".\gslib77\BIN",
                 ):
        self.name = name
        self._work_space = work_space
        self.exe_folder = exe_folder
        self.dataset = dataset()
        self.grid = grid()


    @property
    def work_space(self):
        return self._work_space

    @work_space.setter
    def work_space(self, value):
        self._work_space = value
        self.dataset.ws = value

    @work_space.deleter
    def work_space(self):
        pass


