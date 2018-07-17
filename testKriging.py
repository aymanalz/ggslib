import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from kriging import Kriging

kg = Kriging()

# input dataset
X = np.random.rand(100)
Y = np.random.rand(100)
Z = np.random.rand(100)
Z = Z/Z
primary_value = np.random.rand(100)
dts = kg.dataset
dts.filename = "data2.dat"
dts.x = X
dts.y = Y
dts.z = Z
plt.scatter(X,Y)
dts.primary = primary_value
kg.dataset = dts
#input variogram
vario = kg.variogram()
var1 = dict()
var = {'types':['Spherical'],
        'vario_contribs':[1.0],
        'cor_lengths': [0.6, 0.6, 0.6],
        'sills':[1.0],
        'angles':[0,0,0],
        'nuggets':[0.0]}
vario.var_dict = var
kg.variogram = vario

# input grid
gr = {'lenX':1, #
      'lenY':1,   # Lx
      'lenZ':1,
      'nlays':1,
      'ncols':50,
      'nrows':50,
      'origin':[0,0,0]}
kg.ws = "D:\Workspace\Codes\gslibpy\working_directory"

grid = kg.grid
grid.grid_dict = gr
kg.grid = grid

# run kriging
kg.krg_run()


pass






