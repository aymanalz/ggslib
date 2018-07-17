
import gslibpy as gs
import dataset as ds
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

X = np.random.rand(10)
Y = np.random.rand(10)
Z = np.random.rand(10)

primary_value = np.random.rand(10)
dts = ds.dataset()
dts.filename = "data.dat"
dts.x = X
dts.y = Y
dts.z = Z
dts.primary = primary_value

dts.write_file()
x = 1