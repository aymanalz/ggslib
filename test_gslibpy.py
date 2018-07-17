from gslibpy import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# create gslib project object
gsproj = Gsbase()

# assigne working directory
gsproj.work_space="D:\Workspace\Codes\gslibpy\working_directory"

# Dataset
dts = gsproj.dataset
X = np.random.rand(10)
Y = np.random.rand(10)
Z = np.random.rand(10)
primary_value = np.random.rand(10)
dts.x = X
dts.y = Y
dts.z = Z
dts.primary = primary_value

pass




pass