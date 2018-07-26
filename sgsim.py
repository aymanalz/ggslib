_Auther__ = "Ayman H. Alzraiee"
__Vesrion__ = 1.0
__Date__ = "1-5-2015"

"""
A class to generate GSLIB input files , run GSLIB, and read output.
This class inherent objects from Grid, Variogram, and Dataset classes  
"""


 #TODO: use multiple variograms
 #TODO: More inline Documentation
 #TODO: Check the possibility of using fortran libaries for GSFLIB
 # TODO: Better plotting


import os, sys
import numpy as np
import pandas as pd
from dataset import Dataset
from variogram import Variogram
from grid import Grid
import subprocess
from subprocess import PIPE, STDOUT
import matplotlib.pyplot as plt

class Sgsim(object):
    """

    Sequential Gaussian simulation:
    Help URL: http://www.gslib.com/gslib_help/sgsim.html

    Sgsim is class that uses database class (hard data file), grid class, and variogram class.
    Example: see the ___main__ function at the end of this file

    For help about any object members use [print(objct.MembaerName_), for example to find information about
    sg.sim, you can print sg.nsim_
    """
    def __init__(self, ws = os.getcwd(), par_file='sgsim.par', grid = Grid,
                 dataset=Dataset,
                 variogram = Variogram):
        """

        :param ws: working space
        :param par_file: Parameter file name used by GSLIB, the default is sgsim.par
        :param grid: An object for a grid class
        :param dataset: An object for observation data
        :param variogram: An object from variogram class
        """
        self.ws = ws
        self.par_file = par_file
        self.__read_template()
        self.grid = grid
        self.dataset=dataset
        self.variogram = variogram
        self.exe_folder = ".\gslib77\BIN"


    def write_file(self):
        """
        Write the parameter file, which is the input for sgsim execuatble

        :return: None
        """
        self.update_file_content()
        fname = os.path.join(self.ws, self.par_file)
        fid = open(fname,'w')
        for record in self.file_content_fields:
            for par in record:
                val = getattr(self, par)
                if isinstance(val, str):
                    lin = str(val.strip())+ ' '
                else:
                    lin = str(val) + ' '
                fid.write(lin)
            fid.write('\n')
        fid.close()

    def update_file_content(self):

        """
        Before runing sgsim, this function updates all class members with new values if there is any

        :return:
        """
         # update dataset
        if hasattr(self.dataset, 'filename'):
            self.data_file = os.path.join(self.ws, self.dataset.filename)
        else: # unconditional simulation
            self.data_file = 'nodata.dat'
        self.outfl = os.path.join(self.ws, self.outfl)

        # update grid
        self.update_grid()

        # update variograms
        self.update_variogram()

    def update_variogram(self):
        """


        """
        # we support only one variogram
        vario = self.variogram
        self.n_structures = vario.n_varios
        self.nugget_effect = vario.nugget
        self.it = vario.options[vario.types[0]]
        self.cc = vario.sills[0]
        self.ang1 = vario.angles[0]
        self.ang2 = vario.angles[1]
        self.ang3 = vario.angles[2]
        self.a_hmax = vario.cor_lengths[0]
        self.a_hmin = vario.cor_lengths[1]
        self.a_vert = vario.cor_lengths[2]

    def update_grid(self):
        grid = self.grid
        self.nx = grid.grid_dict['ncols']
        self.ny = grid.grid_dict['nrows']
        self.nz = grid.grid_dict['nlays']
        self.xmn = grid.grid_dict['origin'][0]
        self.ymn = grid.grid_dict['origin'][1]
        self.zmn = grid.grid_dict['origin'][2]
        self.xsiz = grid.grid_dict['lenX']/float(self.nx)
        self.ysiz = grid.grid_dict['lenY']/float(self.ny)
        self.zsiz = grid.grid_dict['lenZ']/float(self.nz)


    def run(self):

        # write datafile
        if hasattr(self.dataset, 'primary'):
          self.dataset.write_file(self.ws)
        self.write_file()

        # run simulation
        FNULL = open(os.devnull, 'w')  # use this if you want to suppress output to stdout from the subprocess
        par_file = self.par_file
        par_file = os.path.join(self.ws, par_file)

        exeF = os.path.join(os.path.dirname(__file__), r'gslib77\BIN')
        exe = os.path.join(exeF, 'SGSIM.exe')
        args = par_file + '\n'
        ris = subprocess.Popen(executable=exe, args="", stdin=PIPE,
                               universal_newlines=True, shell=True)

        com = ris.stdin.write(args)
        ris.stdin.flush()
        ris.wait()

        outfile = self.outfl
        fid = open(outfile, 'r')
        content = fid.readlines()
        fid.close()
        idx = 0
        columns = list()
        outdata = list()
        for record in content:
            if idx == 0:  # header
                pass
            elif idx == 1:
                num_ = int(record.strip())
            elif (idx > 1) and (idx <= 1+num_):
                columns.append(record.strip())
            else:
                curr_ = record.strip().split()
                outdata.append(curr_)
            idx = idx + 1

        self.output = dict()
        outdata = np.array(outdata, dtype="float")
        ncells = self.nx * self.ny * self.nz
        outdata = np.reshape(outdata, (self.nsim, ncells))
        outdata = outdata.transpose()
        self.output['data'] = outdata
        self.output['columns'] = columns

    def __read_template(self):

        cwd = os.path.dirname(os.path.realpath(__file__))
        temp_file = temp_file = os.path.join(cwd, r".\templates\SGSIM.par")
        fid = open(temp_file,'r')
        contents = fid.readlines()
        fid.close()
        file_content = list()
        file_content_fields = list()
        flg = 0
        ii =0
        for rec in contents:
            if rec == 'START OF PARAMETERS:\n':
                flg = 1
                file_content.append(rec)
                st = 'Header' + str(ii)
                file_content_fields.append([st])
                ii = ii + 1
                setattr(self, st, rec)
                continue

            elif flg:
                parts = rec.split("\\")
                values = parts[0].split()
                file_content.append(values)
                names  = parts[1]
                fields_help = names.split("#")
                fields = fields_help[0].split(",")
                fields = [fd.strip() for fd in fields]
                file_content_fields.append(fields)
                helps = fields_help[1]

                idx = 0
                for fd in fields:
                    field_name = fd
                    value = values[idx]
                    try:
                        value = int(value)
                    except:
                        try:
                            value = float(value)
                        except:
                            pass

                    setattr(self, field_name, value)
                    fhelp = field_name + "_"
                    setattr(self, fhelp, helps)
                    idx = idx + 1
            else:
                st = 'Header' + str(ii)
                file_content_fields.append([st])
                ii = ii + 1
                setattr(self, st, rec)

        setattr(self, 'file_content', file_content)
        setattr(self, 'file_content_fields', file_content_fields)

    def plot(self, nodata = -999, realn = 1):

        idx = 0
        data = self.output['data']
        for col in self.output['columns']:
            plt.subplot(1,1,idx+1)
            curr_ = data[:,realn]
            curr_ = curr_.astype('float')
            curr_ = np.reshape(curr_,(self.ny, self.nx))
            x_x = np.linspace(self.xmn,self.grid.grid_dict['lenX'], self.nx)
            y_y = np.linspace(self.ymn, self.grid.grid_dict['lenY'], self.ny)
            xx, yy = np.meshgrid(x_x,y_y)
            loc = np.where(curr_ == nodata)
            curr_[loc] = np.nan
            plt.pcolormesh(xx,yy,curr_,cmap='jet')
            plt.colorbar()
            idx = idx + 1

if __name__ == "__main__":

    # Initialize Sgsim object
    sg = Sgsim(ws = ".\working_directory")

    # grid to be used. you can use diectionary for grid info.
    gr = {'lenX': 50,
          'lenY': 100,
          'lenZ': 1,
          'nlays': 1,
          'ncols': 237,
          'nrows': 134,
          'origin': [0, 0, 0]}

    grid = sg.grid
    grid.grid_dict = gr

    # Generate a random set of observations locations and values
    X = np.random.rand(50)*50.0 # x-coord
    X[0] = X[0] + 60
    Y = np.random.rand(50)*50.0 # y-coord
    Z = np.random.rand(50) # z-coord, for 2D make z constant
    Z = Z / Z
    primary_value = np.random.rand(50)*3.0 # values of measured parameters to be simulated (e.g. Coductivity)
    primary_value[0] = primary_value[0] * 1000

    # Assign artificail data to dataset class, if uncoditional simulation is used this will not used
    dts = sg.dataset()
    dts.filename = "data1.dat"
    dts.x = X
    dts.y = Y
    dts.z = Z
    dts.primary = primary_value
    #sg.dataset = dts    # uncoditional simulation --- dataset is not used

    # input variogram
    vario = sg.variogram()
    var1 = dict()
    var = {'types': ['Spherical'],
           'vario_contribs': [1.0],
           'cor_lengths': [30.0, 100.0, 30.0],
           'sills': [3],
           'angles': [45, 0, 0],
           'nuggets': [0]}
    vario.var_dict = var
    sg.variogram = vario
    sg.max_search_radii_1 = 200
    sg.max_search_radii_2 = 200
    sg.max_search_radii_3 = 200
    sg.max_data_points = 20
    sg.min_data_points = 4

    #sg.grid = grid
    sg.krige_type = 1 # ordinary
    sg.nsim = 10
    sg.run()
    sg.plot(nodata=-999, realn = 2)
    plt.show()
