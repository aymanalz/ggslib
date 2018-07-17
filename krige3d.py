import os, sys
import numpy as np
import pandas as pd
from dataset import Dataset
from variogram import Variogram
from grid import Grid
import subprocess
from subprocess import PIPE, STDOUT
import matplotlib.pyplot as plt

class Kriging(object):
    def __init__(self, ws = os.getcwd(), par_file='krige.par', grid = Grid,
                 dataset=Dataset,
                 variogram = Variogram):


        self.ws = ws
        self.par_file = par_file
        self.__read_template()
        self.grid = grid
        self.dataset=dataset
        self.variogram = variogram
        self.exe_folder = ".\gslib77\BIN"


    def write_file(self):
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

        # update dataset
        self.data_file = os.path.join(self.ws, self.dataset.filename)


        # update grid
        self.update_grid()

        # update variograms
        self.update_variogram()


    def update_variogram(self):

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
        self.dataset.write_file(self.ws)
        self.write_file()

        # run kriging
        FNULL = open(os.devnull, 'w')  # use this if you want to suppress output to stdout from the subprocess
        par_file = self.par_file
        par_file = os.path.join(self.ws, par_file)
        exeF = os.path.join(os.path.dirname(__file__),r'gslib77\BIN')
        exe = os.path.join(exeF, 'KT3D.exe' )
        args = par_file + '\n'
        ris = subprocess.Popen(executable=exe, args="", stdin=PIPE,
                               universal_newlines=True, shell=True)

        com = ris.stdin.write(args)
        ris.stdin.flush()
        ris.wait()

        outfile = self.kriged_output_file
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
        self.output['data'] = np.array(outdata)
        data = self.output['data']
        self.krigr_2d = np.reshape(data[:,0], (self.ny, self.nx))
        self.krigr_2d = self.krigr_2d.astype('float')
        self.map_variance = np.reshape(data[:, 1], (self.ny, self.nx))
        self.map_variance = self.map_variance.astype('float')
        self.output['columns'] = columns


    def __read_template(self):
        cwd = os.path.dirname(os.path.realpath(__file__))
        temp_file = temp_file = os.path.join(cwd,r"templates\KT3D.par" )
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

    def plot(self, nodata = -999):

        idx = 0
        data = self.output['data']
        for col in self.output['columns']:
            plt.subplot(1,2,idx+1)
            curr_ = data[:,idx]
            curr_ = curr_.astype('float')
            curr_ = np.reshape(curr_,( self.ny, self.nx))
            #curr_ = np.flipud(curr_)
            x_x = np.linspace(self.xmn,self.xmn+self.grid.grid_dict['lenX'], self.nx)
            y_y = np.linspace(self.ymn, self.ymn+self.grid.grid_dict['lenY'], self.ny)
            xx, yy = np.meshgrid(x_x,y_y)
            loc = np.where(curr_ == nodata)
            curr_[loc] = np.nan
            plt.pcolormesh(xx,yy,curr_,cmap='jet')
            plt.colorbar()

            idx = idx + 1







if __name__ == "__main__":
    kg = Kriging(ws = ".\working_directory")

    # input dataset
    X = np.random.rand(50)*50.0
    Y = np.random.rand(50)*50.0
    Z = np.random.rand(50)
    Z = Z / Z
    primary_value = np.random.rand(50)*3.0
    dts = kg.dataset()
    dts.filename = "data1.dat"
    dts.x = X
    dts.y = Y
    dts.z = Z
    plt.scatter(X, Y)
    dts.primary = primary_value
    kg.dataset = dts

    # input variogram
    vario = kg.variogram()
    var1 = dict()
    var = {'types': ['Spherical'],
           'vario_contribs': [1.0],
           'cor_lengths': [30.0, 30.0, 30.0],
           'sills': [3],
           'angles': [0, 0, 0],
           'nuggets': [0]}
    vario.var_dict = var
    kg.variogram = vario
    kg.max_search_radii_1 = 200
    kg.max_search_radii_2 = 200
    kg.max_search_radii_3 = 200
    kg.max_data_points = 20
    kg.min_data_points = 4

    # input grid
    gr = {'lenX': 50,  #
          'lenY': 50,  # Lx
          'lenZ': 1,
          'nlays': 1,
          'ncols': 100,
          'nrows': 100,
          'origin': [0, 0, 0]}

    grid = kg.grid
    grid.grid_dict = gr
    kg.grid = grid
    kg.krige_type = 1 # ordinary
    kg.run()

    kg.plot(nodata=-999)
    plt.hold(True)
    plt.scatter(X, Y,color = 'k')
    plt.show()

    x = 1

    pass