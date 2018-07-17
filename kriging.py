import os, sys
import numpy as np
from dataset import Dataset
from variogram import Variogram
from grid import Grid
from grid import Grid
import subprocess
from subprocess import PIPE, STDOUT


class Kriging(object):
    def __init__(self, parfile = "krigpar.par", dataset=Dataset,
                 variogram = Variogram, grid = Grid, work_space=os.getcwd(),
                 exe_folder = ".\gslib77\BIN", nxdis = 1, nydis = 1, nzdis = 1,
                 min_points = 4, max_points = 8, max_per_octant = 0, max_search_raduis = None,
                 search_angles = [0, 0, 0], krig_type = 0, skmean = 0, trendfunction = None,
                 itrend = 0, trend_file = 'extdrift.dat' , GridDatacolumnNo = 1, debug_level = 3,
                 debug_file= 'kt3d.dbg', output_file = 'kr3d.out'):
        self.parfile = parfile
        self.work_space = work_space
        self.exe_folder = exe_folder
        self.dataset = dataset()
        self.grid = grid()
        self.variogram = Variogram
        self.triming_limits = [-1.0e21,   1.0e21]
        self.krig_option = 0 #
        self.options = {0:"grid", 1:"cross", 2:"jack"}
        self.jackknife_file = "xvk.dat"
        self.debug_level = debug_level
        self.debug_file = debug_file
        self.output_file = output_file
        self.nxdis = nxdis # these three variables equal means that point kriging is performed
        self.nydis = nydis
        self.nzdis = nzdis
        self.min_points = min_points # max and min number of points
        self.max_points = max_points
        self.max_per_octant = max_per_octant        #\max per octant (0-> not used)
        if max_search_raduis == None:
            rx = self.grid.lenX * 0.50
            ry = self.grid.lenY * 0.50
            rz = self.grid.lenZ * 0.50
            self.max_search_raduis = [rx, ry, rz]
        else:
            self.max_search_raduis = max_search_raduis
        self.search_angles = search_angles
        self.krig_type = krig_type     #2.302    #0=SK,1=OK,2=non-st SK,3=exdrift
        if skmean == None:
            try:
                self.skmean = np.nanmean(self.dataset.primary)
            except:
                self.skmean = None

        else:
            self.skmean = skmean # simple kriging mean
        if trendfunction == None:
            self.trendfunction = [0, 0, 0, 0, 0, 0, 0, 0, 0] #               \drift: x,y,z,xx,yy,zz,xy,xz,zy
        else:
            self.trendfunction = trendfunction
        self.itrend = itrend #\0, variable; 1, estimate trend
        self.trend_file = trend_file # \gridded file with drift/mean
        self.GridDatacolumnNo = GridDatacolumnNo

    def update_krige_par(self):
        rx = self.grid.lenX * 1.0
        ry = self.grid.lenY * 1.0
        rz = self.grid.lenZ * 1.0
        self.max_search_raduis = [rx, ry, rz]
        self.skmean = np.nanmean(self.dataset.primary)
        pass

    def __str__(self):
        return "To be done ..."

    def to_parfile(self):
        filename = os.path.join(self.work_space, self.parfile)
        self.update_krige_par()
        with open(filename, 'w') as x_file:
            # write header
            x_file.write("                  Parameters for KT3D\n")
            x_file.write("                  *******************\n")
            x_file.write("\n")
            x_file.write("START OF PARAMETERS:\n")

            # write data file
            dtfil = self.dataset.filename + '\n'
            x_file.write(dtfil)

            # columns number
            lin = "1 2 3 4"
            if not(self.dataset.secondary == None):
                lin = lin + " 5 \n"
            else:
                lin = lin + " 0 \n"
            x_file.write(lin)

            # limits
            lin = "%2.2e %2.2e\n" % (self.triming_limits[0],self.triming_limits[1])
            lin = lin.replace('+', '')
            x_file.write(lin)

            # write kriging options
            lin = str(self.krig_option) + '\n'
            x_file.write(lin)
            lin = self.jackknife_file + '\n'
            x_file.write(lin)

            # jacknife columns
            # columns number
            lin = "1 2 3 4"
            if not (self.dataset.secondary == None):
                lin = lin + " 5 \n"
            else:
                lin = lin + " 0 \n"
            x_file.write(lin)

            # debug level
            lin = "%d\n"%(self.debug_level)
            x_file.write(lin)

            # debug file
            lin = self.debug_file + '\n'
            x_file.write(lin)

            # output file
            lin = self.output_file + '\n'
            x_file.write(lin)


            # write grid
            grid = self.grid
            # x
            delx = float(grid.lenX)/grid.ncols
            lin = "%d %f %f \n"%(grid.ncols, grid.origin[0], delx)
            x_file.write(lin)

            #y
            dely = float(grid.lenY)/grid.nrows
            lin = "%d %f %f \n"%(grid.nrows, grid.origin[1], dely)
            x_file.write(lin)

            # z
            delz = float(grid.lenZ)/grid.nlays
            lin = "%d %f %f \n"%(grid.nlays, grid.origin[2], delz)
            x_file.write(lin)

            # x,y and z block discretization
            lin = "%d %d %d\n"%(self.nxdis, self.nydis, self.nzdis )
            x_file.write(lin)

            # min/ max points for kriging
            lin = "%d %d\n"%(self.min_points,self.max_points)
            x_file.write(lin)

            # max_per_octant
            lin = "%d\n"%(self.max_per_octant)
            x_file.write(lin)

            # max search distance: default is 50% of the max distance
            lin = "%f %f %f\n"%(self.max_search_raduis[0], self.max_search_raduis[1], self.max_search_raduis[2])
            x_file.write(lin)

            # angles of search
            lin = "%f %f %f\n"%(self.search_angles[0], self.search_angles[1], self.search_angles[2])
            x_file.write(lin)

            # kriging type
            lin = "%d %f\n"%(self.krig_type, self.skmean)
            x_file.write(lin)

            # trend
            for coeff in self.trendfunction:
                lin = "%d "%(coeff)
                x_file.write(lin)
            x_file.write('\n')

            # estimate trend flag
            lin = "%d \n"%(self.itrend)
            x_file.write(lin)

            # write trend file name
            x_file.write(self.trend_file)
            x_file.write('\n')

            # \  column number in gridded file
            lin = "%d \n"%(self.GridDatacolumnNo)
            x_file.write(lin)

            # write variograms
            lin = "%d %f\n"%(self.variogram.n_varios, self.variogram.nugget)
            x_file.write(lin)

            var_index = 1
            for ivar in np.arange(self.variogram.n_varios):
                variance = self.variogram.sills[var_index-1]
                ang = self.variogram.angles
                corr_lengths = self.variogram.cor_lengths
                lin = "%d %f %f %f %f\n"%(var_index, variance, ang[0], ang[1], ang[2])
                x_file.write(lin)
                lin = "         %f  %f  %f\n"%(corr_lengths[0], corr_lengths[1], corr_lengths[2])
                x_file.write(lin)
                var_index = var_index + 1

        x_file.close()


    def from_parfile(self):
        pass

    def krg_run(self):

        # write datafile
        self.dataset.write_file(self.ws)

        # write kriging parfile
        self.to_parfile()

        # run kriging
        FNULL = open(os.devnull, 'w')  # use this if you want to suppress output to stdout from the subprocess
        par_file = self.parfile
        par_file = os.path.join(self.work_space, par_file)
        exe = os.path.join(self.exe_folder, 'KT3D.exe' )
        args = par_file + '\n'

        #########################
        ris = subprocess.Popen(executable=exe, args="", stdin=PIPE,
                               universal_newlines=True, shell=True)
        com = ris.stdin.write(args)
        ris.stdin.flush()
        ris.wait()

        ########################




        pass








class Cokriging(Kriging):
    pass

class Ikriging(object):
    pass