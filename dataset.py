import os,sys
import numpy as np
import matplotlib.pyplot as plt


class Dataset(object):

    def __init__(self, cooment = None, x= None, y= None,
                 z = None, primary = None, secondary = None,
                 decluster_weight=None, ws = os.getcwd(), fname = "data01.dat"):

        self.comment = "Data input File"
        self.x = x
        self.y = y
        self.z = z
        self.primary = primary
        self.secondary = secondary
        self.decluster_weight = decluster_weight
        self.ws = ws
        self.fname = os.path.basename(fname)


    def write_file(self,ws):
        """
            This function write data to file
            :param data:
            :return:
        """
        self.check()
        file_name = os.path.join(os.path.abspath(ws), self.filename)
        #file_name = self.filename
        siz = self.dtable.shape

        with open(file_name, 'w') as x_file:
            if self.comment:
                x_file.write(self.comment)
                x_file.write("\n")

            stt = "%d\n"%(6)
            x_file.write(stt)

            x_file.write("Xlocation \n")
            x_file.write("Ylocation \n")
            x_file.write("Zlocation \n")
            x_file.write("Primary \n")
            x_file.write("Secondary \n")
            x_file.write("Decluster Weight \n")

            for row in np.arange(siz[0]):
                st=''
                for col in np.arange(siz[1]):
                    st = st + "%10.5f " % (self.dtable[row,col])
                st = st + "\n"
                x_file.write(st)

        x_file.close()

    def check(self):
        fields = ['x','y','z','primary','secondary','decluster_weight']
        dtable = np.array([])
        n_data = max([len(self.x) , len(self.y) , len(self.z)])

        for fd in fields :
            var = getattr(self, fd)
            if var is not None:
                self.var = var
                var = self.con_data()
                setattr(self, fd, var)
                if dtable.shape[0]==0:
                    dtable = np.reshape(var,(len(var),1))
                else:
                    var = np.reshape(var, (len(var), 1))
                    dtable = np.hstack((dtable,var))
            else:
                if fd == 'secondary':
                    dd = np.zeros((n_data,1))
                    dtable = np.hstack((dtable, dd))
                elif fd == 'decluster_weight':
                    dd = np.ones((n_data,1))/n_data
                    dtable = np.hstack((dtable, dd))





        delattr(self, "var")
        self.dtable = dtable

    def con_data(self):

        # make sure all variables are 1d numpy array
        var = self.var
        isarray = type(var) is np.ndarray
        isalist = isinstance(var, list)

        if isarray:
            var = var.ravel()
            return var
        elif isalist:
            var = np.array(var,dtype=float)
            return var
        else:
            raise ValueError('Data must be lists or numpy arrays')



    def data_plot(self):
        pass

    def from_datafile(self):
        pass

    def from_dic(self):
        pass

    @property
    def datafile(self):
        self.fname = os.path.basename(self.fname)
        self._datfile = os.path.join(self.ws,self.fname)
        return self._datfile

    @datafile.setter
    def datafile(self, dtfile):
        if not (os.path.basename(dtfile) == dtfile):
            self.fname = os.path.basename(dtfile)
            self.ws = os.path.dirname(dtfile)
        self._datfile = os.path.join(self.ws, self.fname)









