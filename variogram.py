import numpy as np

class Variogram(object):

    def __init__(self, n_varios = 1, types = ["spherical"], vario_contribs = [1],
                 cor_lengths = None, sills = None, angles = np.array([0,0,0]),
                 nugget = 0):
        """
        Azimuth angle :  is the angle measured clockwise from north-south
        dip angle : measured in negative degress from the horizontal plan.
        """
        self.n_varios = n_varios
        self.types = types # a list of variograms types
        self.options = {'Spherical':1,'Exponential':2, \
                        'Gaussian':3, 'Power':4, 'Hole Effect':5}
        self.vario_contribs = vario_contribs
        self.cor_lengths = cor_lengths
        self.sills = sills
        self.angles = angles
        self.nugget = nugget
        self._var_dict = None

    @property
    def var_dict(self):
        vr_dic = {'types':['Spherical'],
        'vario_contribs':[1.0],
        'cor_lengths': [0.5],
        'sills':[1.0],
        'angles':[0,0,0],
        'nugget':0.0}

        for key in vr_dic.keys():
            vr_dic[key] = getattr(self,key)
        self._var_dict = vr_dic
        return self._var_dict


    @var_dict.setter
    def var_dict(self,vr_dic):
        for key in vr_dic.keys():
            setattr(self,key,vr_dic[key])
        self._var_dict = vr_dic
        return self._var_dict




class Emperical_variogram(object):

    def __init__(self, parfile = "gamv.par", dataset = None, nvar = 1, var_cols = None,
                 limits = np.array([-1e21, 1e21]), outfile = "gamv.out", nlags = 10,
                 lag_sep_dist=None, lag_tolerance = 0.5, Num_directions = 1,
                 var_geom = {'azm': [0],'atol':[20],'bandwh':[1.0], 'dip' : [0.0], 'dtol':[1], 'bandwd':[1]},
                 varhead_tail={'ivtail': [1], 'ivhead': [1], 'ivtype': [1]}):

        self.parfile = parfile
        self.dataset = None
        self.nvar = 1
        self.var_cols = None
        self.limits = limits
        self.outfile = outfile
        self.nlags = nlags
        self.lag_sep_dist = lag_sep_dist
        self.lag_tolerance = lag_tolerance #lag
        self.Num_directions = 1
        # for each directions we should define the geometry of the search
        self.var_geom = var_geom
        # for each variogram define head, tail, and variogram type
        self.varhead_tail = varhead_tail
        self.standardize_sills = True





    def calc_exp_variogram(self):
        pass

    def write_exp_vario_file(self):
        fn = self.filename
        with open(fn, 'w') as file:
            file.write("START OF PARAMETERS:\n")
            # data file
            stt = self.datafile + " \n"
            file.write(stt)

            # Columns for x, y, z

            stt = self.calc_exp_variogram()





    def fit_variogram(self):
        x = 1
