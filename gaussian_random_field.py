
__author__ = "Ayman H. Alzraiee"
__Date__ = "04/2/2018"
__version__ = 2.0

"""
Pure python class for geostatistics including 3D kriging and 3D Guassian Random fields.
Features:    
    
To Do:
    * Cokriging
    * Multi-point Geostatistics 
    * Indicator Kriging
    * Indicator simulations
    * Unstructured Grids
     
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy import linalg
import dask

def mp(a, b):
    if b.ndim == 1:
        return linalg.blas.dgemv(alpha=1.0, a=a, x=b)
    else:
        return linalg.blas.dgemm(alpha=1.0, a=a, b=b)

class Covariance(object):
    """"
    Currently supports only exponential covariance....
    """
    def __init__(self, type = "exponential", scale = 0.5, sigma = 1.0):
        self.options = ['exponential']
        self.type = type
        self.scale = scale
        self.sigma = sigma

    def calc_covariance(self, r):
        if self.type == "exponential":
            cov = self.sigma * np.exp((-r/self.scale))
        return cov

    def cal_sparse_covariance(self, r, max_dist):
        """
        Assuming that locations at large distances (r > max_dist) are probably uncorrelated, it is possible to
        sparsify the covaraince matrix allowing faster calculations.
        :param r: a distance matrix
               max_dist: is
        :return: cov sparse covriance matrix
        """
        pass


class Gaussian_Field(Covariance):
    def __init__(self, data = [], sample_locations = [], block_raduis = None, search_radius = None,
                 max_points = None):
        Covariance.__init__(self)
        self.data = data
        self.sample_locations = sample_locations
        self.nrealizations = 1
        self.block_radius = block_raduis
        self.search_radius = search_radius
        self.max_points = max_points

    def decompose_covariance(self, covm):
        U, s, V = linalg.svd(covm, full_matrices=True)
        nu = U.shape[0]
        nv = V.shape[0]
        s = np.power(s, 0.5)
        if nu == nv:
            s = np.diag(s)
        else:
            s = np.vstack( (np.diag(s), np.zeros((nu-nv, nv))))

        US = linalg.blas.dgemm(alpha=1.0, a=U, b=s)
        s = None
        U = None
        L = linalg.blas.dgemm(alpha=1.0, a=US, b=V)
        return L

    def seq_simulation(self):
        nd = np.shape(self.data)[0]
        ns = np.shape(self.sample_locations)[0]
        sim_flg = np.zeros(ns)
        sim_ensemble = np.zeros((ns, self.nrealizations))
        while True:
            unsample_loc = np.where(sim_flg == 0)[0]
            if not len(unsample_loc)> 0:
                break
            center_loc = np.random.choice(unsample_loc)
            x0, y0 = self.sample_locations[center_loc,:]
            dx2 = np.power(self.sample_locations[:,0] - x0, 2.0)
            dy2 = np.power(self.sample_locations[:,1] - y0, 2.0)
            rad = np.power(dx2+dy2, 0.5)
            # find nodes around (x0,y0) that will be estimated
            nodes_block = np.where(np.logical_and(rad<= self.block_radius, sim_flg==0))[0]
            if len(nodes_block) > self.max_points:
                pass
            # get data within search radius
            nodes_data = np.where(np.logical_and(rad<=self.search_radius,
                                    sim_flg==1))[0]
            curr_data = np.hstack((self.sample_locations[nodes_data,:],sim_ensemble[nodes_data,:]))
            curr_unsampled_loc = self.sample_locations[nodes_block,:]
            sim = self.gen_real(curr_data, curr_unsampled_loc)
            sim_ensemble[nodes_block,:] = sim
            plt.imshow(sim_ensemble[:, 0].reshape(60, 60))
            sim_flg[nodes_block] = 1
            xx = 1






    def gen_real(self, data, unsamp_loc):
        nd = np.shape(data)[0]
        ns = np.shape(unsamp_loc)[0]
        if len(data) > 0:  # conditional sim
            data_loc = data[:, 0:-1]
            try:
                all_locations = np.vstack((data_loc, unsamp_loc))
            except:
                pass
        else:
            all_locations = unsamp_loc
        r = cdist(all_locations, all_locations)
        cov_mat = self.calc_covariance(r)
        if len(data)>0:
            L11 = self.decompose_covariance(cov_mat[0:nd, 0:nd])
            L21 = self.decompose_covariance(cov_mat[:, 0:nd])
            L11in = np.linalg.pinv(L11)
            L11 = None
            r = None
            y1 = data[:, 2]
            mean = mp(mp(L21, L11in), y1)
            L21 = None
            L11in = None
            mean = mean[:, np.newaxis]
        else:
            mean = 0.0

        L22 = self.decompose_covariance(cov_mat)
        cov_mat = None
        z = np.random.randn(ns + nd, int(self.nrealizations))
        fluct = mp(L22, z)
        gs_field = mean + fluct
        gs_field = gs_field[nd:]
        return gs_field
    def krige_1(self):
        pass

    def krige_seq(self):
        nd = np.shape(self.data)[0]
        ns = np.shape(self.sample_locations)[0]
        sim_flg = np.zeros(ns)
        sim_ensemble = np.zeros((ns, 1))
        while True:
            unsample_loc = np.where(sim_flg == 0)[0]
            if not len(unsample_loc) > 0:
                break
            center_loc = np.random.choice(unsample_loc)
            x0, y0 = self.sample_locations[center_loc, :]
            dx2 = np.power(self.sample_locations[:, 0] - x0, 2.0)
            dy2 = np.power(self.sample_locations[:, 1] - y0, 2.0)
            rad = np.power(dx2 + dy2, 0.5)
            # find nodes around (x0,y0) that will be estimated
            nodes_block = np.where(np.logical_and(rad <= self.block_radius, sim_flg == 0))[0]
            if len(nodes_block) > self.max_points:
                pass
            # get data within search radius
            nodes_data = np.where(np.logical_and(rad <= self.search_radius,
                                                 sim_flg == 1))[0]
            curr_data = np.hstack((self.sample_locations[nodes_data, :], sim_ensemble[nodes_data, :]))
            curr_unsampled_loc = self.sample_locations[nodes_block, :]
            sim = self.krige_block(curr_data, curr_unsampled_loc)
            sim_ensemble[nodes_block, :] = sim
            plt.imshow(sim_ensemble[:, 0].reshape(60, 60))
            sim_flg[nodes_block] = 1
            xx = 1

    def krige_block(self, data, unsamp_loc):
        nd = np.shape(data)[0]
        ns = np.shape(unsamp_loc)[0]
        if len(data) > 0:  # conditional sim
            data_loc = data[:, 0:-1]
            try:
                all_locations = np.vstack((data_loc, unsamp_loc))
            except:
                pass
        else:
            all_locations = unsamp_loc
        r = cdist(all_locations, all_locations)
        cov_mat = self.calc_covariance(r)
        if len(data) > 0:
            L11 = self.decompose_covariance(cov_mat[0:nd, 0:nd])
            L21 = self.decompose_covariance(cov_mat[:, 0:nd])
            L11in = np.linalg.pinv(L11)
            L11 = None
            r = None
            y1 = data[:, 2]
            mean = mp(mp(L21, L11in), y1)
            L21 = None
            L11in = None
            mean = mean[:, np.newaxis]
        else:
            mean = 0.0
        return mean


    def conditional_simulation(self):
        nd = np.shape(self.data)[0]
        ns = np.shape(self.sample_locations)[0]
        if len(self.data)>0: # conditional sim
            data_loc = self.data[:,0:-1]
            all_locations = np.vstack((data_loc, self.sample_locations))
        else:
            all_locations = self.sample_locations
        r = cdist(all_locations, all_locations)
        cov_mat = self.calc_covariance(r)
        L11 = self.decompose_covariance(cov_mat[0:nd, 0:nd])
        L22 = self.decompose_covariance(cov_mat[nd:, nd:])
        L21 = self.decompose_covariance(cov_mat[nd:, 0:nd])
        L11in = np.linalg.pinv(L11)
        L11 = None
        r = None; cov_mat = None
        y1 = self.data[:,2]
        mean = mp(mp(L21, L11in), y1)
        L21 = None
        L11in = None
        z = np.random.randn(ns, int(self.nrealizations))
        fluct = mp(L22, z)
        mean = mean[:,np.newaxis]
        gs_field = mean + fluct
        return gs_field

    def unconditional_simulation(self):
        nd = np.shape(self.data)[0]
        ns = np.shape(self.sample_locations)[0]
        if len(self.data)>0: # conditional sim
            data_loc = self.data[:,0:-1]
            all_locations = np.vstack((data_loc, self.sample_locations))
        else:
            all_locations = self.sample_locations
        r = cdist(all_locations, all_locations)
        cov_mat = self.calc_covariance(r)
        L = self.decompose_covariance(cov_mat)
        r = None; cov_mat = None
        if self.nrealizations > 1:
            z = np.random.randn(ns,int(self.nrealizations))
            sim_reralization =  linalg.blas.dgemm(alpha = 1.0, a = L, b = z )
        else:
            z = np.random.randn(ns, int(self.nrealizations))
            sim_reralization = linalg.blas.dgemv(alpha=1.0, a=L, x=z)
        return sim_reralization



if __name__ == "__main__":
    gf = Gaussian_Field()
    gf.scale = 30.0
    gf.sigma = 1.0
    gf.max_points = 500
    gf.block_radius = gf.scale*1.0
    gf.search_radius = gf.scale*2.0
    gf.nrealizations = 1
    nx = 60
    ny = 60
    nd = 1000

    xy = 100 * np.random.rand(nd, 2)
    #val = np.random.rand(nd, 1)
    val = np.ones((nd, 1))
    gf.data = np.hstack((xy, val))
    x = np.linspace(0, 100, nx)
    y = np.linspace(0, 100, ny)
    xv, yv = np.meshgrid(x,y)
    sample_locations = np.vstack((xv.flatten(), yv.flatten())).transpose()
    gf.sample_locations = sample_locations
    #gf.seq_simulation()
    #gf.conditional_simulation()
    gf.krige_seq()
    pass





    pass



