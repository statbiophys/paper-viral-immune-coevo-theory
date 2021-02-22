
"""
TAKES THE OUTPUT OF CLUSTERING SCRIPT AND INFERS THE LINEAGES PERSISTENCE TIME

  *  Copyright (C) 2021 Jacopo Marchi
  * 
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *  You should have received a copy of the GNU General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  * 
"""

import matplotlib
matplotlib.use('Agg')
import os
import glob

import sys
sys.path.append('..')
from lib.mppaper import *
import lib.mpsetup as mpsetup
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from scipy.optimize import fsolve
from scipy.special import factorial
from scipy.spatial.distance import pdist
import itertools
from scipy.spatial.distance import squareform
import shutil
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from sklearn.preprocessing import StandardScaler
import gc
import math
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import IncrementalPCA
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from scipy import interpolate
from scipy.interpolate import splrep, splev
from scipy import stats
from sklearn.neighbors import KDTree
from scipy.optimize import curve_fit
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from statsmodels.tsa.stattools import acf


def get_cmap(N):
    ''' Returns a function that maps each index in 0, 1, ...
        N-1 to a distinct RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    #scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap='jet')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color




def find_nearest_idx(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx


def avg_diff_lengths_1d(arrays_list): # given a list of (n_reals) 1D arrays, compute the average accounting for length differences
    #if len(arrays_list)>0:
    maxlen=max([a.size for a in arrays_list])
    arr = np.ma.empty((maxlen,len(arrays_list)))
    arr.mask = True
    
    for i, a in enumerate(arrays_list):
        arr[:a.shape[0],i] = a
    
    print arr.mean(axis = 1)
    
    return arr.mean(axis = 1)
#    else:
#	return np.array([0.])


def avg_diff_lengths_1d_count(arrays_list): # given a list of (n_reals) 1D arrays, compute the average accounting for length differences
    #if len(arrays_list)>0:
    maxlen=max([a.size for a in arrays_list])
    arr = np.ma.empty((maxlen,len(arrays_list)))
    arr.mask = True
    
    for i, a in enumerate(arrays_list):
        arr[:a.shape[0],i] = a
    
    print arr.count(axis = 1)
    
    return arr.count(axis = 1)
#    else:
#	return np.array([0.])

    

def avg_diff_lengths_1d_stderr(arrays_list): # given a list of (n_reals) 1D arrays, compute the average accounting for length differences
    #if len(arrays_list)>0:
    maxlen=max([a.size for a in arrays_list])
    arr = np.ma.empty((maxlen,len(arrays_list)))
    arr.mask = True
    
    for i, a in enumerate(arrays_list):
        arr[:a.shape[0],i] = a
    stderrs=arr.std(axis = 1) / np.sqrt(arr.count(axis = 1))
    print stderrs
    
    return stderrs
#    else:
#	return np.array([0.])

    

    
def max_delta_i_fct(times, delta_t_thr):
    
    for i in range(times.size-2,-1,-1): 
	
	delta_t=np.mean(times[i:] - times[:times.size - i])
	if delta_t <= delta_t_thr:
	    return i
    return times.size



def align_times(arrays_list): # return offset indexes aligned according to times
    first_times= np.array([a[0] for a in arrays_list])
    index_min_first_time=np.argmin(first_times)
    ref_times=arrays_list[index_min_first_time]
    alignment_indexes =  np.array([find_nearest_idx(ref_times, a[0]) for a in arrays_list]) # indexes that I need to align the first times, wrt the array with minimum first time
    return alignment_indexes
    


def avg_diff_lengths_1d_aligned(arrays_list, alignment_indexes): # given a list of (n_reals) 1D arrays, compute the average accounting for length differences, after aligning according to times
    
    
    indmax=np.amax(alignment_indexes)
    #maxlen=max([a.size for a in arrays_list])
    maxlen=max([a.size + alignment_indexes[i]  for i, a in enumerate(arrays_list)])
    arr = np.ma.empty((maxlen,len(arrays_list)))
    arr.mask = True
    
    for i, a in enumerate(arrays_list):
        arr[alignment_indexes[i]:a.shape[0] + alignment_indexes[i],i] = a

    
    return arr.mean(axis = 1)

col_dict=colors.cnames

def label_line(line, label, x, y, color='0.5'):
    """Add a label to a line, at the proper angle.

    Arguments
    ---------
    line : matplotlib.lines.Line2D object,
    label : str
    x : float
        x-position to place center of text (in data coordinated
    y : float
        y-position to place center of text (in data coordinates)
    color : str
    size : float
    """
    xdata, ydata = line.get_data()
    x1 = xdata[0]
    x2 = xdata[-1]
    y1 = ydata[0]
    y2 = ydata[-1]

    ax = line.get_axes()
#    text = ax.annotate(label, xy=(x, y), xytext=(-10, 0),
#                       textcoords='offset points',
#                       color=color,
#                       horizontalalignment='left',
#                       verticalalignment='bottom')
 
    text = ax.annotate(label, xy=(x, y), xytext=(1, -6),
                       textcoords='offset points',
                       color=color,
                       horizontalalignment='right',
                       verticalalignment='bottom')

    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))

    rise = (sp2[1] - sp1[1])
    run = (sp2[0] - sp1[0])

    slope_degrees = np.degrees(np.arctan2(rise, run))
    text.set_rotation(slope_degrees)
    return text

#		    #ax.annotate(r'$\propto r^2$',
#			#xy=(powlaw_x[powlaw_x.size/2], powlaw_y[powlaw_y.size/2]*4./3), xycoords='data', rotation= np.angle((powlaw_x[-1]/powlaw_x[0])*(xmin1/xmax1) + 1j * (powlaw_y[-1]/powlaw_y[0])*(ymin1/ymax1), deg=True),
#		    ax.annotate(r'$\propto r^2$',
#			xy=(powlaw_x[powlaw_x.size/2], 2*powlaw_y[powlaw_y.size/2]*4./3), xycoords='data', rotation= np.angle((powlaw_x[-1]/powlaw_x[0])**(1./(xmax1- xmin1)) + 1j * (powlaw_y[-1]/powlaw_y[0])**(1./(ymax1- ymin1)), deg=True), rotation_mode='anchor',
#			horizontalalignment='right', verticalalignment='bottom', color='grey')
#		    #

def mystderr(v):
    
    return np.std(v)/np.sqrt(v.shape[0])

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between_ND_pos(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def angle_between(v1, v2):
    
    #return np.arctan2(v1[1], v1[0]) - np.arctan2(v2[1], v2[0])
    
    
    # rotate v1 in reff of v2
    v2_perpdir=np.array([v2[1], -v2[0]])
    
    Q = np.squeeze(np.array([v2, v2_perpdir])).T
    #print 'Transformed to original system with \n Q={}'.format(Q)
    #print 'Orthogonality check \n {}'.format(Q.dot(Q.T))
    #print 'Orthogonality check \n {}'.format(np.dot(Q.T, X[0,:2]) == np.linalg.solve(Q, X[0,:2]))
    #
    #print np.dot(Q.T, points_in_clust[0,:])
    #print np.linalg.solve(Q, points_in_clust[0,:])
    
    v1_newbasis=np.dot(Q.T, np.asarray(v1).T).T
    return np.arctan2(v1_newbasis[1], v1_newbasis[0]) 



def acf_decay(x, autocorr):
    return np.exp(-x/autocorr)

def rw_msd_fct_time(x, diff, c):
    return x*diff + c # const for sampling noise high frequency

def ballistic_msd_constvel(x, p_l):
    return 2*x*p_l*(1 - (p_l/x)*( 1 - np.exp(-x/p_l)) )

    

def ballistic_msd_constvel_renorm(x, p_l, ltot): # not sure about this
    return 2*x*p_l*(ltot**2)*(1 - (p_l/x)*( 1 - np.exp(-x/p_l)) )


def ballistic_msd_constvel_fct_time(x, p_l, v): # in fact this is equal to ballistic_msd_constvel_renorm, I keep them separated for readability
    return 2*x*p_l*(v**2)*(1 - (p_l/x)*( 1 - np.exp(-x/p_l)) )


def ballistic_msd_varvel_fct_time_nobeta(x, p_l, v, corr_l, var_v): # 
    return 2*x*p_l*(v**2)*(1 - (p_l/x)*( 1 - np.exp(-x/p_l)) )  +  2*x*corr_l*(var_v)*(1 - (corr_l/x)*( 1 - np.exp(-x/corr_l)) )


def ballistic_msd_varvel_fct_time(x, p_l, v, beta, var_v): # beta is autocorreltion rate for velocity
    
    corr_l=1./(beta + 1./p_l)
    return 2*x*p_l*(v**2)*(1 - (p_l/x)*( 1 - np.exp(-x/p_l)) )  +  2*x*corr_l*(var_v)*(1 - (corr_l/x)*( 1 - np.exp(-x/corr_l)) )


    

def ballistic_msd_varvel(x, p_l, beta, var_v_o_vsq): # beta is autocorreltion rate for velocity over average speed
    
    corr_l=1./(beta + 1./p_l)
    return 2*x*p_l*(1 - (p_l/x)*( 1 - np.exp(-x/p_l)) )  +  2*x*corr_l*(var_v_o_vsq)*(1 - (corr_l/x)*( 1 - np.exp(-x/corr_l)) )


    

def func_saturate(x, a,  c):
    return a*(1 -  np.exp(-c * x) )



def func_powerlaw(x, m, c):
    return  x*m  +  c


def func_powerlaw_neut(x, c):
    return  c*x**-1. 

def func_powerlaw_neut_lin(x, c):
    return  -x + c 

def truncated_powerlaw(x): #ks test wants the cdf, not the pdf!!
    #SFS_alive_viruses_alltime[SFS_alive_viruses_alltime>0].size/np.sum(1./x)
    
    xmin=np.amin(SFS_alive_viruses_alltime[SFS_alive_viruses_alltime>0])
    
    if np.isscalar(x):
        if x<=xmin:
            return 0.
        if x>=xmin and x<=1:
            return 1 - (np.log(x)/np.log(xmin))
        else:
            return 1.
    else:
        
        print "in kstest"
        
        print xmin
        res = np.zeros(x.shape)
        
        #res=np.where((x>=np.amin(SFS_alive_viruses_alltime[SFS_alive_viruses_alltime>0])) & (x<=1), -1./(np.log(np.amin(SFS_alive_viruses_alltime[SFS_alive_viruses_alltime>0]))*x), 0.)
        
        
        #res=np.where((x>=np.amin(SFS_alive_viruses_alltime[SFS_alive_viruses_alltime>0])) & (x<=1), 1./x, 0.)
        res=np.where((x>=xmin) & (x<=1), 1 - (np.log(x)/np.log(xmin)), 0.)
        res=np.where((x<=1), res, 1.)
        #res=np.where((x>=np.amin(SFS_alive_viruses_alltime[SFS_alive_viruses_alltime>0])) & (x<=1), 1./x, 0.)
        
        #return res/simps(res, x)
        return res
	
	
def secder_double_smooth(end2end_coarse_gr_list_reals_mean, path_l_serie_list_reals_mean, window):
    
	
    end2end_coarse_gr_list_smooth =np.convolve(end2end_coarse_gr_list_reals_mean, np.ones((window,))/window, mode='valid')
    path_l_serie_list_smooth      =np.convolve(path_l_serie_list_reals_mean, np.ones((window,))/window, mode='valid')
    
    
    h_deriv=path_l_serie_list_smooth[1:] - path_l_serie_list_smooth[:-1]
    
    firstder_smooth=(end2end_coarse_gr_list_smooth[1:] -  end2end_coarse_gr_list_smooth[:-1])/h_deriv
    
    
    
    path_l_serie_list_smooth_smooth      =np.convolve(path_l_serie_list_smooth[:-1], np.ones((window,))/window, mode='valid')
    firstder_smooth      =np.convolve(firstder_smooth, np.ones((window,))/window, mode='valid')
    
    h_deriv=path_l_serie_list_smooth_smooth[1:] - path_l_serie_list_smooth_smooth[:-1]
    
    #secder_smooth=(firstder_smooth[1:] -  firstder_smooth[:-1])/(h_deriv[:-1])
    secder_smooth =(firstder_smooth[1:] -  firstder_smooth[:-1])/(h_deriv)
    
    #end2end_der_path_l = (-1*end2end_coarse_gr_list_reals_mean[:-4]+16*end2end_coarse_gr_list_reals_mean[1:-3]-30*end2end_coarse_gr_list_reals_mean[2:-2]+16*end2end_coarse_gr_list_reals_mean[3:-1]-1*end2end_coarse_gr_list_reals_mean[4:])/(12*1.0*h_deriv**2)
    
    #path_l_serie_list_reals_mean=path_l_serie_list_reals_mean[2:-2]
    
    return [secder_smooth, path_l_serie_list_smooth_smooth[:-1]]
    

def vel_fromparams(recog_width, diff_const, mem_points, F_0, vir_number):
    return (diff_const**(2/3.))*((mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)**(1/3.))*((24*np.log((diff_const**(1/3.)) * ((mem_points*(np.exp(F_0/mem_points) -1.)/recog_width)**(2/3.)) * vir_number))**(1/3.)) 
    

def sig_fromparams(recog_width, diff_const, mem_points, F_0, vir_number):
    return (diff_const**(1/3.))*((mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)**(-1/3.))*((24*np.log((diff_const**(1/3.)) * ((mem_points*(np.exp(F_0/mem_points) -1.)/recog_width)**(2/3.)) * vir_number))**(1/6.)) 
    

def sel_fromparams(recog_width, mem_points, F_0):
    return (mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)
    

def vtauor_fromparams(recog_width, mem_points, F_0):
    return (1./(np.exp(F_0/mem_points)-1.))
    

def vtau_fromparams(recog_width, mem_points, F_0):
    return (recog_width/(np.exp(F_0/mem_points)-1.))
    

def tau_fromparams(mem_points, host_number, vir_number):
    return (mem_points*host_number)/vir_number.astype(float)
    
        
    

def R_persl_fromparams(recog_width, mem_points, F_0):
     
    return mem_points*recog_width/(sel_fromparams(recog_width, mem_points, F_0)*recog_width + mem_points) 
    
    

#dir_in='data/1d_{progr_tol}_{maxiters}_{algo}/opt.dat'.format(progr_tol=intra_host_neut.progr_tol, maxiters=int(intra_host_neut.maxiters), algo=intra_host_neut.algo)
#dir_in=intra_host_neut.file_o
#dir_io='../../contagion_simulation_results/data_neutral/DIM_2_people_number_10000_maxinfections_500001_alpha_4d000000_mu_1d000000_recog_width_5d000000_sigma_0d500000_f_m_0d100000_Rnot_mode_rec_force' # directory with input files
dir_io=sys.argv[1] # directory with input files
dir_in_tot='{inp}/realizations'.format(inp=dir_io) # directory with input files
#dir_in = dir_in.replace(".", "d")

dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots
#dir_out_plots = dir_out_plots.replace(".", "d")
   
def get_cmap(N):
    ''' Returns a function that maps each index in 0, 1, ...
        N-1 to a distinct RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    #scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap='jet')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

    
    
param_file='{inp}/parameters_backup.dat'.format(inp=dir_io)

params = np.genfromtxt(param_file, dtype="f8,f8,f8,i8,i8,i8,i8, |S10, |S10, |S10, i8, i8, f8", names=['mu','rec_width','jump_size', 'ppl_num', 'maxinfections','save_full_time','n_real','initial_condition','fake_initial_condition','phylo_subsample','mem_points','t_init','F0'])
in_cond=params['initial_condition']
t_init=params['t_init']

print params
print params.shape
n_real=params['n_real']
num_ppl=params['ppl_num']
recog_width=params['rec_width']
mu=params['mu']
mem_points=params['mem_points']
F0=params['F0']
latt_sp=1
traj_anal=True

#if n_real >10:
#    n_real=10

print n_real   

thisfigsize = figsize
thisfigsize[1] *= 0.75

# get data

for real in np.arange(1,n_real+1):
    dir_in='{inp}/realization_{real}'.format(inp=dir_in_tot,real=real) # directory with input files
    dir_out_plots='{inp}/realization_{real}/pers_l'.format(inp=dir_out_plots_tot,real=real) # directory with output plots
    dir_out_frames_zoom='{inp}/frames_zoom'.format(inp=dir_out_plots) # directory with output plots
    #dir_in_frames='{inp}/frames_npz_compr'.format(inp=dir_in,real=real) # directory with input files
    dir_in_frames='{inp}/frames/'.format(inp=dir_in) # directory with input files
    

    file_exploded='{inp}/expl_file.txt'.format(inp=dir_in)
    file_extinct='{inp}/extinct_file.txt'.format(inp=dir_in)
    
     
    timelast=0
    file_in_avg_npz_compr='{inp}/evo_mean_stats_real_{real}.dat'.format(inp=dir_in, real=real)
    
    if os.path.isfile(file_in_avg_npz_compr):
        
        
        data=[]
        with open(file_in_avg_npz_compr) as f:
            lines=f.readlines()
            for line in lines:
                if not line.startswith("#"):
                    line=line.decode('utf-8','ignore').encode("utf-8")
                    myarray = np.fromstring(line, dtype=float, sep=' ')
                    data.append(myarray)
                    #print(myarray)
            print len(data)
        
        data = np.asarray(data)
        print data.shape
    
                
        print data
        print data.ndim
        
        if data.ndim > 1:    
        
    
            
            timelast=data[-1,0]
            
    if not os.path.exists(dir_out_plots):
        os.makedirs(dir_out_plots)	
       
       
      
    
    #CLUSTER TRACKING! USINH OPT KVAR
    print "CLUSTER TRACKING! USING OPT KVAR"
    
    avg_clust_size = 0.
    std_clust_size = 0.
    SIG_SPACE = 0.
    avg_var_displ_x=0
 

    file_in='{inp}/global_features_{real}.txt'.format(inp=dir_in, real=real)
     
    if os.path.isfile(file_in):
     
        data_glob = np.loadtxt(file_in)
        
        print data_glob
        print data_glob.shape
        
        if data_glob.ndim >1:
            data_glob=data_glob[-1]

        print data_glob
        print data_glob.shape

        #~ data = np.array([avg_I, std_I, var_log_I, avg_fit_avg, var_fit_avg, avg_fit_var, var_fit_var, avg_fit_nose, var_fit_nose, avg_vel, var_vel, avg_num_IS_coords, avg_num_vir_coords, avg_mean_displ_x, avg_mean_displ_y, avg_var_displ_x, avg_var_displ_y, avg_mean_displ_tot, diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, avg_sig_est, var_sig_est])
    
        avg_I                            = data_glob[0]
        avg_var_displ_x                  = data_glob[15]
        
        D_mod=mu*(avg_var_displ_x/2.) # simulations units!
        SIG_SPACE    = sig_fromparams(recog_width, D_mod, mem_points, F0, avg_I*num_ppl)
     
     
        avg_clust_size = SIG_SPACE


    
    D_mod=mu*(avg_var_displ_x/2.) # simulations units!
    
    # ~ rec_width                                     = jump_size/rec_width # /(jump_size), units space already in false rec width, so here keep simulations units
    R_PERS_L     = R_persl_fromparams(recog_width, mem_points, F0)# /(jump_size) # units space
    PERS_TIME=R_PERS_L**2./(4.*D_mod)
    VEL_Nteo          = vel_fromparams(recog_width, D_mod, mem_points, F0, avg_I*num_ppl)
    
    timescale_system = R_PERS_L/VEL_Nteo

    print mu, recog_width, mem_points, F0, avg_var_displ_x, avg_I
    print VEL_Nteo
    print R_PERS_L, timescale_system, PERS_TIME
      
    file_in='{inp}/space_clustering_fct_time_real_{real}.dat.npz'.format(inp=dir_in, real=real)
    
    if os.path.isfile(file_in):
     
        z = np.load(file_in)    
        data_vir_time = z['data']
        data_vir_time_allclust = z['data_allclust']
        data_vir_time_grad = z['data_allclust_grads']
        
         #~ np.savez_compressed(file_out,data = data, data_allclust = data_allclust, data_allclust_grads = data_allclust_grads, frac_1cl=frac_1cl, times_sharp_turn_split_alltimes_allclust= times_sharp_turn_split_alltimes_allclust)
       
       
         
        vel_clust_split_alltimes_allclust        =  data_vir_time_allclust[0,:]
        var_parall_clust_split_alltimes_allclust =  data_vir_time_allclust[1,:]
        var_perp_clust_split_alltimes_allclust   =  data_vir_time_allclust[2,:]
        var_tot_clust_split_alltimes_allclust    =  data_vir_time_allclust[3,:]
        
        
        print data_vir_time_allclust.shape
        print data_vir_time_allclust[0,:].mean()
        
        
          
        # ~ avg_vel_clust_split_alltimes_allclust         = vel_clust_split_alltimes_allclust         .mean()           
        # ~ avg_var_parall_clust_split_alltimes_allclust  = var_parall_clust_split_alltimes_allclust  .mean()           
        # ~ avg_var_perp_clust_split_alltimes_allclust    = var_perp_clust_split_alltimes_allclust    .mean()           
        # ~ avg_var_tot_clust_split_alltimes_allclust     = var_tot_clust_split_alltimes_allclust     .mean()           
        # ~ var_vel_clust_split_alltimes_allclust         = vel_clust_split_alltimes_allclust         .var()   
        # ~ var_var_parall_clust_split_alltimes_allclust  = var_parall_clust_split_alltimes_allclust  .var()   
        # ~ var_var_perp_clust_split_alltimes_allclust    = var_perp_clust_split_alltimes_allclust    .var()   
        # ~ var_var_tot_clust_split_alltimes_allclust     = var_tot_clust_split_alltimes_allclust     .var()  
        
        avg_clust_size = (3.*np.sqrt(var_tot_clust_split_alltimes_allclust)).mean() # used to sparsify trajectories for fast population size fluctuations
        std_clust_size = (3.*np.sqrt(var_tot_clust_split_alltimes_allclust)).std()


  
        print var_tot_clust_split_alltimes_allclust
        print var_tot_clust_split_alltimes_allclust.shape
        
    print avg_clust_size
    print std_clust_size
    print SIG_SPACE
    
    
    # TIME THRESHOLDS TO SPLIT TRAJECTORIES FOR MSD, I WANT ENOUGH SUBTRAJECTORIES STATISTICS, BUT ALSO LONG ENOUGH TRACES WITHRESPECT TO T AND SMOOTHING WINDOW
    smax=min(max(timelast/13, 30000), 200000)
    splitsize_typ =  min(max(20*timescale_system, 30000), smax)
    # ~ splitsize_typ =  100000
    window_typ = 5
    time_thr   = min(2*timescale_system, splitsize_typ/2)
    # ~ time_thr   = min(1000, splitsize_typ/3)
    
    print timelast, timescale_system
    print time_thr, splitsize_typ
    
    # ~ sys.exit()
  
    
    hist_tdiff_trajclade                =[]
    hist_end2end_trajclade              =[]
    hist_length_trajclade               =[]
    hist_end2end_o_length_trajclade     =[]
    
    

    xs_turns =[]
    ys_turns =[]
    
    
    hist_nofull_time_turns              =[]
    hist_nofull_time_turns_o_ttot_traj  =[]  
    
    hist_time_turns              =[]
    hist_time_turns_o_ttot_traj  =[]  
    hist_num_turns               =[]   
    hist_num_turns_o_tclade               =[]   
           

    
    angles_subseg_fct_time_list                     =[]
    end2end_coarse_gr_fct_time_list                     =[]
    time_diff_traj_list                                  =[]
    time_diff_o_ttot_traj_list             =[]
                                           
    end2end_coarse_gr_fct_time_list_binstat=[]
    time_diff_traj_list_binstat            =[]
    time_diff_o_ttot_traj_list_binstat     =[]
    
    end2end_coarse_gr_list              =[]
    path_l_serie_list                   =[]
    path_l_o_ltot_serie_list            =[]
    end2end_coarse_gr_o_path_l_list     =[]
    
    
    end2end_coarse_gr_list_binstat              =[]
    path_l_serie_list_binstat                   =[]
    path_l_o_ltot_serie_list_binstat            =[]
   
    angles_subseg_fct_time_list_binstat=[]
    time_diff_angles_traj_list_binstat=[]
    
    time_list_binstat=[]
    angles_vs_time_list_binstat=[]
    tlags_autocorr_list_binstat=[]
    angles_autocorr_list_binstat=[]
    vel_autocorr_list_binstat=[]
    
    
    time_raw_list_binstat            =[]
    angles_vs_time_raw_list_binstat  =[]
    tlags_autocorr_raw_list_binstat =[]
    angles_autocorr_raw_list_binstat=[]
    
    data_track_tot     =[]
    
    last_centr_ID=0
    
    xs_all_centr_all_time =[] # list of the lists of all xs a centroid had for all times, in ID order
    ys_all_centr_all_time =[] # list of the lists of all ys a centroid had for all times, in ID order


    # DATA, OUTPUT OF CLUSTERING SCRIPT
    outfile_track='{inp}/clusters_track_real_{real}.dat'.format(inp=dir_in, real=real)

    
        
    if os.path.isfile(outfile_track) and traj_anal:
        
        data_track_tot	    = np.loadtxt(outfile_track)
        
        print data_track_tot.shape
        
        # data_son = np.array([time, state, ID_son, t_birth_son,  centr_son[0], centr_son[1], num_vir_son, size_son, ID_fath_son, centr_fath_son[0], centr_fath_son[1], num_vir_fath_son, size_fath_son])
                
        if data_track_tot.shape[0]>1:
    
            fathers_split= np.unique(data_track_tot[ (data_track_tot[:,1]>1),8])
            
            print fathers_split
            print fathers_split.shape
                    
            for i_centr, id_new in enumerate(np.unique(data_track_tot[:,2])):
                
                # e2e fct time
                    
                track_id_orig=data_track_tot[data_track_tot[:,2]==id_new,:]
                # ~ print id_new, track_id_orig.shape
                
                split_end=(id_new in fathers_split)
                
                
                x=track_id_orig[:,4]
                y=track_id_orig[:,5]
                t_traj = track_id_orig[:,0] #- track_id_orig[:,3]
                
                
                

                
                x_angles= x[::1]
                y_angles= y[::1]
            
            
                clust_dir_x=(x_angles[1:] -  x_angles[:-1])/np.sqrt((x_angles[1:] - x_angles[:-1])**2 + (y_angles[1:] - y_angles[:-1])**2)
                
                clust_dir_y=(y_angles[1:] - y_angles[:-1])/np.sqrt((x_angles[1:] - x_angles[:-1])**2 + (y_angles[1:] - y_angles[:-1])**2)
                
                if clust_dir_x.size>10:
                    
                    #clust_dir_angles=np.angle(clust_dir_x + 1j * clust_dir_y, deg=False) 
                    clust_dir_angles= np.unwrap(np.asarray([angle_between([dir_x , clust_dir_y[i_dir_x]],[clust_dir_x[0] , clust_dir_y[0]]) for i_dir_x, dir_x in enumerate(clust_dir_x)]))
                    clust_dir_angles= clust_dir_angles - clust_dir_angles[0]
                    
                    
                    t_traj_angles= t_traj[1::1] - t_traj[0]
                    
                    #~ print clust_dir_angles
                    # ~ print clust_dir_angles.shape
                    # ~ print t_traj_angles.shape
                    
                    
                
                    
                    autocorr_span=min(1000, int(clust_dir_angles.size*4./5))
                    
                    # ~ print autocorr_span
                    
                    autocorr_ang=acf(clust_dir_angles, fft=True, nlags=autocorr_span)
                    autocorr_tlag = t_traj_angles[np.arange(autocorr_span + 1)] - t_traj_angles[0]
                    
                
    
                    # ~ time_raw_list_binstat .extend(t_traj.tolist())        
                    # ~ angles_vs_time_raw_list_binstat .extend(clust_dir_angles.tolist())        
                    tlags_autocorr_raw_list_binstat .extend(autocorr_tlag.tolist())        
                    angles_autocorr_raw_list_binstat.extend(autocorr_ang.tolist()) 
                    
                    print "print angles timeseries"
    
                    print        clust_dir_angles
    
                
                
                
            print "CALCULATING PERSISTENCE LENGTH"
            
            print data_track_tot[:,3] - data_track_tot[:,0]
            
            print data_track_tot[ (data_track_tot[:,1]>1),8]
            
            # ids that go extinct after less than 1000 cycles:
            ids_ext_long=data_track_tot[((data_track_tot[:,0] - data_track_tot[:,3]) < 1000) & (data_track_tot[:,1]==0),2] 
            
            print ids_ext_long
            print ids_ext_long.shape
            
            ids=data_track_tot[:,2]
            
            mask = np.isin(ids, ids_ext_long)
            
            
            
            data_track_ext_long= data_track_tot[~mask,:]
            
            data_split=data_track_ext_long[ (data_track_ext_long[:,1]>1),:]
            
            fathers_split= data_track_ext_long[ (data_track_ext_long[:,1]>1),8]
            
            fathers_split_unique, fathers_split_unique_idx, fathers_split_unique_counts = np.unique(fathers_split, return_index=True, return_counts=True) 
            
            fathers_split_unique=fathers_split_unique[np.argsort(fathers_split_unique_idx)]
            fathers_split_unique_counts=fathers_split_unique_counts[np.argsort(fathers_split_unique_idx)]
            fathers_split_unique_idx=fathers_split_unique_idx[np.argsort(fathers_split_unique_idx)]
            
            print fathers_split
            print fathers_split.shape
            
            print fathers_split_unique
            print fathers_split_unique_counts
            print fathers_split_unique_idx
            print fathers_split_unique.shape
            
            
            idxs_1branch=fathers_split_unique_idx[ fathers_split_unique_counts==1 ]
            
            print idxs_1branch
            print idxs_1branch.shape
            
            data_track_ext_long_orig=data_track_ext_long.copy()
            
            # REMOVE SHORT LINEAGES AND JOIN CONSECUTIVE SINGLE LINEAGES
            
            #data_son = np.array([time, state, ID_son, t_birth_son,  centr_son[0], centr_son[1], num_vir_son, size_son, ID_fath_son, centr_fath_son[0], centr_fath_son[1], num_vir_fath_son, size_fath_son])
            for i, idx in enumerate(fathers_split_unique_idx):
                print i, idx, data_split[idx,0], data_split[idx,8], data_split[idx,2], (idx in idxs_1branch)
        
                #father_data=data_track_ext_long[np.argmax(data_track_ext_long_orig[:,2]==data_split[idx,8]), :]
                new_id=data_track_ext_long[np.argmax(data_track_ext_long_orig[:,2]==data_split[idx,8]), 2]
                father_data=data_track_ext_long[data_track_ext_long[:,2]==new_id, :]
                
                print new_id, data_split[idx,8]
                
                #~ print father_data[:,0]
                #~ print father_data[:,2]
                #~ print father_data[:,3]
                print np.amin(father_data[:,0] - father_data[:,3])
                father_data=father_data[np.argmax(father_data[:,0] - father_data[:,3])]
                
                print father_data[0] - father_data[3]
        
                
                if idx in idxs_1branch: #substitute id of the son with the one of the father
                #data_track_ext_long[data_track_ext_long[:,8]==fathers_split_unique[i], 2]=fathers_split_unique[i]
                
                    print "joining single branch"
                    
                    
                    
                    data_track_ext_long[(data_track_ext_long_orig[:,8]==data_split[idx,8]) & (data_track_ext_long_orig[:,1]>1), 1]=1 # as if it continued tracking 
                    #data_track_ext_long[data_track_ext_long[:,8]==data_split[idx,8], 2]=data_split[idx,8]
                    
                    prev_id=data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 2]
                    #print prev_id
                    data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 2] =father_data[2]
                    data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 3] =father_data[3]
                    data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 9] =father_data[9]
                    data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 10]=father_data[10]
                    data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 11]=father_data[11]
                    data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 12]=father_data[12]
                    
                    
                    data_track_ext_long[data_track_ext_long_orig[:,8]==data_split[idx,8], 8] =father_data[8] # after this cannot be recognized anymore
                    #print prev_id
                    
                    
                
                if (father_data[0] - father_data[3] <1000) and not (idx in idxs_1branch) and not (idx == 1):
                    
                    print "joining short split"
                    father_father_data=data_track_ext_long[np.argmax(data_track_ext_long[:,2]==father_data[8]), :]
            
                    data_track_ext_long[(data_track_ext_long[:,2]==new_id) & (data_track_ext_long[:,1]>1), 1]=1 # as if it continued tracking 
                    data_track_ext_long[data_track_ext_long[:,2]==new_id, 3] =father_father_data[3]
                    data_track_ext_long[data_track_ext_long[:,2]==new_id, 9] =father_father_data[9]
                    data_track_ext_long[data_track_ext_long[:,2]==new_id, 10]=father_father_data[10]
                    data_track_ext_long[data_track_ext_long[:,2]==new_id, 11]=father_father_data[11]
                    data_track_ext_long[data_track_ext_long[:,2]==new_id, 12]=father_father_data[12]
                    data_track_ext_long[data_track_ext_long[:,2]==new_id, 8] =father_father_data[8] # after this cannot be recognized anymore
                    data_track_ext_long[data_track_ext_long[:,2]==new_id, 2] =father_father_data[2]
            
            
            
                
            print data_track_ext_long[:,2]
            print data_track_ext_long_orig[:,2]
            print np.unique(data_track_ext_long_orig[:,2])
            ids_new=np.unique(data_track_ext_long[:,2])
            
            print ids_new
            print ids_new.shape
            print data_track_ext_long.shape
            
        
            print "calculating end2end" # I could do this on raw data in fact. Quick fix in case
        
            #data_son = np.array([time, state, ID_son, t_birth_son,  centr_son[0], centr_son[1], num_vir_son, size_son, ID_fath_son, centr_fath_son[0], centr_fath_son[1], num_vir_fath_son, size_fath_son])
            fathers_split= np.unique(data_track_ext_long[ (data_track_ext_long[:,1]>1),8])
            
            print fathers_split
            print fathers_split.shape
            
                
            if data_track_tot.shape[0]>1:
                    
                    
                    
                
                for i_centr, id_new in enumerate(np.unique(data_track_ext_long[:,2])): # cycle on lineages
                
                        
                    # e2e fct time
                        
                    track_id_orig=data_track_ext_long[data_track_ext_long[:,2]==id_new,:]
                    print id_new, track_id_orig.shape
                    
                    split_end=(id_new in fathers_split)
                    
                    
                    x      = track_id_orig[:,4]
                    y      = track_id_orig[:,5]
                    t_traj = track_id_orig[:,0] #- track_id_orig[:,3]
                    
                    tdiff_clade=t_traj[-1] - t_traj[0]
                
                
                    
                    splitsize = min(splitsize_typ, tdiff_clade) # 
                    
                    
                    nseg=int(tdiff_clade/(splitsize+1)) +1
                    
                    print splitsize, tdiff_clade, nseg
                    
                    #~ print t_traj
                    
                    if x.size>3:
                    
                        dist_thresh_angles=2.*(avg_clust_size + std_clust_size)
                        idxs_angles_sparse= []
                        
                        i_x=0
                        
                        # SPARSIFY TRAJECTORY FOR FAST POPULATION SIZE FLUCTUATIONS
                        
                        while i_x < x.size :
                            idxs_angles_sparse.append(i_x)
                            dists=np.sqrt((x[i_x:] - x[i_x])**2 + (y[i_x:] - y[i_x])**2)
                            dists_previous=np.zeros_like(dists) + 2.*dist_thresh_angles
                            if len(idxs_angles_sparse)>1:
                                dists_previous=np.sqrt((x[i_x:] - x[idxs_angles_sparse[len(idxs_angles_sparse)-2]])**2 + (y[i_x:] - y[idxs_angles_sparse[len(idxs_angles_sparse)-2]])**2)
                            
                            
                            i_incr=np.argmax((dists>0.) & (dists>=dist_thresh_angles) & (dists_previous>=2.*dist_thresh_angles)) # stops at the first True...
                            i_x= i_x + i_incr
                            
                            print i_x,  len(idxs_angles_sparse), dist_thresh_angles
                            print dists
                            print dists_previous
                            
                            if i_incr==0: # ...but returns 0 if no True
                                break
                    
                        idxs_angles_sparse=np.asarray(idxs_angles_sparse)
                        
                        x_angles_sparse= x[idxs_angles_sparse].copy()
                        y_angles_sparse= y[idxs_angles_sparse].copy()
                        times_angles_sparse= t_traj[idxs_angles_sparse].copy()
                        
                        print x_angles_sparse.size, x.size
                        print times_angles_sparse
                       
                        t_traj_angles= times_angles_sparse[1::1] - times_angles_sparse[0]
                        
                        
                        x_angles= x_angles_sparse.copy()
                        y_angles= y_angles_sparse.copy()
                    
                        if x_angles_sparse.size>3:
    
                            clust_dir_x=(x_angles[1:] -  x_angles[:-1])/np.sqrt((x_angles[1:] - x_angles[:-1])**2 + (y_angles[1:] - y_angles[:-1])**2)
                            
                            clust_dir_y=(y_angles[1:] - y_angles[:-1])/np.sqrt((x_angles[1:] - x_angles[:-1])**2 + (y_angles[1:] - y_angles[:-1])**2)
                            
                            
                            #clust_dir_angles=np.angle(clust_dir_x + 1j * clust_dir_y, deg=False) 
                            clust_dir_angles_orig= np.unwrap(np.asarray([angle_between([dir_x , clust_dir_y[i_dir_x]],[clust_dir_x[0] , clust_dir_y[0]]) for i_dir_x, dir_x in enumerate(clust_dir_x)]))
                            clust_dir_angles_orig= clust_dir_angles_orig - clust_dir_angles_orig[0]
                            
                            
                            
                            # smooth since it's a derivative
                            window=np.min([window_typ, clust_dir_angles_orig.size])
                            # ~ clust_dir_angles =np.convolve(clust_dir_angles_orig, np.ones((window,))/window, mode='same')
                            # ~ w=np.blackman(window)
                            w=np.ones((window,))/window
                            
                            clust_dir_angles =np.convolve(clust_dir_angles_orig, w/w.sum(), mode='same')
                            time_diff_angles=times_angles_sparse[:-1] - times_angles_sparse[0]
        
                            
                            #~ print clust_dir_angles
                            #~ print clust_dir_angles.shape
                            #~ print t_traj_angles.shape
                                
                            angles_vs_time_raw_list_binstat .extend(clust_dir_angles_orig.tolist())    
                            time_list_binstat           .extend(t_traj_angles.tolist())        
                            angles_vs_time_list_binstat .extend(clust_dir_angles.tolist())        
                            
        
                            print "print angles timeseries"
                            
                            print        clust_dir_angles
                            print        clust_dir_angles.shape
                            print        clust_dir_angles_orig.shape
                             
                            # ~ sys.exit()
                        
                            
                            print "fct time"
                            
                            if x_angles_sparse.size>3*window: # NEED TIME TRACE LONG W.R.T. SMOOTHING WINDOW
                                
                                
                                for i_seg in range(0,nseg): # cycle on subsegments
                                    time_start=t_traj[0] + i_seg*splitsize
                                    time_end= t_traj[0] + min((i_seg+1)*splitsize, tdiff_clade)
                                    
                                    
                                    new_x_coarse_gr=x[(t_traj>=time_start) & (t_traj<time_end)].copy()
                                    new_y_coarse_gr=y[(t_traj>=time_start) & (t_traj<time_end)].copy()
                                    new_t_traj_coarse_gr=t_traj[(t_traj>=time_start) & (t_traj<time_end)].copy()
                                    
                                    
                                    clust_dir_angles_subseg=clust_dir_angles[(times_angles_sparse[:-1]>=time_start) & (times_angles_sparse[:-1]<time_end)].copy()
                                    time_diff_angles_subseg=time_diff_angles[(times_angles_sparse[:-1]>=time_start) & (times_angles_sparse[:-1]<time_end)].copy()
                                    
                                    if clust_dir_angles_subseg.size > 1:
                                        
                                        
                                        end2end_coarse_gr_fct_time= (new_x_coarse_gr[1:] - new_x_coarse_gr[0])**2 + (new_y_coarse_gr[1:] - new_y_coarse_gr[0])**2 # all the R^2 in the reduced segments
                                
                                        #~ print end2end_coarse_gr_fct_time.shape
                                        
                                        
                                        time_diff=new_t_traj_coarse_gr[1:] - new_t_traj_coarse_gr[0]
                                        # ~ time_diff_angles=new_t_traj_coarse_gr_angles[1:] - new_t_traj_coarse_gr_angles[0]
                                        
                                        time_diff_angles_subseg=time_diff_angles_subseg[1:] - time_diff_angles_subseg[0]
                                        clust_dir_angles_subseg= (clust_dir_angles_subseg[1:] - clust_dir_angles_subseg[0])**2
                                        
                                        
                                        
                                        angles_subseg_fct_time_list.append(clust_dir_angles_subseg)        
                                        end2end_coarse_gr_fct_time_list.append(end2end_coarse_gr_fct_time)        
                                        time_diff_traj_list.append(time_diff)         
                                        time_diff_o_ttot_traj_list.append(time_diff/tdiff_clade)       
                                        
                                        end2end_coarse_gr_fct_time_list_binstat.extend(end2end_coarse_gr_fct_time.tolist())        
                                        time_diff_traj_list_binstat.extend(time_diff.tolist())         
                                        time_diff_o_ttot_traj_list_binstat.extend((time_diff/tdiff_clade).tolist())       
                                        
                                        
                                        time_diff_angles_traj_list_binstat.extend(time_diff_angles_subseg.tolist())         # NEW  
                                        angles_subseg_fct_time_list_binstat.extend(clust_dir_angles_subseg.tolist())        
                                        
                                        
                                        #~ end2end_coarse_gr_fct_time_list_reals_binstat.extend(end2end_coarse_gr_fct_time.tolist())        
                                        #~ time_diff_o_ttot_traj_list_reals_binstat.extend((time_diff/tdiff_clade).tolist())       
                                        
                                        #~ time_diff_traj_list_reals_binstat.extend(time_diff.tolist())         
                                    
                                hist_tdiff_trajclade.append(tdiff_clade)         
                                #~ hist_tdiff_trajclade_reals.append(tdiff_clade)         
                            
                                
                                
                            print "angles"
                            
                            
                            vel_clust = np.sqrt((x[:-1] - x[1:])**2 + (y[:-1] - y[1:])**2)/(t_traj[1:] - t_traj[:-1])
        
                        
                            print x.shape
                            print split_end
                            
                            
                            if clust_dir_angles.size>10:
                                    
                                autocorr_span=min(1000, int(clust_dir_angles.size*4./5))
                                
                                print autocorr_span
                                
                                autocorr_ang=acf(clust_dir_angles, fft=True, nlags=autocorr_span)
                                autocorr_tlag = t_traj_angles[np.arange(autocorr_span + 1)] - t_traj_angles[0]
                                
                                
                                autocorr_vel=acf(vel_clust, fft=True, nlags=autocorr_span)
                                
                            
                                tlags_autocorr_list_binstat .extend(autocorr_tlag.tolist())        
                                angles_autocorr_list_binstat.extend(autocorr_ang.tolist())        
                                vel_autocorr_list_binstat   .extend(autocorr_vel.tolist()) 
                                
                                
                                
                    
                    
                
                hist_tdiff_trajclade                      =np.asarray(hist_tdiff_trajclade                      )
                hist_end2end_trajclade                      =np.asarray(hist_end2end_trajclade                      )
                hist_length_trajclade                       =np.asarray(hist_length_trajclade                       )
                hist_end2end_o_length_trajclade             =np.asarray(hist_end2end_o_length_trajclade             )
                
                
                
                
                
                hist_nofull_time_turns              =np.asarray(hist_nofull_time_turns              )
                hist_nofull_time_turns_o_ttot_traj  =np.asarray(hist_nofull_time_turns_o_ttot_traj      )
                
                hist_time_turns              =np.asarray(hist_time_turns              )
                hist_time_turns_o_ttot_traj  =np.asarray(hist_time_turns_o_ttot_traj      )
                hist_num_turns               =np.asarray(hist_num_turns                   )
                hist_num_turns_o_tclade               =np.asarray(hist_num_turns_o_tclade                   )
                
                
                avg_end2end_coarse_gr_fct_time_list           = avg_diff_lengths_1d(end2end_coarse_gr_fct_time_list)  
                avg_end2end_coarse_gr_fct_time_list_count           = avg_diff_lengths_1d_count(end2end_coarse_gr_fct_time_list)  
                avg_time_diff_traj_list                       = avg_diff_lengths_1d(time_diff_traj_list            )  
                avg_time_diff_o_ttot_traj_list                       = avg_diff_lengths_1d(time_diff_o_ttot_traj_list            )  
                std_end2end_coarse_gr_fct_time_list           = avg_diff_lengths_1d_stderr(end2end_coarse_gr_fct_time_list) 
                
                #avg_end2end_coarse_gr_fct_time_list_count           = avg_end2end_coarse_gr_fct_time_list_count.compressed()
                avg_end2end_coarse_gr_fct_time_list           = avg_end2end_coarse_gr_fct_time_list.compressed()
                std_end2end_coarse_gr_fct_time_list           = std_end2end_coarse_gr_fct_time_list.compressed()
                avg_time_diff_traj_list                       = avg_time_diff_traj_list.compressed()
                avg_time_diff_o_ttot_traj_list                       = avg_time_diff_o_ttot_traj_list.compressed()
                
                
                # ~ avg_angles_subseg_fct_time_list           = avg_diff_lengths_1d(angles_subseg_fct_time_list)  
                # ~ std_angles_subseg_fct_time_list           = avg_diff_lengths_1d_stderr(angles_subseg_fct_time_list) 
                
                # ~ avg_angles_subseg_fct_time_list           = avg_angles_subseg_fct_time_list.compressed()
                # ~ std_angles_subseg_fct_time_list           = std_angles_subseg_fct_time_list.compressed()
                
                
                avg_end2end_coarse_gr_list           = np.array([0.])
                avg_path_l_serie_list                = np.array([0.])
                avg_path_l_o_ltot_serie_list         = np.array([0.])
                avg_end2end_coarse_gr_o_path_l_list  = np.array([0.])
                
                 
                # MSD STATISTICS FROM SUBTRAJECTORIES
                
                print "binned mean MSD fct time"
                
                
                #nbins=int((avg_time_diff_traj_list[-1] - avg_time_diff_traj_list[0])/(avg_time_diff_traj_list[1] - avg_time_diff_traj_list[0] ) )#*3
                
                bindiffs=np.asarray(time_diff_angles_traj_list_binstat)[1:] - np.asarray(time_diff_angles_traj_list_binstat)[:-1] 

                # ~ nbins=int((np.amax(time_diff_angles_traj_list_binstat)  - np.amin(time_diff_angles_traj_list_binstat))/np.mean(bindiffs[np.abs(bindiffs)<30000]) )#*3
                nbins=int((np.amax(time_diff_angles_traj_list_binstat)  - np.amin(time_diff_angles_traj_list_binstat))/np.mean(bindiffs[bindiffs>0]) )#*3
                
                print nbins
                print np.mean(bindiffs[bindiffs>0])
                
                # ~ sys.exit()
                
                
                avg_angles_subseg_fct_time_list_count    , _, _       = stats.binned_statistic(time_diff_angles_traj_list_binstat, angles_subseg_fct_time_list_binstat, 'count', bins=nbins)  
                avg_angles_subseg_fct_time_list    , _, _       = stats.binned_statistic(time_diff_angles_traj_list_binstat, angles_subseg_fct_time_list_binstat, 'mean', bins=nbins)  
                std_angles_subseg_fct_time_list    , _, _       = stats.binned_statistic(time_diff_angles_traj_list_binstat, angles_subseg_fct_time_list_binstat, mystderr, bins=nbins)  
                avg_time_diff_angles_traj_list         , _, _       = stats.binned_statistic(time_diff_angles_traj_list_binstat, time_diff_angles_traj_list_binstat, 'mean', bins=nbins)  
                
                print time_diff_angles_traj_list_binstat
                print angles_subseg_fct_time_list_binstat
                
                print avg_angles_subseg_fct_time_list_count
                print avg_time_diff_angles_traj_list
                print bindiffs
                # ~ print np.mean(bindiffs[np.abs(bindiffs)<50000]), np.mean(bindiffs), 
                
                # ~ sys.exit()
                
                nbins=int((np.amax(time_diff_traj_list_binstat)  - np.amin(time_diff_traj_list_binstat))/(avg_time_diff_traj_list[1] - avg_time_diff_traj_list[0] ) )#*3
                print nbins
                
                avg_time_diff_o_ttot_traj_list  , _, _       = stats.binned_statistic(time_diff_o_ttot_traj_list_binstat, time_diff_o_ttot_traj_list_binstat, 'mean', bins=nbins)   
                avg_end2end_coarse_gr_fct_time_list_bin_l_o_ltot    , _, _       = stats.binned_statistic(time_diff_o_ttot_traj_list_binstat, end2end_coarse_gr_fct_time_list_binstat, 'mean', bins=nbins) 
                
                
                
                
                nbins=int((np.amax(tlags_autocorr_list_binstat) - np.amin(tlags_autocorr_list_binstat))/2000 )
                
                print nbins
                
                avg_angles_autocorr_list_count    , _, _       = stats.binned_statistic(tlags_autocorr_list_binstat, angles_autocorr_list_binstat, 'count', bins=nbins)  
                avg_angles_autocorr_list    , _, _       = stats.binned_statistic(tlags_autocorr_list_binstat, angles_autocorr_list_binstat, 'mean', bins=nbins)  
                avg_vel_autocorr_list    , _, _       = stats.binned_statistic(tlags_autocorr_list_binstat, vel_autocorr_list_binstat, 'mean', bins=nbins)  
                avg_tlags_autocorr_list    , _, _       = stats.binned_statistic(tlags_autocorr_list_binstat, tlags_autocorr_list_binstat, 'mean', bins=nbins)  
                
                
                
                nbins=int((np.amax(tlags_autocorr_raw_list_binstat) - np.amin(tlags_autocorr_raw_list_binstat))/2000 )
                
                print nbins
                
                avg_angles_autocorr_raw_list_count    , _, _       = stats.binned_statistic(tlags_autocorr_raw_list_binstat, angles_autocorr_raw_list_binstat, 'count', bins=nbins)  
                avg_angles_autocorr_raw_list    , _, _       = stats.binned_statistic(tlags_autocorr_raw_list_binstat, angles_autocorr_raw_list_binstat, 'mean', bins=nbins)  
                avg_tlags_autocorr_raw_list    , _, _       = stats.binned_statistic(tlags_autocorr_raw_list_binstat, tlags_autocorr_raw_list_binstat, 'mean', bins=nbins)  
                
                
                
                
                
                
                avg_end2end_coarse_gr_list_count           = np.array([0.]) 
                std_end2end_coarse_gr_list           = np.array([0.]) 
                avg_end2end_coarse_gr_list          = np.array([0.]) 
                avg_path_l_serie_list            = np.array([0.]) 
                
                
                avg_path_l_o_ltot_serie_list      = np.array([0.]) 
                avg_end2end_coarse_gr_list_bin_l_o_ltot         = np.array([0.]) 
                        
                        
                        
                        

        else:
		    
            
            avg_angles_subseg_fct_time_list_count         = np.array([0.])  
            avg_angles_subseg_fct_time_list        = np.array([0.]) 
            std_angles_subseg_fct_time_list          = np.array([0.]) 
            avg_time_diff_angles_traj_list               = snp.array([0.]) 
            
            avg_angles_autocorr_list_count           = np.array([0.]) 
            avg_angles_autocorr_list              = np.array([0.]) 
            avg_vel_autocorr_list              = np.array([0.]) 
            avg_tlags_autocorr_list              = np.array([0.]) 
            avg_angles_autocorr_raw_list_count           = np.array([0.]) 
            avg_angles_autocorr_raw_list             = np.array([0.]) 
            avg_tlags_autocorr_raw_list              = np.array([0.]) 
            
            std_end2end_coarse_gr_fct_time_list           = np.array([0.]) 
            avg_end2end_coarse_gr_fct_time_list           = np.array([0.]) 
            avg_time_diff_traj_list                       = np.array([0.]) 
            
            avg_end2end_coarse_gr_list_count           = np.array([0.]) 
            std_end2end_coarse_gr_list           = np.array([0.]) 
            avg_end2end_coarse_gr_list           = np.array([0.]) 
            avg_path_l_serie_list                = np.array([0.]) 
            avg_path_l_o_ltot_serie_list         = np.array([0.]) 
            avg_end2end_coarse_gr_o_path_l_list  = np.array([0.]) 
            avg_end2end_coarse_gr_list_bin_l_o_ltot  = np.array([0.])
              
            
            
            avg_end2end_coarse_gr_fct_time_list_count  = np.array([0.])
            
            avg_time_diff_o_ttot_traj_list = np.array([0.])
            avg_end2end_coarse_gr_fct_time_list_bin_l_o_ltot = np.array([0.])
            
        
	    
	    
    else:
        
        
        
            
        
        file_out='{inp}/space_clustering_pers_l_real_{real}.dat'.format(inp=dir_in, real=real)
        
        #np.savez_compressed(file_out, hist_end2end_trajclade = hist_end2end_trajclade, hist_length_trajclade = hist_length_trajclade, hist_end2end_o_length_trajclade = hist_end2end_o_length_trajclade, avg_end2end_coarse_gr_list = avg_end2end_coarse_gr_list, avg_end2end_coarse_gr_list_count = avg_end2end_coarse_gr_list_count, avg_path_l_serie_list = avg_path_l_serie_list, avg_path_l_o_ltot_serie_list = avg_path_l_o_ltot_serie_list, avg_end2end_coarse_gr_list_bin_l_o_ltot = avg_end2end_coarse_gr_list_bin_l_o_ltot, avg_end2end_coarse_gr_fct_time_list_count = avg_end2end_coarse_gr_fct_time_list_count, avg_end2end_coarse_gr_fct_time_list = avg_end2end_coarse_gr_fct_time_list, avg_time_diff_traj_list = avg_time_diff_traj_list, avg_time_diff_o_ttot_traj_list = avg_time_diff_o_ttot_traj_list, avg_end2end_coarse_gr_fct_time_list_bin_l_o_ltot = avg_end2end_coarse_gr_fct_time_list_bin_l_o_ltot, hist_tdiff_trajclade = hist_tdiff_trajclade, path_l_serie_list_binstat = path_l_serie_list_binstat, end2end_coarse_gr_list_binstat = end2end_coarse_gr_list_binstat, path_l_o_ltot_serie_list_binstat = path_l_o_ltot_serie_list_binstat, time_diff_traj_list_binstat = time_diff_traj_list_binstat, end2end_coarse_gr_fct_time_list_binstat = end2end_coarse_gr_fct_time_list_binstat, time_diff_o_ttot_traj_list_binstat = time_diff_o_ttot_traj_list_binstat, hist_nofull_time_turns = hist_nofull_time_turns, hist_nofull_time_turns_o_ttot_traj = hist_nofull_time_turns_o_ttot_traj, hist_time_turns = hist_time_turns, hist_time_turns_o_ttot_traj = hist_time_turns_o_ttot_traj, hist_num_turns_o_tclade = hist_num_turns_o_tclade, hist_num_turns = hist_num_turns)#, hist_time_trajclade_turndet_reals = hist_time_trajclade_turndet_reals, hist_num_turns_wsplit_reals = hist_num_turns_wsplit_reals
        
        
        
        z = np.load(file_out)
        
        avg_angles_subseg_fct_time_list_count = z['avg_angles_subseg_fct_time_list_count']
        avg_angles_subseg_fct_time_list       = z['avg_angles_subseg_fct_time_list      ']
        std_angles_subseg_fct_time_list       = z['std_angles_subseg_fct_time_list      ']
        avg_time_diff_angles_traj_list        = z['avg_time_diff_angles_traj_list       ']
        
        avg_end2end_coarse_gr_list_count = z['avg_end2end_coarse_gr_list_count']
        std_end2end_coarse_gr_list = z['std_end2end_coarse_gr_list']
        avg_end2end_coarse_gr_list = z['avg_end2end_coarse_gr_list']
        avg_path_l_serie_list = z['avg_path_l_serie_list']
        avg_path_l_o_ltot_serie_list = z['avg_path_l_o_ltot_serie_list']
        avg_end2end_coarse_gr_list_bin_l_o_ltot = z['avg_end2end_coarse_gr_list_bin_l_o_ltot']
        std_end2end_coarse_gr_fct_time_list = z['std_end2end_coarse_gr_fct_time_list']
        avg_end2end_coarse_gr_fct_time_list = z['avg_end2end_coarse_gr_fct_time_list']
        avg_end2end_coarse_gr_fct_time_list_count = z['avg_end2end_coarse_gr_fct_time_list_count']
        avg_time_diff_traj_list = z['avg_time_diff_traj_list']
        
        avg_angles_autocorr_list_count = z['avg_angles_autocorr_list_count']
        avg_angles_autocorr_list       = z['avg_angles_autocorr_list      ']
        avg_vel_autocorr_list          = z['avg_vel_autocorr_list         ']
        avg_tlags_autocorr_list        = z['avg_tlags_autocorr_list       ']
        avg_angles_autocorr_raw_list_count   = z['avg_angles_autocorr_raw_list  ']
        avg_angles_autocorr_raw_list   = z['avg_angles_autocorr_raw_list  ']
        avg_tlags_autocorr_raw_list    = z['avg_tlags_autocorr_raw_list   ']
        
    
    

    ## plot  avg_end2end_coarse_gr_list
    
    fig = plt.figure(figsize=thisfigsize)
    grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.22,
                 wspace=0.4, hspace=0.35) 
    labeled_axes = []
    ax = plt.Subplot(fig, grid[0, 0])
    fig.add_subplot(ax)
    labeled_axes.append(ax)
    ax.plot(np.arange(len(angles_vs_time_raw_list_binstat)), angles_vs_time_raw_list_binstat, linestyle='-', color='g', label="noisy data")
    #ax.plot(avg_path_l_serie_list, fit_spline, linestyle='-', color='r', label="fitted spline")
    ax.plot(np.arange(len(angles_vs_time_list_binstat)), angles_vs_time_list_binstat, linestyle='-', color='b', label="smooth")
    ax.set_xlabel('time ')
    ax.set_ylabel('angle')
    ax.xaxis.labelpad = axis_labelpad
    ax.yaxis.labelpad = axis_labelpad
    ax.legend(frameon=False, ncol=1, columnspacing=0.5, handletextpad=0.2,
          loc='upper right', bbox_to_anchor=(2.15, 1.18))
    mpsetup.despine(ax) 
    
    
    #### finish figure ####
    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
             xycoords=('axes fraction'), fontweight = 'bold')
#    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
    out_file='{out}/split_traj_angle_fct_time_{real}.png'.format(out=dir_out_plots, real=real)
#    print out_file
    fig.savefig(out_file)
    


    ## plot  avg_end2end_coarse_gr_list
    
    fig = plt.figure(figsize=thisfigsize)
    grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.22,
                 wspace=0.4, hspace=0.35) 
    labeled_axes = []
    ax = plt.Subplot(fig, grid[0, 0])
    fig.add_subplot(ax)
    labeled_axes.append(ax)
    #ax.plot(avg_path_l_serie_list, fit_spline, linestyle='-', color='r', label="fitted spline")
    ax.plot(np.arange(len(angles_vs_time_list_binstat) - 1), np.asarray(angles_vs_time_list_binstat)[1:] - np.asarray(angles_vs_time_list_binstat)[:-1], linestyle='-', color='b', label="smooth")
    ax.set_xlabel('time ')
    ax.set_ylabel('angle')
    ax.xaxis.labelpad = axis_labelpad
    ax.yaxis.labelpad = axis_labelpad
    ax.legend(frameon=False, ncol=1, columnspacing=0.5, handletextpad=0.2,
          loc='upper right', bbox_to_anchor=(2.15, 1.18))
    mpsetup.despine(ax) 
    
    
    #### finish figure ####
    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
             xycoords=('axes fraction'), fontweight = 'bold')
#    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
    out_file='{out}/split_traj_angle_diff_fct_time_{real}.png'.format(out=dir_out_plots, real=real)
#    print out_file
    fig.savefig(out_file)
    



    
    fig.clf()
    plt.close('all')
    
        
        
        
        
        
    
    file_out='{inp}/space_clustering_pers_l_real_{real}.dat'.format(inp=dir_in, real=real)
    
    np.savez_compressed(file_out, hist_end2end_trajclade = hist_end2end_trajclade, hist_length_trajclade = hist_length_trajclade, hist_end2end_o_length_trajclade = hist_end2end_o_length_trajclade, avg_end2end_coarse_gr_list = avg_end2end_coarse_gr_list, avg_end2end_coarse_gr_list_count = avg_end2end_coarse_gr_list_count, avg_path_l_serie_list = avg_path_l_serie_list, avg_path_l_o_ltot_serie_list = avg_path_l_o_ltot_serie_list, avg_end2end_coarse_gr_list_bin_l_o_ltot = avg_end2end_coarse_gr_list_bin_l_o_ltot, avg_end2end_coarse_gr_fct_time_list_count = avg_end2end_coarse_gr_fct_time_list_count, avg_end2end_coarse_gr_fct_time_list = avg_end2end_coarse_gr_fct_time_list, avg_time_diff_traj_list = avg_time_diff_traj_list, avg_time_diff_o_ttot_traj_list = avg_time_diff_o_ttot_traj_list, avg_end2end_coarse_gr_fct_time_list_bin_l_o_ltot = avg_end2end_coarse_gr_fct_time_list_bin_l_o_ltot, hist_tdiff_trajclade = hist_tdiff_trajclade, path_l_serie_list_binstat = path_l_serie_list_binstat, end2end_coarse_gr_list_binstat = end2end_coarse_gr_list_binstat, path_l_o_ltot_serie_list_binstat = path_l_o_ltot_serie_list_binstat, time_diff_traj_list_binstat = time_diff_traj_list_binstat, end2end_coarse_gr_fct_time_list_binstat = end2end_coarse_gr_fct_time_list_binstat, time_diff_o_ttot_traj_list_binstat = time_diff_o_ttot_traj_list_binstat, hist_nofull_time_turns = hist_nofull_time_turns, hist_nofull_time_turns_o_ttot_traj = hist_nofull_time_turns_o_ttot_traj, hist_time_turns = hist_time_turns, hist_time_turns_o_ttot_traj = hist_time_turns_o_ttot_traj, hist_num_turns_o_tclade = hist_num_turns_o_tclade, hist_num_turns = hist_num_turns, avg_angles_subseg_fct_time_list = avg_angles_subseg_fct_time_list, avg_angles_subseg_fct_time_list_count = avg_angles_subseg_fct_time_list_count, std_angles_subseg_fct_time_list = std_angles_subseg_fct_time_list, avg_time_diff_angles_traj_list = avg_time_diff_angles_traj_list, angles_subseg_fct_time_list_binstat = angles_subseg_fct_time_list_binstat, angles_vs_time_raw_list_binstat = angles_vs_time_raw_list_binstat, angles_vs_time_list_binstat = angles_vs_time_list_binstat, avg_angles_autocorr_list_count = avg_angles_autocorr_list_count, avg_angles_autocorr_list = avg_angles_autocorr_list, avg_vel_autocorr_list = avg_vel_autocorr_list, avg_tlags_autocorr_list = avg_tlags_autocorr_list, avg_angles_autocorr_raw_list_count = avg_angles_autocorr_raw_list_count, avg_angles_autocorr_raw_list = avg_angles_autocorr_raw_list, avg_tlags_autocorr_raw_list = avg_tlags_autocorr_raw_list, angles_autocorr_list_binstat = angles_autocorr_list_binstat, vel_autocorr_list_binstat = vel_autocorr_list_binstat, tlags_autocorr_list_binstat = tlags_autocorr_list_binstat, angles_autocorr_raw_list_binstat = angles_autocorr_raw_list_binstat, tlags_autocorr_raw_list_binstat = tlags_autocorr_raw_list_binstat, std_end2end_coarse_gr_list = std_end2end_coarse_gr_list, std_end2end_coarse_gr_fct_time_list = std_end2end_coarse_gr_fct_time_list)#, hist_time_trajclade_turndet_reals = hist_time_trajclade_turndet_reals, hist_num_turns_wsplit_reals = hist_num_turns_wsplit_reals
    
    
    
     





    print 'FITTING' # fit MSD to appropriate function
          
    
            
            
    #data = np.array([times, n_clusters_dbscan_list, n_clusters_dbscan_est_list, n_clusters_dbscan_opt_kvar_list, times, eps_dbscan_est_list, eps_dbscan_opt_kvar_list,  size_clusters_dbscan_list, size_clusters_dbscan_est_list, size_clusters_dbscan_opt_kvar_list, density_clusters_dbscan_list, density_clusters_dbscan_est_list, density_clusters_dbscan_opt_kvar_list, n_clusters_dbscan_opt_dispersion_list, n_clusters_dbscan_opt_dispersion_CV_list, n_clusters_dbscan_opt_kCV_list, mean_dispersion_dbscan_list, mean_dispersion_dbscan_est_list, mean_dispersion_dbscan_opt_list, mean_dispersion_CV_dbscan_list, mean_dispersion_CV_dbscan_est_list, mean_dispersion_CV_dbscan_opt_list, mean_kCV_dbscan_list, mean_kCV_dbscan_est_list, mean_kCV_dbscan_opt_list, mean_kvar_dbscan_list, mean_kvar_dbscan_est_list, mean_kvar_dbscan_opt_list, inter_clusters_dist_dbscan_opt_kvar_list, size_clusters_split, inter_clusters_dist_split, density_clusters_split, n_clusters_split, num_vir_clusters_split, num_ext_tot, num_split_tot, num_ext_per_cl_tot, num_split_per_cl_tot, num_vir_o_nclust_tot_list, vel_clust_split_alltimes, var_parall_clust_split_alltimes, var_perp_clust_split_alltimes, var_tot_clust_split_alltimes])
            
        
    file_out='{inp}/space_clustering_fct_time_real_{real}.dat.npz'.format(inp=dir_in, real=real)
    
    #np.savez_compressed(file_out,data = data, frac_1cl=frac_1cl)
            
               
    #~ np_load_old = np.load

    #~ # modify the default parameters of np.load
    #~ np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)
   
    z = np.load(file_out, allow_pickle=True)
    #~ np.load = np_load_old
    
    
    #~ data_clust=np.asarray(z['data'])
    data_clust=z['data_allclust'] 
    
    print [a.size for a in data_clust]
    
    print data_clust
    print data_clust.shape
    
    if data_clust.ndim>1:
     
        frac_1cl=z['frac_1cl']
    
        #if data_clust.shape[1]>55:
        
        
        vel_clust_split_alltimes        =  data_clust[0,:]
        var_parall_clust_split_alltimes =  data_clust[1,:]
        var_perp_clust_split_alltimes   =  data_clust[2,:]
        var_tot_clust_split_alltimes    =  data_clust[3,:]
        
        
        
        avg_vel_clust_split_alltimes=vel_clust_split_alltimes.mean()
        var_vel_clust_split_alltimes=vel_clust_split_alltimes.var()
        
        
        avg_var_parall_clust_split_alltimes=var_parall_clust_split_alltimes.mean()
        
        
        avg_var_perp_clust_split_alltimes=var_perp_clust_split_alltimes.mean()
        
        
        avg_var_tot_clust_split_alltimes=var_tot_clust_split_alltimes.mean()
        
            
            
    
    
    avg_frac_1cl                        = frac_1cl.mean()
    #avg_hist_avg_end2end_o_length_trajclade = hist_avg_end2end_o_length_trajclade.mean()
    avg_hist_length_trajclade = hist_length_trajclade.mean()
    avg_hist_length_sq_trajclade = (hist_length_trajclade**2).mean()
    
    #avg_hist_tdiff_trajclade = hist_tdiff_trajclade.mean()
    
    print avg_vel_clust_split_alltimes, var_vel_clust_split_alltimes, var_vel_clust_split_alltimes/(avg_vel_clust_split_alltimes**2)



    
    
    
    
    if not traj_anal:
        
        
        
        std_end2end_coarse_gr_list = z['std_end2end_coarse_gr_list']
        avg_end2end_coarse_gr_list = z['avg_end2end_coarse_gr_list']
        avg_path_l_serie_list = z['avg_path_l_serie_list']
        avg_path_l_o_ltot_serie_list = z['avg_path_l_o_ltot_serie_list']
        avg_end2end_coarse_gr_list_bin_l_o_ltot = z['avg_end2end_coarse_gr_list_bin_l_o_ltot']
        
        
        avg_end2end_coarse_gr_list_count = z['avg_end2end_coarse_gr_list_count']
        
     






    
       
    
    
    avg_angles_autocorr_list_orig       = avg_angles_autocorr_list    .copy() 
    avg_vel_autocorr_list_orig       = avg_vel_autocorr_list       .copy() 
    avg_angles_autocorr_list_count_orig       = avg_angles_autocorr_list_count.copy() 
    avg_tlags_autocorr_list_orig       = avg_tlags_autocorr_list .copy() 
    avg_angles_autocorr_raw_list_orig       = avg_angles_autocorr_raw_list.copy() 
    avg_tlags_autocorr_raw_list_orig       = avg_tlags_autocorr_raw_list .copy() 
       
    avg_angles_subseg_fct_time_list_orig       = avg_angles_subseg_fct_time_list      .copy()    
    std_angles_subseg_fct_time_list_orig       = std_angles_subseg_fct_time_list      .copy()    
    avg_time_diff_angles_traj_list_orig        = avg_time_diff_angles_traj_list       .copy()    
    avg_angles_subseg_fct_time_list_count_orig       = avg_angles_subseg_fct_time_list_count.copy()   
     
    std_end2end_coarse_gr_fct_time_list_orig       = std_end2end_coarse_gr_fct_time_list.copy()    
    avg_end2end_coarse_gr_fct_time_list_orig       = avg_end2end_coarse_gr_fct_time_list.copy()    
    avg_time_diff_traj_list_orig                   = avg_time_diff_traj_list.copy()                
                               
                               
                               
    avg_end2end_coarse_gr_list_count_orig                = avg_end2end_coarse_gr_list_count.copy()             
    std_end2end_coarse_gr_list_orig                = std_end2end_coarse_gr_list.copy()             
    avg_end2end_coarse_gr_list_orig                = avg_end2end_coarse_gr_list.copy()             
    avg_path_l_serie_list_orig                     = avg_path_l_serie_list.copy()                  
    avg_path_l_o_ltot_serie_list_orig              = avg_path_l_o_ltot_serie_list.copy()           
    #avg_end2end_coarse_gr_o_path_l_list_orig       = avg_end2end_coarse_gr_o_path_l_list.copy()      
                                                     
    avg_end2end_coarse_gr_list_bin_l_o_ltot_orig   = avg_end2end_coarse_gr_list_bin_l_o_ltot.copy()
    
    
    
    
    avg_angles_autocorr_list      = avg_angles_autocorr_list[(np.isfinite(avg_angles_autocorr_list_orig)) & (np.isfinite(avg_tlags_autocorr_list_orig)) & (avg_tlags_autocorr_list_orig>0)]     
    avg_vel_autocorr_list         = avg_vel_autocorr_list[(np.isfinite(avg_vel_autocorr_list_orig)) & (np.isfinite(avg_tlags_autocorr_list_orig)) & (avg_tlags_autocorr_list_orig>0)]     
    avg_angles_autocorr_list_count       = avg_angles_autocorr_list_count[(np.isfinite(avg_angles_autocorr_list_count_orig)) & (np.isfinite(avg_tlags_autocorr_list_orig)) & (avg_tlags_autocorr_list_orig>0)]     
    avg_tlags_autocorr_list       = avg_tlags_autocorr_list[(np.isfinite(avg_tlags_autocorr_list_orig)) & (np.isfinite(avg_tlags_autocorr_list_orig)) & (avg_tlags_autocorr_list_orig>0)]     
    avg_angles_autocorr_raw_list  = avg_angles_autocorr_raw_list[(np.isfinite(avg_angles_autocorr_raw_list_orig)) & (np.isfinite(avg_tlags_autocorr_raw_list_orig)) & (avg_tlags_autocorr_raw_list_orig>0)]     
    avg_tlags_autocorr_raw_list   = avg_tlags_autocorr_raw_list[(np.isfinite(avg_angles_autocorr_raw_list_orig)) & (np.isfinite(avg_tlags_autocorr_raw_list_orig)) & (avg_tlags_autocorr_raw_list_orig>0)]     
    
    # ~ print avg_angles_autocorr_raw_list_orig.size
    # ~ print avg_tlags_autocorr_raw_list_orig.size
    
    # ~ print avg_angles_autocorr_raw_list.size
    # ~ print avg_tlags_autocorr_raw_list.size
    
    # ~ sys.exit()
    
         
    avg_end2end_coarse_gr_list_count             = avg_end2end_coarse_gr_list_count_orig[(np.isfinite(avg_end2end_coarse_gr_list_count)) & (np.isfinite(avg_path_l_serie_list_orig)) & (avg_path_l_serie_list_orig>0)]          
    std_end2end_coarse_gr_list             = std_end2end_coarse_gr_list_orig[(np.isfinite(std_end2end_coarse_gr_list_orig)) & (np.isfinite(avg_path_l_serie_list_orig)) & (avg_path_l_serie_list_orig>0)]          
    avg_end2end_coarse_gr_list             = avg_end2end_coarse_gr_list_orig[(np.isfinite(avg_end2end_coarse_gr_list_orig)) & (np.isfinite(avg_path_l_serie_list_orig)) & (avg_path_l_serie_list_orig>0)]          
    avg_path_l_serie_list                  = avg_path_l_serie_list_orig[(np.isfinite(avg_end2end_coarse_gr_list_orig)) & (np.isfinite(avg_path_l_serie_list_orig)) & (avg_path_l_serie_list_orig>0)]
              
    avg_path_l_o_ltot_serie_list           = avg_path_l_o_ltot_serie_list_orig[(np.isfinite(avg_end2end_coarse_gr_list_bin_l_o_ltot_orig)) & (np.isfinite(avg_path_l_o_ltot_serie_list_orig)) & (avg_path_l_o_ltot_serie_list_orig>0)]          
    avg_end2end_coarse_gr_list_bin_l_o_ltot= avg_end2end_coarse_gr_list_bin_l_o_ltot_orig[(np.isfinite(avg_end2end_coarse_gr_list_bin_l_o_ltot_orig)) & (np.isfinite(avg_path_l_o_ltot_serie_list_orig)) & (avg_path_l_o_ltot_serie_list_orig>0)] 
    
    
    ####
    
    # ~ time_thr=10000
    
    if time_thr < 2*window_typ*np.mean(avg_time_diff_angles_traj_list[1:] - avg_time_diff_angles_traj_list[:-1]):
        time_thr = 2*window_typ*np.mean(avg_time_diff_angles_traj_list[1:] - avg_time_diff_angles_traj_list[:-1])
  
    
    std_angles_subseg_fct_time_list= std_angles_subseg_fct_time_list_orig[(np.isfinite(std_angles_subseg_fct_time_list_orig)) & (np.isfinite(avg_time_diff_angles_traj_list_orig)) & (avg_time_diff_angles_traj_list_orig>time_thr)]
    
    avg_angles_subseg_fct_time_list_count= avg_angles_subseg_fct_time_list_count[(np.isfinite(avg_angles_subseg_fct_time_list_count_orig)) & (np.isfinite(avg_time_diff_angles_traj_list_orig)) & (avg_time_diff_angles_traj_list_orig>time_thr)]
    
    avg_angles_subseg_fct_time_list= avg_angles_subseg_fct_time_list_orig[(np.isfinite(avg_angles_subseg_fct_time_list_orig)) & (np.isfinite(avg_time_diff_angles_traj_list_orig)) & (avg_time_diff_angles_traj_list_orig>time_thr)]
    
    avg_time_diff_angles_traj_list= avg_time_diff_angles_traj_list[(np.isfinite(avg_time_diff_angles_traj_list_orig)) & (np.isfinite(avg_time_diff_angles_traj_list_orig)) & (avg_time_diff_angles_traj_list_orig>time_thr)]
    
    print avg_angles_subseg_fct_time_list_count_orig
    print avg_angles_subseg_fct_time_list_count
    print avg_angles_subseg_fct_time_list_count_orig.size, avg_angles_subseg_fct_time_list_count.size
    print time_thr
    # ~ print np.mean(bindiffs[bindiffs>0]), np.mean(avg_time_diff_angles_traj_list[1:] - avg_time_diff_angles_traj_list[:-1])
    # ~ sys.exit()

    
    #################
    
    
            
            
    print 'FILTER'
    # filter few traces
    #avg_avg_end2end_coarse_gr_fct_time_list_count
    print avg_end2end_coarse_gr_fct_time_list_count
    print avg_end2end_coarse_gr_fct_time_list_count.shape
    # this filters out those realizationsa with too little time, at least 3 realizarions of 100k cycles, or 3 shorter clades
    
    
    # ~ mask_time=(avg_angles_subseg_fct_time_list_count > np.amax(avg_angles_subseg_fct_time_list_count)/100.) & (avg_angles_subseg_fct_time_list_count > 8)
    mask_time=(avg_angles_subseg_fct_time_list_count > np.amax(avg_angles_subseg_fct_time_list_count) - np.amax(avg_angles_subseg_fct_time_list_count)/3 ) & (avg_angles_subseg_fct_time_list_count > 9)

    avg_angles_subseg_fct_time_list        = avg_angles_subseg_fct_time_list      [mask_time]         # [:-1]
    std_angles_subseg_fct_time_list        = std_angles_subseg_fct_time_list      [mask_time]         # [:-1]
    avg_time_diff_angles_traj_list         = avg_time_diff_angles_traj_list [mask_time]         # [:-1]
    avg_angles_subseg_fct_time_list_count  = avg_angles_subseg_fct_time_list_count[mask_time]         # [:-1]   
    
    mask_time=(avg_end2end_coarse_gr_fct_time_list_count > np.amax(avg_end2end_coarse_gr_fct_time_list_count)/100.) & (avg_end2end_coarse_gr_fct_time_list_count > 2)
    
    std_end2end_coarse_gr_fct_time_list  = std_end2end_coarse_gr_fct_time_list[mask_time] # [:-1]     
    avg_end2end_coarse_gr_fct_time_list  = avg_end2end_coarse_gr_fct_time_list[mask_time] # [:-1]     
    avg_time_diff_traj_list  = avg_time_diff_traj_list[mask_time]                         # [:-1]
    
    avg_end2end_coarse_gr_fct_time_list_count = avg_end2end_coarse_gr_fct_time_list_count[mask_time]  
    
        
    print avg_end2end_coarse_gr_fct_time_list.shape
    
  
      
    print 'FILTERED'
    
    
    
    # fit with const vel ballistic RW as in Peruani, 2 parameters. Non const vel would require fitting 4 parameters 
    
    
    # OUTDATED
    
    pers_l_avg_end2end_fit              =0
    err_pers_l_avg_end2end_fit               =0
    chisq_avg_end2end_fit               =0
    pers_l_avg_end2end_bin_l_o_ltot_fit =0
    pers_l_avg_end2end_fit_varvel_fixvel   =0
    beta_avg_end2end_fit_varvel_fixvel     =0
    err_pers_l_avg_end2end_fit_varvel_fixvel =0
    err_beta_avg_end2end_fit_varvel_fixvel   =0
    chisq_avg_end2end_fit_varvel_fixvel   =0
    #print avg_path_l_serie_list
    #print avg_end2end_coarse_gr_list
    #
    #
    #print avg_path_l_o_ltot_serie_list
    #print avg_end2end_coarse_gr_list_bin_l_o_ltot
    
    
    
    
    # autocorr
    
    
    const_angles_autocorr_raw   = 0
    const_angles_autocorr   = 0
    const_vel_autocorr   = 0
    
     
    #print avg_time_diff_traj_list
    
    
    if avg_angles_autocorr_raw_list.size >10:
        
        print 'fit autocorr'
        
        try:
        
             # acf_decay(x, autocorr)
            sol_acf_decay_angles_autocorr_raw, pcov_acf_decay_angles_autocorr_raw = curve_fit(acf_decay, avg_tlags_autocorr_raw_list, avg_angles_autocorr_raw_list, p0 = np.asarray([10000]), bounds=(0.0000001, np.inf), max_nfev=10000)
        
            print sol_acf_decay_angles_autocorr_raw, pcov_acf_decay_angles_autocorr_raw
            const_angles_autocorr_raw         =    sol_acf_decay_angles_autocorr_raw[0]
            
        
             # acf_decay(x, autocorr)
            sol_acf_decay_angles_autocorr, pcov_acf_decay_angles_autocorr = curve_fit(acf_decay, avg_tlags_autocorr_list, avg_angles_autocorr_list, p0 = np.asarray([10000]), bounds=(0.0000001, np.inf), max_nfev=10000)
        
            print sol_acf_decay_angles_autocorr, pcov_acf_decay_angles_autocorr
            const_angles_autocorr         =    sol_acf_decay_angles_autocorr[0]
            
        
             # acf_decay(x, autocorr)
            sol_acf_decay_vel_autocorr, pcov_acf_decay_vel_autocorr = curve_fit(acf_decay, avg_tlags_autocorr_list, avg_vel_autocorr_list, p0 = np.asarray([100]), bounds=(0.0000001, np.inf), max_nfev=10000)
        
            print sol_acf_decay_vel_autocorr, pcov_acf_decay_vel_autocorr
            const_vel_autocorr         =    sol_acf_decay_vel_autocorr[0]
            
        
        
        except ValueError:
            print 'Error in fitting'		
        except RuntimeError:
            print 'Exceeded max iterations in fitting'		
    
    
    
    
    
    
    
    
    # fct time
    
    
    
    diff_avg_angles_subseg_fit_rw_fct_time   = 0
    err_diff_avg_angles_subseg_fit_rw_fct_time   = 0
    chisq_avg_angles_subseg_fit_rw_fct_time   = 0
    pers_l_avg_end2end_fit_constvel_fct_time_fixvel   = 0
    err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel   = 0
    chisq_avg_end2end_fit_constvel_fct_time_fixvel   = 0
    
    max_clade_time   = 0
    
    #print avg_time_diff_traj_list
    
    print avg_angles_subseg_fct_time_list.size
    
    # FIT LINEAR ANGLES MSD
    
    
    if avg_angles_subseg_fct_time_list.size >3:
        
        max_clade_time=avg_time_diff_traj_list_orig[-1]
        print 'fit time const vel'
        
        try:
        
            
            xData = avg_time_diff_angles_traj_list
            yData=avg_angles_subseg_fct_time_list
            ystd=np.where(std_angles_subseg_fct_time_list>0, std_angles_subseg_fct_time_list , np.amin(std_angles_subseg_fct_time_list[std_angles_subseg_fct_time_list>0])) 
            #~ print ystd     
             # rw_msd_fct_time(x, diff)
             
            func =rw_msd_fct_time
            # ~ sol, pcov = curve_fit(func, xData, yData, sigma = ystd, absolute_sigma = False, p0 = np.asarray([0.1]), bounds=(0., np.inf), max_nfev=10000, method='trf')#, loss='cauchy'
            sol, pcov = curve_fit(func, xData, yData, sigma = ystd, absolute_sigma = False, p0 = np.asarray([0.1, 0]), bounds=(0., np.inf), max_nfev=10000, method='trf')#, loss='cauchy'
            
            perr = np.sqrt(np.diag(pcov))
            
            for i in range(len(sol)):
                print sol[i], ' +- ', perr[i] 
                        
            #print sol_rw_msd_fct_time, pcov_rw_msd_fct_time
            diff_avg_angles_subseg_fit_rw_fct_time         =    sol[0]
            err_diff_avg_angles_subseg_fit_rw_fct_time         =     perr[0]
            print "pers time ", 2./diff_avg_angles_subseg_fit_rw_fct_time, ' +- ', (2./diff_avg_angles_subseg_fit_rw_fct_time)*err_diff_avg_angles_subseg_fit_rw_fct_time/diff_avg_angles_subseg_fit_rw_fct_time
            
                
            # prepare confidence level curves
            nstd = 2. # to draw 2-sigma intervals
            popt_up = sol + nstd * perr
            popt_dw = sol - nstd * perr
            
            fit_up_avg_angles_subseg_fit_rw_fct_time = func(xData, *popt_up) 
            fit_dw_avg_angles_subseg_fit_rw_fct_time = func(xData, *popt_dw)  
             
            
            
            # goodness of fit
            
            modelPredictions_avg_angles_subseg_fit_rw_fct_time = func(xData, *sol) 
            
            absError = modelPredictions_avg_angles_subseg_fit_rw_fct_time  - yData
            
            SE = np.square(absError) # squared errors
            MSE = np.mean(SE) # mean squared errors
            RMSE  = np.sqrt(MSE) # Root Mean Squared Error, RMSE
            Rsquared  = 1.0 - (np.var(absError) / np.var(yData))
            print 'RMSE: ', RMSE
            print 'R-squared: ', Rsquared
            
            chisq_avg_angles_subseg_fit_rw_fct_time = np.sum(SE/(ystd**2))/float(yData.size) # chisq
            print 'Chi-squared: ', chisq_avg_angles_subseg_fit_rw_fct_time
               
            
            #v=avg_vel_clust_split_alltimes

        except ValueError:
            print 'Error in fitting'		
        except RuntimeError:
            print 'Exceeded max iterations in fitting'		
    
        
        
     # FIT POSITION MSD FROM PERUANI PRL, ALTERNATIVE INFERENCE
    if avg_end2end_coarse_gr_fct_time_list.size >10:
        
        max_clade_time=avg_time_diff_traj_list_orig[-1]
        print 'fit time const vel'
        
        try:
        
        
        
           
            xData = avg_time_diff_traj_list
            yData=avg_end2end_coarse_gr_fct_time_list
            ystd=std_end2end_coarse_gr_fct_time_list   
            
            #~ xData = time_diff_traj_list_binstat
            #~ yData=end2end_coarse_gr_fct_time_list_binstat
            
            #~ func =  lambda x, p_l: ballistic_msd_constvel_fct_time(x, p_l, avg_vel_clust_split_alltimes)
            #~ p0 = np.asarray([100.])

            func =  ballistic_msd_constvel_fct_time # also fit speed
            
            
            p0 = np.asarray([2./diff_avg_angles_subseg_fit_rw_fct_time,avg_vel_clust_split_alltimes])
            
            if diff_avg_angles_subseg_fit_rw_fct_time==0:
                p0 = np.asarray([1000000,avg_vel_clust_split_alltimes])

            sol, pcov = curve_fit( func, xData, yData, sigma = ystd, absolute_sigma = False, p0 = p0, bounds=(0.00000001, np.inf), max_nfev=10000, method='trf')#, f_scale=np.mean(ystd), method='trf', loss='cauchy', loss='soft_l1'

        
            perr = np.sqrt(np.diag(pcov))
            
            for i in range(len(sol)):
                print sol[i], ' +- ', perr[i] 
            
            #print sol_ballistic_msd_constvel_fct_time_fixvel, pcov_ballistic_msd_constvel_fct_time_fixvel
            
            pers_l_avg_end2end_fit_constvel_fct_time_fixvel          =    sol[0]
            err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel          =    perr[0]
            
          
            # prepare confidence level curves
            nstd = 2. # to draw 2-sigma intervals
            popt_up = sol + nstd * perr
            popt_dw = sol - nstd * perr
    
            #~ fit_up_avg_end2end_fit_constvel_fct_time_fixvel = ballistic_msd_constvel_fct_time(xData, popt_up[0], avg_vel_clust_split_alltimes) 
            #~ fit_dw_avg_end2end_fit_constvel_fct_time_fixvel = ballistic_msd_constvel_fct_time(xData, popt_dw[0], avg_vel_clust_split_alltimes) 
             
    
            fit_up_avg_end2end_fit_constvel_fct_time_fixvel = func(xData, *popt_up) 
            fit_dw_avg_end2end_fit_constvel_fct_time_fixvel = func(xData, *popt_dw) 
             
            
            
            # goodness of fit
            
            #~ modelPredictions_avg_end2end_fit_constvel_fct_time_fixvel = ballistic_msd_constvel_fct_time(xData, pers_l_avg_end2end_fit_constvel_fct_time_fixvel, avg_vel_clust_split_alltimes)
            modelPredictions_avg_end2end_fit_constvel_fct_time_fixvel = func(xData, *sol) 
            
            absError = modelPredictions_avg_end2end_fit_constvel_fct_time_fixvel  - yData
            
            SE = np.square(absError) # squared errors
            MSE = np.mean(SE) # mean squared errors
            RMSE  = np.sqrt(MSE) # Root Mean Squared Error, RMSE
            Rsquared  = 1.0 - (np.var(absError) / np.var(yData))
            print 'RMSE: ', RMSE
            print 'R-squared: ', Rsquared
            
            chisq_avg_end2end_fit_constvel_fct_time_fixvel = np.sum(SE/(ystd**2))/float(yData.size) # chisq
            print 'Chi-squared: ', chisq_avg_end2end_fit_constvel_fct_time_fixvel
               
        
        except ValueError:
            print 'Error in fitting'		
        except RuntimeError:
            print 'Exceeded max iterations in fitting'		
    
    
    
    
    
    
    pers_l_avg_end2end_fit_varvel_fct_time_fixvel   = 0
    beta_avg_end2end_fit_varvel_fct_time_fixvel= 0
    err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel= 0
    err_beta_avg_end2end_fit_varvel_fct_time_fixvel  = 0
    chisq_avg_end2end_fit_varvel_fct_time_fixvel  = 0
    
    # FIT POSITION MSD ACCOUNTING FOR SPEED FLUCTUATIONS FROM PERUANI PRL, OTHER ALTERNATIVE INFERENCE
    if avg_end2end_coarse_gr_fct_time_list.size >10:
        print 'fit time var vel'
        
        try:
            
            #sol_ballistic_msd_varvel_fct_time_fixvel, pcov_ballistic_msd_varvel_fct_time_fixvel = curve_fit( lambda x, p_l: ballistic_msd_varvel_fct_time_nobeta(x, p_l, avg_vel_clust_split_alltimes, p_l, var_vel_clust_split_alltimes), avg_time_diff_traj_list, avg_end2end_coarse_gr_fct_time_list, p0 = np.asarray([100.]), bounds=(0.0001, np.inf), max_nfev=10000)
            
            
            xData = avg_time_diff_traj_list
            yData=avg_end2end_coarse_gr_fct_time_list
            ystd=std_end2end_coarse_gr_fct_time_list   
            
            #~ func= lambda x, p_l, beta: ballistic_msd_varvel_fct_time(x, p_l, avg_vel_clust_split_alltimes, beta, var_vel_clust_split_alltimes)
            #~ p0 = np.asarray([10000., 0.001])
            
            func= lambda x, p_l, v, beta: ballistic_msd_varvel_fct_time(x, p_l, v, beta, var_vel_clust_split_alltimes)
            
            p0 = np.asarray([2./diff_avg_angles_subseg_fit_rw_fct_time, avg_vel_clust_split_alltimes, 0.001])
            
            if diff_avg_angles_subseg_fit_rw_fct_time==0:
                p0 = np.asarray([100000., avg_vel_clust_split_alltimes, 0.001])
            
            #~ func= ballistic_msd_varvel_fct_time
            #~ p0 = np.asarray([10000., avg_vel_clust_split_alltimes, 0.001, var_vel_clust_split_alltimes])
            
            sol, pcov = curve_fit( func, xData, yData, sigma = ystd, absolute_sigma = False, p0 = p0, bounds=(0.0000001, np.inf), max_nfev=10000, method='trf') # , loss='soft_l1', loss='cauchy'
        
        
            perr = np.sqrt(np.diag(pcov))
            
            for i in range(len(sol)):
                print sol[i], ' +- ', perr[i] 
            
            #print sol_ballistic_msd_varvel_fct_time_fixvel, pcov_ballistic_msd_varvel_fct_time_fixvel
            
            pers_l_avg_end2end_fit_varvel_fct_time_fixvel          =    sol[0]
            beta_avg_end2end_fit_varvel_fct_time_fixvel        =    sol[1]
            err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel          =  perr[0]
            err_beta_avg_end2end_fit_varvel_fct_time_fixvel          =    perr[1]
            
            
            #velcorr_l_avg_end2end_fit_varvel_fct_time_fixvel       =    sol_ballistic_msd_varvel_fct_time_fixvel[1] - sol_ballistic_msd_varvel_fct_time_fixvel[0]
            
            
            
          
            # prepare confidence level curves
            
            perr[1]/=1000 # hack to keep cl narrow in the plot.
            nstd = 2. # to draw 2-sigma intervals
            popt_up = sol + nstd * perr
            popt_dw = np.maximum(np.zeros_like(sol) + 0.000001, sol - nstd * perr)
    
            #~ fit_up_avg_end2end_fit_varvel_fct_time_fixvel = ballistic_msd_varvel_fct_time(xData, popt_up[0], avg_vel_clust_split_alltimes, popt_up[1], var_vel_clust_split_alltimes)  
            #~ fit_dw_avg_end2end_fit_varvel_fct_time_fixvel =ballistic_msd_varvel_fct_time(xData, popt_dw[0], avg_vel_clust_split_alltimes, max(0, popt_dw[1]), var_vel_clust_split_alltimes) 
            
            #popt_dw[1] = max(0, popt_dw[1])
            fit_up_avg_end2end_fit_varvel_fct_time_fixvel = func(xData, *popt_up)  
            fit_dw_avg_end2end_fit_varvel_fct_time_fixvel =func(xData, *popt_dw) 
             
            
            
            # goodness of fit
            
            #~ modelPredictions_avg_end2end_fit_varvel_fct_time_fixvel = ballistic_msd_varvel_fct_time(xData, pers_l_avg_end2end_fit_varvel_fct_time_fixvel, avg_vel_clust_split_alltimes, beta_avg_end2end_fit_varvel_fct_time_fixvel, var_vel_clust_split_alltimes) 
            modelPredictions_avg_end2end_fit_varvel_fct_time_fixvel = func(xData, *sol) 
            
            absError = modelPredictions_avg_end2end_fit_varvel_fct_time_fixvel  - yData
            
            SE = np.square(absError) # squared errors
            MSE = np.mean(SE) # mean squared errors
            RMSE  = np.sqrt(MSE) # Root Mean Squared Error, RMSE
            Rsquared  = 1.0 - (np.var(absError) / np.var(yData))
            print 'RMSE: ', RMSE
            print 'R-squared: ', Rsquared
            
            chisq_avg_end2end_fit_varvel_fct_time_fixvel = np.sum(SE/(ystd**2))/float(yData.size) # chisq
            print 'Chi-squared: ', chisq_avg_end2end_fit_varvel_fct_time_fixvel
               
            
            print 'compatible with rescaled fit space?'
            print pers_l_avg_end2end_fit_varvel_fixvel/avg_vel_clust_split_alltimes, beta_avg_end2end_fit_varvel_fixvel*avg_vel_clust_split_alltimes

        
        
            
        except ValueError:
            print 'Error in fitting'		
        except RuntimeError:
            print 'Exceeded max iterations in fitting'		
    


    print "theoretical"

    print mu, recog_width, mem_points, F0, avg_var_displ_x
    print VEL_Nteo, avg_vel_clust_split_alltimes
    print R_PERS_L, R_PERS_L/VEL_Nteo, PERS_TIME
    
    
    print "turns"
    print hist_nofull_time_turns.mean(), hist_nofull_time_turns_o_ttot_traj.mean()
    print hist_num_turns_o_tclade.mean(), hist_num_turns.mean(), hist_time_turns.mean(), hist_time_turns_o_ttot_traj.mean()
    
    
    # save observables
    
    
    dir_out_data=dir_in


    data = np.array([pers_l_avg_end2end_fit, err_pers_l_avg_end2end_fit, chisq_avg_end2end_fit, pers_l_avg_end2end_fit_varvel_fixvel, beta_avg_end2end_fit_varvel_fixvel, err_pers_l_avg_end2end_fit_varvel_fixvel, err_beta_avg_end2end_fit_varvel_fixvel, chisq_avg_end2end_fit_varvel_fixvel, diff_avg_angles_subseg_fit_rw_fct_time, err_diff_avg_angles_subseg_fit_rw_fct_time, chisq_avg_angles_subseg_fit_rw_fct_time, pers_l_avg_end2end_fit_constvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel, chisq_avg_end2end_fit_constvel_fct_time_fixvel, pers_l_avg_end2end_fit_varvel_fct_time_fixvel, beta_avg_end2end_fit_varvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel, err_beta_avg_end2end_fit_varvel_fct_time_fixvel, chisq_avg_end2end_fit_varvel_fct_time_fixvel])
    
    
    file_out='{inp}/global_features_pers_l_{real}.txt'.format(inp=dir_out_data, real=real)

    with open(file_out,'a') as f_handle:
        np.savetxt(f_handle, data, fmt='%15.15f', newline=" ")
        f_handle.write("\n")

    
    
    
    
    dir_out_plots='{inp}/realization_{real}/pers_l_fitting'.format(inp=dir_out_plots_tot,real=real) # directory with output plots
    
    
    if not os.path.exists(dir_out_plots):
        os.makedirs(dir_out_plots)	
        
        

    window=5
    
    print avg_end2end_coarse_gr_fct_time_list.size
    
    if avg_end2end_coarse_gr_fct_time_list.size >10:
        
        
    
        
        ## plot  avg_end2end_coarse_gr_list
        
        fig = plt.figure(figsize=thisfigsize)
        grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.22,
                     wspace=0.4, hspace=0.35)
        labeled_axes = []
        ax = plt.Subplot(fig, grid[0, 0])
        fig.add_subplot(ax)
        labeled_axes.append(ax)
        ax.errorbar(avg_time_diff_traj_list, avg_end2end_coarse_gr_fct_time_list, yerr=std_end2end_coarse_gr_fct_time_list, linestyle='-', color='g', elinewidth=0.1, label="data")

        ax.plot(avg_time_diff_traj_list, modelPredictions_avg_end2end_fit_constvel_fct_time_fixvel, linestyle='-', color='b', label="fit")
        
        ax.plot(avg_time_diff_traj_list, modelPredictions_avg_end2end_fit_varvel_fct_time_fixvel, linestyle='-', color='r', label="fit varvel")
        ax.fill_between(avg_time_diff_traj_list, fit_up_avg_end2end_fit_constvel_fct_time_fixvel, fit_dw_avg_end2end_fit_constvel_fct_time_fixvel, color='b', alpha=.25, label="2-sigma interval")
        
        ax.fill_between(avg_time_diff_traj_list, fit_up_avg_end2end_fit_varvel_fct_time_fixvel, fit_dw_avg_end2end_fit_varvel_fct_time_fixvel, color='r', alpha=.1, label="2-sigma interval varvel")
        
        ax.set_xlabel('segment length')
        ax.set_ylabel('R squared')
        ax.xaxis.labelpad = axis_labelpad
        ax.yaxis.labelpad = axis_labelpad
        ax.legend(frameon=False, ncol=1, columnspacing=0.5, handletextpad=0.2,
              loc='upper right', bbox_to_anchor=(2.15, 1.18))
        mpsetup.despine(ax) 
        
        
        #### finish figure ####
        labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                 xycoords=('axes fraction'), fontweight = 'bold')
        #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
        #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
        out_file='{out}/split_traj_avg_end2end_fct_time_fit_all.pdf'.format(out=dir_out_plots)
        #    print out_file
        fig.savefig(out_file)
        
        
          
    
        
        ## plot  avg_end2end_coarse_gr_list
        
        fig = plt.figure(figsize=thisfigsize)
        grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.22,
                     wspace=0.4, hspace=0.35)
        labeled_axes = []
        ax = plt.Subplot(fig, grid[0, 0])
        fig.add_subplot(ax)
        labeled_axes.append(ax)
        
            
        ax.errorbar(avg_time_diff_angles_traj_list, avg_angles_subseg_fct_time_list, yerr=std_angles_subseg_fct_time_list, linestyle='-', color='g', elinewidth=0.1, label="data")
        #errorbar(x, y0, yerr=noise, xerr=0, hold=True, ecolor='k', fmt='none', label='data')
    
        ax.fill_between(avg_time_diff_angles_traj_list, fit_up_avg_angles_subseg_fit_rw_fct_time, fit_dw_avg_angles_subseg_fit_rw_fct_time, color='b', alpha=.25, label="2-sigma interval")
                 
        ax.plot(avg_time_diff_angles_traj_list, modelPredictions_avg_angles_subseg_fit_rw_fct_time, linestyle='-', color='b', label="fit")
        
        
        
        ax.set_xlabel('segment length')
        ax.set_ylabel('R squared angle')
        ax.xaxis.labelpad = axis_labelpad
        ax.yaxis.labelpad = axis_labelpad
        ax.legend(frameon=False, ncol=1, columnspacing=0.5, handletextpad=0.2,
              loc='upper right', bbox_to_anchor=(2.15, 1.18))
        mpsetup.despine(ax) 
        
        
        #### finish figure ####
        labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                 xycoords=('axes fraction'), fontweight = 'bold')
        #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
        #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
        out_file='{out}/split_traj_avg_angle_msd_fct_time_fit_all.pdf'.format(out=dir_out_plots)
        #    print out_file
        fig.savefig(out_file)
        
        
    
        
        
    gc.collect()


