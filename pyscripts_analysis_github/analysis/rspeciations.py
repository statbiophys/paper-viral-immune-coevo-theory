
"""
TAKES THE OUTPUT OF CLUSTERING SCRIPT AND COUNTS THE NUMBER OF SPECIATIONS EVENT FOR DIFFERENT SPECIATIONS THRESHOLDS

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
    

def explore_tree_id(xs_all_clade, ys_all_clade, times_all_clade, ancestor_id, data_track_tot,  time_max = None):
    # recursively explore tree downstream, and add  id to a list and maximum time to other list 
    
    sons_clade = np.unique(data_track_tot[data_track_tot[:,8]==ancestor_id,2])
    
    print sons_clade.size
    if sons_clade.size >0: 
        
        for son in sons_clade:
            
            new_lineage = data_track_tot[data_track_tot[:,2] ==son,:]
            
            print son, new_lineage.shape
            times_all_clade.extend(new_lineage[:,0].tolist())
            xs_all_clade.extend(new_lineage[:,4].tolist())
            ys_all_clade.extend(new_lineage[:,5].tolist())
            
            time_max_clade= np.amax(new_lineage[:,0])
            
            if time_max is None or time_max > time_max_clade :
                
                ancestor_id = son
                explore_tree_id(xs_all_clade, ys_all_clade, times_all_clade, ancestor_id, data_track_tot,  time_max)
        
        
        

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
# ~ traj_anal=True

#if n_real >10:
#    n_real=10

print n_real   

thisfigsize = figsize
thisfigsize[1] *= 0.75

# get data

for real in np.arange(1,n_real+1):
    dir_in='{inp}/realization_{real}'.format(inp=dir_in_tot,real=real) # directory with input files
    dir_out_plots='{inp}/realization_{real}/rspec'.format(inp=dir_out_plots_tot,real=real) # directory with output plots
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
    
                
        # ~ print data
        print data.ndim
        
        if data.ndim > 1:    
        
    
            
            timelast=data[-1,0]
            
    if not os.path.exists(dir_out_plots):
        os.makedirs(dir_out_plots)	

    file_in='{inp}/space_clustering_fct_time_real_{real}.dat.npz'.format(inp=dir_in, real=real)
    
    if os.path.isfile(file_in):
     
        z = np.load(file_in)    
        data_vir_time = z['data']
        data_vir_time_allclust = z['data_allclust']
        data_vir_time_grad = z['data_allclust_grads']
        
         #~ np.savez_compressed(file_out,data = data, data_allclust = data_allclust, data_allclust_grads = data_allclust_grads, frac_1cl=frac_1cl, times_sharp_turn_split_alltimes_allclust= times_sharp_turn_split_alltimes_allclust)
       
        
        times_clust                                    = data_vir_time[:,0 ] #- data_vir_time[0,0 ] # THIS WAS WRONG IN CLUSTERING SCRIPT, RATES OVERESTIMATED
        n_clusters_split                         = data_vir_time[:,33]
      
      
    
    
    # data, output of clustering script
    outfile_track='{inp}/clusters_track_real_{real}.dat'.format(inp=dir_in, real=real)

    
        
    if os.path.isfile(outfile_track):
        
        data_track_tot	    = np.loadtxt(outfile_track)
        
        print data_track_tot.shape
        
        # data_son = np.array([time, state, ID_son, t_birth_son,  centr_son[0], centr_son[1], num_vir_son, size_son, ID_fath_son, centr_fath_son[0], centr_fath_son[1], num_vir_fath_son, size_fath_son])
                
        if data_track_tot.shape[0]>1:
    
    
            
            all_lineages= np.unique(data_track_tot[ :,2])
            
            time_tot_lineages=0
            # PRECOMPUTE THE TOTAL LINEAGES TIME
            for lin in all_lineages:
                
                lin_times=data_track_tot[data_track_tot[:,2]==lin, 0]
                
                lin_duration=lin_times[-1] - lin_times[0]
                
                time_tot_lineages+=lin_duration

            
            fathers_split= data_track_tot[ (data_track_tot[:,1]>1),8]
            data_split= data_track_tot[ (data_track_tot[:,1]>1),:]
            
            fathers_split_unique, fathers_split_unique_idx, fathers_split_unique_counts = np.unique(fathers_split, return_index=True, return_counts=True) 
            
            fathers_split_unique=fathers_split_unique[np.argsort(fathers_split_unique_idx)]
            fathers_split_unique_counts=fathers_split_unique_counts[np.argsort(fathers_split_unique_idx)]
            fathers_split_unique_idx=fathers_split_unique_idx[np.argsort(fathers_split_unique_idx)]
            
            
            
            # ~ print fathers_split
            print fathers_split.shape
            
            # ~ print fathers_split_unique
            # ~ print fathers_split_unique_counts
            # ~ print fathers_split_unique_idx
            print fathers_split_unique.shape
            
            
            # ~ idxs_1branch=fathers_split_unique_idx[ fathers_split_unique_counts==1 ]
            
            # ~ print idxs_1branch
            # ~ print idxs_1branch.shape
            
            thresholds = np.asarray([0.1 *recog_width, 0.5 *recog_width, 1 *recog_width, 10 *recog_width])
            
            
            
            num_split_tot       =   [[] for i in range(len(thresholds) + 1)]
            num_split_per_cl_tot=   [[] for i in range(len(thresholds) + 1)]
            
            
            
            #data_son = np.array([time, state, ID_son, t_birth_son,  centr_son[0], centr_son[1], num_vir_son, size_son, ID_fath_son, centr_fath_son[0], centr_fath_son[1], num_vir_fath_son, size_fath_son])
            for i, idx in enumerate(fathers_split_unique_idx): # CYCLE ON LINEAGE SPLITTING EVENTS COMPARING THE TWO RESULTING LINEAGES
                
                
                n_sons_active = data_split[idx,1]
                
                father = fathers_split_unique[i]
                
                ids_sons= np.unique(data_split[data_split[:,8]==father,2])
                
                # ~ print n_sons_active, ids_sons.shape
                if n_sons_active>1:  
                    
                    clades=[]
                    clades_sizes=[]
                    for id_new in ids_sons:
                        # ~ print "id ", id_new
                        track_id_orig=data_track_tot[data_track_tot[:,2]==id_new,:]
                        # ~ print track_id_orig.shape[0]
                        clades.append(track_id_orig)
                        clades_sizes.append(track_id_orig.shape[0])
                    
                    
                    clades_to_comp=[]
                    clades_sizes=np.asarray(clades_sizes)
                    if clades_sizes.size>2:
                        inds = np.argpartition(clades_sizes, 2)[-2:]
                        
                    else:    
                        inds=range(len(clades))
                        
                        
                    
                    for i in inds:
                        clades_to_comp.append(clades[i])
                        
                    if len(clades_to_comp) !=2:
                        print "error, too many clades ", len(clades_to_comp) 
                        sys.exit()
                        
                        
                    times_all= [ c[:,0] for c in clades_to_comp]
                    
                    # ~ print len(times_all), times_all[0].shape
                    
                    # ~ print times_all[0]
                    # ~ print times_all[1]
                    
                    times_common_idx  = np.nonzero(np.in1d(times_all[0], times_all[1]))[0]
                    times_common_idx2 = np.nonzero(np.in1d(times_all[1], times_all[0]))[0]
                    
                    
                    centr_clade1 =  clades_to_comp [0] [times_common_idx, 4:6]
                    centr_clade2 =  clades_to_comp [1] [times_common_idx2, 4:6]
                    
                    dists_clades= np.sqrt(np.sum((centr_clade1 - centr_clade2)**2., axis=1))
                    # ~ print dists_clades.shape 
                    
                    for i_thr, thr in enumerate(thresholds): # CYCLE ON THRESHOLDS
                        
                        if np.any(dists_clades>thr): # FIRST CHECK IF EVENT HAPPENS ON THE TWO CONSIDERED LINEAGES....
                        
                            idx_spec = np.argmax(dists_clades>thr)
                            
                            time_ev= times_all[0][times_common_idx[idx_spec]]
                            
                            num_split_tot       [0].append(time_ev)
                            num_split_per_cl_tot[0].append(time_ev)
                            
                            
                            for it2 in range(len(thresholds)):
                                toadd=0
                                if it2==i_thr:
                                    toadd=1 
                                
                                if np.all(times_clust != time_ev):
                                    print "Error, time not found"
                                    sys.exit()
                                
                                idx_time = np.nonzero(times_clust == time_ev)[0][0]
                                num_lin_now = n_clusters_split[idx_time] 
                                
                                num_split_tot       [it2+1].append(toadd)
                                num_split_per_cl_tot[it2+1].append(toadd)#/float(num_lin_now))
                
                            print "found event first couple"
                            # ~ print num_split_tot
                            # ~ print num_split_per_cl_tot
                            print thr
                
                            
                        else: # .... IF IT DOESN'T CYCLE ON ALL POSSIBLE DESCENDANT LINEAGES OF EACH ORIGINAL LINEAGE
                            
                            times_last_clades = [ c[-1,0] for c in clades_to_comp ] 
                            
                            ind_short_clade = times_last_clades.index(min(times_last_clades))
                            print "check descendants"
                            # ~ print ind_short_clade
                            
                            short_clade = clades_to_comp[ind_short_clade].copy() 
                            long_clade = clades_to_comp[ind_short_clade].copy() 
                            
                            # ~ print short_clade[-1, 1]
                            
                            if short_clade[-1, 1] > 0: 
                                
                                ancestor_id= short_clade[-1, 2]
                                
                                xs_all_clade_short    = []
                                ys_all_clade_short    = []
                                times_all_clade_short = []
                                
                                explore_tree_id(xs_all_clade_short, ys_all_clade_short, times_all_clade_short, ancestor_id, data_track_tot)
                                
                                # ~ print times_all_clade_short
                                # ~ print len(times_all_clade_short)
                                
                                if len(times_all_clade_short)>0:
                                    
                                    ancestor_id= long_clade[-1, 2]
                                    time_max_short = max(times_all_clade_short)
                                    
                                    xs_all_clade_long    = []
                                    ys_all_clade_long    = []
                                    times_all_clade_long = []
                                    
                                    explore_tree_id(xs_all_clade_long, ys_all_clade_long, times_all_clade_long, ancestor_id, data_track_tot, time_max_short)
                                    
                                    times_common_all= np.intersect1d(times_all_clade_long, times_all_clade_short)[::-1]
                                    
                                    for t_comp in times_common_all:
                                        
                                        # ~ print t_comp, times_common_all.size
                                        
                                        idx_time1 = np.nonzero(np.asarray(times_all_clade_long ) == t_comp)[0]
                                        idx_time2 = np.nonzero(np.asarray(times_all_clade_short) == t_comp)[0]
                                    
                                        for couple_idxs in itertools.product(idx_time1, idx_time2):
                                            i1=couple_idxs[0]
                                            i2=couple_idxs[1]
                                            dist= np.sqrt((xs_all_clade_long[i1] - xs_all_clade_short[i2] )**2. + (ys_all_clade_long[i1] - ys_all_clade_short[i2] )**2.)
                                            
                                            if dist > thr:
                                            
                                                time_ev= t_comp
                                                
                                                num_split_tot       [0].append(time_ev)
                                                num_split_per_cl_tot[0].append(time_ev)
                                                
                                                
                                                for it2 in range(len(thresholds)):
                                                    toadd=0
                                                    if it2==i_thr:
                                                        toadd=1 
                                                    
                                                    if np.all(times_clust != time_ev):
                                                        print "Error, time not found"
                                                        sys.exit()
                                                    
                                                    idx_time = np.nonzero(times_clust == time_ev)[0][0]
                                                    num_lin_now = n_clusters_split[idx_time] 
                                                    
                                                    num_split_tot       [it2+1].append(toadd)
                                                    num_split_per_cl_tot[it2+1].append(toadd)#/float(num_lin_now))
                                                  
                                                print "found event"
                                                print num_split_tot
                                                print num_split_per_cl_tot
                                                print thr
                                                break
                                                
                                        else:
                                            # Continue if the inner loop wasn't broken.
                                            continue
                                        # Inner loop was broken, break the outer.
                                        break
                        
                    #####################################################
                    
                    
                    
                    
                    
                    
                    
             
             
             
            num_split_tot        =  np.asarray(num_split_tot        ).T
            num_split_per_cl_tot =  np.asarray(num_split_per_cl_tot ).T
            
            
            
            print num_split_tot
            print num_split_tot.shape
            print num_split_per_cl_tot.shape
            
            num_split_tot=num_split_tot[np.argsort(num_split_tot[:,0]),:]
            num_split_per_cl_tot=num_split_per_cl_tot[np.argsort(num_split_per_cl_tot[:,0]),:]
            
            rate_num_split_tot                            = np.cumsum(num_split_tot[:,1:], axis=0)       [-1]/float( times_clust[-1] - data_vir_time[0,0 ])      
            # ~ rate_num_split_per_cl_tot                     = np.cumsum(num_split_per_cl_tot[:,1:], axis=0)[-1]/float( times_clust[-1] - data_vir_time[0,0 ]) 
            rate_num_split_per_cl_tot                     = np.cumsum(num_split_per_cl_tot[:,1:], axis=0)[-1]/float( time_tot_lineages) # RATE PER LINEAGE

            print num_split_tot[:,1:]
            print rate_num_split_tot
            print rate_num_split_tot.shape
            
            
            print time_tot_lineages
            print times_clust[-1] - data_vir_time[0,0 ]
            
            
            

            # ~ data = np.array([frac_1cl, avg_n_clusters_dbscan_opt_kvar_list, avg_size_clusters_split, avg_n_clusters_split, avg_num_vir_clusters_split, avg_num_coords_clusters_split, avg_num_ext_tot, avg_num_split_tot, avg_num_ext_per_cl_tot, avg_num_split_per_cl_tot, avg_num_vir_o_nclust_tot_list, avg_vel_clust_split_alltimes_allclust, avg_var_parall_clust_split_alltimes_allclust, avg_var_perp_clust_split_alltimes_allclust, avg_var_tot_clust_split_alltimes_allclust, avg_inter_clusters_dist_split, max_inter_clusters_dist_split, max_inter_clusters_dist_max_split, var_vel_clust_split_alltimes_allclust, var_var_parall_clust_split_alltimes_allclust, var_var_perp_clust_split_alltimes_allclust, var_var_tot_clust_split_alltimes_allclust, rate_num_ext_tot, rate_num_split_tot, rate_num_ext_per_cl_tot, rate_num_split_per_cl_tot, num_sharp_turns, rate_sharp_turns, avg_shortmem_err_ys, avg_shortmem_err_svsx, frac_shortmem_err_ys, frac_shortmem_err_svsx,    diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, diff_ang_fit_rw_fct_time, err_diff_ang_fit_rw_fct_time, chisq_ang_fit_rw_fct_time, diff_ang_short_fit_rw_fct_time, err_diff_ang_short_fit_rw_fct_time, chisq_ang_short_fit_rw_fct_time ])
            data = np.concatenate([rate_num_split_tot,rate_num_split_per_cl_tot])
            
            
            file_out='{inp}/global_features_rspeciation_{real}.txt'.format(inp=dir_in, real=real)
        
            with open(file_out,'w') as f_handle:
                np.savetxt(f_handle, data, fmt='%15.15f', newline=" ")
                f_handle.write("\n")
                
                
                
            
            ## plot  num_split_tot num_ext_tot cum as fct of time
            
            fig = plt.figure(figsize=thisfigsize)
            grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.22,
                         wspace=0.4, hspace=0.35)
            labeled_axes = []
            ax = plt.Subplot(fig, grid[0, 0])
            fig.add_subplot(ax)
            labeled_axes.append(ax)
            ax.plot(num_split_per_cl_tot[:,0], np.cumsum(num_split_per_cl_tot[:,1]), linestyle='-', color='b', label='0d1')
            ax.plot(num_split_per_cl_tot[:,0], np.cumsum(num_split_per_cl_tot[:,2]), linestyle='-', color='g', label='0d5')
            ax.plot(num_split_per_cl_tot[:,0], np.cumsum(num_split_per_cl_tot[:,3]), linestyle='-', color='r', label='1')
            ax.plot(num_split_per_cl_tot[:,0], np.cumsum(num_split_per_cl_tot[:,4]), linestyle='-', color='m', label='10')
            ax.set_xlabel('time (y)')
            ax.set_ylabel('cumulative number of events')
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
            out_file='{out}/split_clusters_ext_birth_per_cl_cumulative_{real}.png'.format(out=dir_out_plots, real=real)
            #    print out_file
            fig.savefig(out_file)
            
                       
            
        
