
"""
COLLECTS FINAL SUMMARY STATISTICS TO PLOT IN THE MANUSCRIPT

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
import scipy.integrate as integrate
import sys
sys.path.append('..')
from lib.mppaper import *
import lib.mpsetup as mpsetup
import shutil




    


dir_io=sys.argv[1] # directory with input files
dir_in_tot='{inp}/realizations'.format(inp=dir_io) # directory with input files
#dir_in = dir_in.replace(".", "d")

dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots
#dir_out_plots = dir_out_plots.replace(".", "d")
  
      
    
    
param_file='{inp}/parameters_backup.dat'.format(inp=dir_io)










if os.path.isfile(param_file):
    
    
    #    param_file<< "#  1 <mu> 	2 <recognition width>	3 <sigma>	4 <eps>		5 <I>		6 <people number>	7 <coverage save rate>		8 <fraction infected>	9 <number infections>	10 <full_frame save time>	11 <n_reals>;
    params = np.genfromtxt(param_file, dtype="f8,f8,f8,i8,i8,i8,i8, |S10, |S10, |S10, i8, i8, f8", names=['mu','rec_width','jump_size', 'ppl_num', 'maxinfections','save_full_time','n_real','initial_condition','fake_initial_condition','phylo_subsample','mem_points','t_init','F0'])
    in_cond=params['initial_condition']
    
    ppl_num=params['ppl_num']
    F0=params['F0']
    mem_points=params['mem_points']
    rec_width=params['rec_width']
    mu=params['mu']
    #~ jump_size=params['jump_size']
    
    
        
    #~ mem_points                                    =-1
    #~ F0                                            =-1
    #~ ppl_num                                       =-1
    #~ mu                                            =-1
    #~ rec_width                                     =-1
    
    
    
    
    time_tot_run                                   = -1
    probext                                        = -1
    probexpl                                       = -1
    numext                                        = -1
    numexpl                                       = -1
    avg_I                                        = 0
    std_I                                        = 0
    var_log_I                                    = 0
    avg_fit_avg                                  = 0
    var_fit_avg                                  = 0
    avg_fit_var                                  = 0
    var_fit_var                                  = 0
    #~ avg_fit_nose                              = 0
    #~ var_fit_nose                              = 0
    #~ avg_vel                                   = 0
    #~ var_vel                                   = 0
    avg_num_IS_coords                            = 0
    avg_num_vir_coords                           = 0
    avg_mean_displ_x                             = 0
    avg_mean_displ_y                             = 0
    avg_var_displ_x                              = 0
    avg_var_displ_y                              = 0
    avg_mean_displ_tot                           = 0
    diff_maxfit_fit_rw_fct_time                  = 0
    err_diff_maxfit_fit_rw_fct_time              = 0
    chisq_maxfit_fit_rw_fct_time                 = 0
    avg_sig_est                                  = 0
    var_sig_est                                  = 0
    avg_time_switch                              = 0
    avg_time_fullconv                            = 0
    max_fitn_abserr_switch_max_rel               = 0
    avg_fitn_abserr_switch_max_rel               = 0
    #~ avg_conv_abserr_switch_max_rel            = 0
    #~ avg_conv_err_switch_max_rel               = 0
    frac_1cl                                     = 0
    #~ avg_n_clusters_dbscan_opt_kvar_list       = 0
    avg_size_clusters_split                      = 0
    avg_n_clusters_split                         = 0
    avg_num_vir_clusters_split                   = 0
    avg_num_coords_clusters_split                = 0
    avg_num_ext_tot                              = 0
    avg_num_split_tot                            = 0
    avg_num_ext_per_cl_tot                       = 0
    avg_num_split_per_cl_tot                     = 0
    avg_num_vir_o_nclust_tot_list                = 0
    avg_vel_clust_split_alltimes_allclust        = 0
    avg_var_parall_clust_split_alltimes_allclust = 0
    avg_var_perp_clust_split_alltimes_allclust   = 0
    avg_var_tot_clust_split_alltimes_allclust    = 0
    avg_inter_clusters_dist_split                = 0
    max_inter_clusters_dist_split                = 0
    max_inter_clusters_dist_max_split            = 0
    var_vel_clust_split_alltimes_allclust        = 0
    var_var_parall_clust_split_alltimes_allclust = 0
    var_var_perp_clust_split_alltimes_allclust   = 0
    var_var_tot_clust_split_alltimes_allclust    = 0
    rate_num_ext_tot                             = 0
    rate_num_split_tot                           = 0
    rate_num_ext_per_cl_tot                      = 0
    rate_num_split_per_cl_tot                    = 0
    num_sharp_turns                              = 0
    rate_sharp_turns                             = 0
    avg_shortmem_err_ys                          = 0
    avg_shortmem_err_svsx                        = 0
    frac_shortmem_err_ys                         = 0
    frac_shortmem_err_svsx                       = 0
    diff_maxfit_fit_rw_fct_time                  = 0
    err_diff_maxfit_fit_rw_fct_time              = 0
    chisq_maxfit_fit_rw_fct_time                 = 0
    #~ diff_ang_fit_rw_fct_time                  = 0
    #~ err_diff_ang_fit_rw_fct_time              = 0
    #~ chisq_ang_fit_rw_fct_time                 = 0
    #~ diff_ang_short_fit_rw_fct_time            = 0
    #~ err_diff_ang_short_fit_rw_fct_time        = 0
    #~ chisq_ang_short_fit_rw_fct_time           = 0
    timelast                                     = -1
    NUM_VIR                                      = -1
    TAU                                          = -1
    VEL                                          = -1
    S                                            = -1
    SIG_SPACE                                    = -1
    NOSE_SPACE                                   = -1
    R_PERS_L                                     = -1
    
    avg_grad_fittest_split_alltimes_allclust = 0    
    avg_shortmem_gradx_bulk_split_alltimes_allclust = 0 

    avg_S_est                                  = 0

    pers_l_avg_end2end_fit                                = -1
    err_pers_l_avg_end2end_fit                            = -1
    chisq_avg_end2end_fit                                 = -1
    pers_l_avg_end2end_fit_varvel_fixvel                  = -1
    beta_avg_end2end_fit_varvel_fixvel                    = -1
    err_pers_l_avg_end2end_fit_varvel_fixvel              = -1
    err_beta_avg_end2end_fit_varvel_fixvel                = -1
    chisq_avg_end2end_fit_varvel_fixvel                   = -1
    diff_avg_angles_subseg_fit_rw_fct_time                = -1
    err_diff_avg_angles_subseg_fit_rw_fct_time            = -1
    chisq_avg_angles_subseg_fit_rw_fct_time               = -1
    pers_l_avg_end2end_fit_constvel_fct_time_fixvel       = -1
    err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel   = -1
    chisq_avg_end2end_fit_constvel_fct_time_fixvel        = -1
    pers_l_avg_end2end_fit_varvel_fct_time_fixvel         = -1
    beta_avg_end2end_fit_varvel_fct_time_fixvel           = -1
    err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel     = -1
    err_beta_avg_end2end_fit_varvel_fct_time_fixvel       = -1
    chisq_avg_end2end_fit_varvel_fct_time_fixvel          = -1
    
    
    
    rate_num_split_tot_0d1            =0
    rate_num_split_tot_0d5            =0
    rate_num_split_tot_1              =0
    rate_num_split_tot_10             =0
    rate_num_split_per_cl_tot_0d1     =0
    rate_num_split_per_cl_tot_0d5     =0
    rate_num_split_per_cl_tot_1       =0
    rate_num_split_per_cl_tot_10      =0
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    print params
    print params.shape
    n_real =params['n_real']
    print n_real
    
    
    thisfigsize = figsize
    thisfigsize[1] *= 0.75
      
    
    
    real=1
    
    dir_in='{inp}/realization_{real}'.format(inp=dir_in_tot,real=real) # directory with input files
    
    file_in='{inp}/IC_trav_wave.txt'.format(inp=dir_in, real=real)
     
    if os.path.isfile(file_in):
     
        data_travwave = np.loadtxt(file_in)
        
        NUM_VIR      = data_travwave[0]  
        TAU          = data_travwave[1]
        VEL          = data_travwave[2]
        S            = data_travwave[3]
        SIG_SPACE    = data_travwave[4]
        NOSE_SPACE   = data_travwave[5]
        R_PERS_L     = data_travwave[6]
        
    
    file_in_avg_npz_compr='{inp}/evo_mean_stats_real_{real}.dat'.format(inp=dir_in, real=real)
    
    if os.path.isfile(file_in_avg_npz_compr):
        
        #data = np.loadtxt(file_in_avg_npz_compr)
        
        #with open(file_in_avg_npz_compr, encoding='utf-8', errors='ignore') as f:
        #    #content = f.read().splitlines()
        #    data = np.loadtxt(f)
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
        
    
            time = data[:,0]
            
            time_tot_run=time[-1]
            timelast=time[-1]
            
            time_mod=time[time>10000] # throw first 100 years
            
            
            
            if time_mod.size>0:
                sec_thr=min(100000, time_mod[time_mod.size/4])
            
                time_mod=time_mod[time_mod>sec_thr]
            
                #~ print sec_thr, time.size
            
        
        
            time_ss = time_mod[:]
            
            
            print time_ss.size
            
            if time_ss.size > 3000: # at least 30000 cycles after transient, otherwise parameter does not exist
                
                time_fin=time_ss[-1]
                
                time_tot_run=time_fin
                           
                data= data[time>sec_thr,:]
 
                
                num_x       = data[:,1]
                num_y       = data[:,2]
                area       = data[:,3]
                num_IS_coords       = data[:,4]
                num_vir_coords       = data[:,5]
                num_IS_tot       = data[:,6]
                num_vir_tot       = data[:,7]
                avg_fitn       = data[:,8]
                var_fitn       = data[:,9]
                start_x       = data[:,10]
                start_y       = data[:,11]
                vir_x          = data[:,12] # 
                vir_y          = data[:,13] # 
                fit_nose          = data[:,14]- avg_fitn # 
                fit_nose_bare          = data[:,14]  # 
                
                num_IS_coords_upd       = data[:,15]
                num_IS_tot_upd       = data[:,16]
                
                #"25 mean_displ_x_ "   <<setw(30)<<"26 mean_displ_y_ "   <<setw(30)<<"27 mean_displ_tot_ "   <<setw(30)<<"28 var_displ_x_ "   <<setw(30)<<"29 var_displ_y_ "   <<setw(30)<<"30 count_displ_ "   <<setw(30)<<"31 x1_max_fitn "   <<setw(30)<<"32 y1_max_fitn "
                
                mean_displ_x     = data[:,25]
                mean_displ_y     = data[:,26]
                mean_displ_tot   = data[:,27]
                var_displ_x      = data[:,28]
                var_displ_y      = data[:,29]
                count_displ      = data[:,30]
                x_max_fitn       = data[:,31]
                y_max_fitn       = data[:,32]
                
                sig_est=np.sqrt((x_max_fitn - vir_x)**2 + (y_max_fitn - vir_y)**2)
                    
                S_est=fit_nose[sig_est>0]/sig_est[sig_est>0]
                   
                avg_S_est  = S_est.mean()
           
                file_in='{inp}/cumulative_events_rerun_{real}.txt'.format(inp=dir_io, real=real)
                 
                if os.path.isfile(file_in):
                 
                    data_ext = np.loadtxt(file_in)
                    print data_ext
                    
                    if data_ext.ndim > 1  and data_ext.size > 1: 
                        data_ext = data_ext[data_ext[:,2]>0,:]
                        
                        data_ext_trans = data_ext[data_ext[:,0]<=sec_thr,:]
                        data_ext       = data_ext[data_ext[:,0]>sec_thr,:]
                    elif data_ext.size > 1: 
                        data_ext = data_ext[data_ext[2]>0]
                        
                        data_ext_trans = data_ext[data_ext[0]<=sec_thr]
                        data_ext       = data_ext[data_ext[0]>sec_thr]
            
                    if data_ext_trans.ndim > 1 and data_ext_trans.size > 0:  # last rerun in transient 
                        data_ext_trans=data_ext_trans[-1]
                        
                    numext   = 0
                    numexpl  = 0
                    probext  = 0
                    probexpl = 0
                          
                    if data_ext.ndim > 1 and data_ext.size > 1:    
                        tprec    = data_ext[:,0]                    
                        trun     = data_ext[:,2]                    
                        numext   = data_ext[:,3]                    
                        numexpl  = data_ext[:,4]                    
                        probext  = data_ext[:,5]                    
                        probexpl = data_ext[:,6]
                    elif data_ext.size > 1:
                        tprec    = data_ext[0]                    
                        trun     = data_ext[2]                    
                        numext   = data_ext[3]                    
                        numexpl  = data_ext[4]                    
                        probext  = data_ext[5]                    
                        probexpl = data_ext[6]
                        
                    print "data ext"
                    print data_ext
                    print data_ext.size
                    print data_ext.ndim
                        
                    print probext
    
                    if data_ext.size > 1:
                        
                        tfirst        =0
                        numextfirst   =0
                        numexplfirst  =0
                        
                        if data_ext_trans.size > 0:  # last rerun in transient 
                            tfirst        = data_ext_trans[2] 
                            numextfirst   = data_ext_trans[3] 
                            numexplfirst  = data_ext_trans[4] 
                             
                        
                        
                        if  np.ndim(trun) == 0:
                            time_tot_run = trun 
                            probext      = probext 
                            probexpl     = probexpl
                            numext   = numext  
                            numexpl  = numexpl                             
                        else:
                            time_tot_run = trun[-1] 
                            probext      = probext [-1] 
                            probexpl     = probexpl[-1] 
                            numext       = numext  [-1]
                            numexpl      = numexpl [-1]
                    
                        time_tot_run -= tfirst      
                        numext       -= numextfirst 
                        numexpl      -= numexplfirst
                    
                        
                        if probext!=1 and probexpl!=1: 
                            time_tot_run+= time_fin - tprec[-1]
                        
                    print probext
    
            
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
                    std_I                            = data_glob[1]
                    var_log_I                        = data_glob[2]
                    avg_fit_avg                      = data_glob[3]
                    var_fit_avg                      = data_glob[4]
                    avg_fit_var                      = data_glob[5]
                    var_fit_var                      = data_glob[6]
                    #~ avg_fit_nose                  = data_glob[7]
                    #~ var_fit_nose                  = data_glob[8]
                    #~ avg_vel                       = data_glob[9]
                    #~ var_vel                       = data_glob[10]
                    avg_num_IS_coords                = data_glob[11]
                    avg_num_vir_coords               = data_glob[12]
                    avg_mean_displ_x                 = data_glob[13]
                    avg_mean_displ_y                 = data_glob[14]
                    avg_var_displ_x                  = data_glob[15]
                    avg_var_displ_y                  = data_glob[16]
                    avg_mean_displ_tot               = data_glob[17]
                    diff_maxfit_fit_rw_fct_time      = data_glob[18]
                    err_diff_maxfit_fit_rw_fct_time  = data_glob[19]
                    chisq_maxfit_fit_rw_fct_time     = data_glob[20]
                    avg_sig_est                      = data_glob[21]
                    var_sig_est                      = data_glob[22]
                
    
                
                
                
    
            
                file_in='{inp}/global_features_bench_fitupd_{real}.txt'.format(inp=dir_in, real=real)
                 
                if os.path.isfile(file_in):
                 
                    data_glob = np.loadtxt(file_in)                
        
                    if data_glob.ndim >1:
                        data_glob=data_glob[-1]

                     #~ data = np.array([avg_time_switch, avg_time_fullconv, max_fitn_abserr_switch_max_rel, avg_fitn_abserr_switch_max_rel, avg_conv_abserr_switch_max_rel, avg_conv_err_switch_max_rel])
                                            
                    avg_time_switch                = data_glob[0]
                    avg_time_fullconv              = data_glob[1]
                    max_fitn_abserr_switch_max_rel = data_glob[2]
                    avg_fitn_abserr_switch_max_rel = data_glob[3]
                    #~ avg_conv_abserr_switch_max_rel = data_glob[4]
                    #~ avg_conv_err_switch_max_rel    = data_glob[5]                  
        
    
    
    
    
    
    
    
    
    
    
                file_in='{inp}/global_features_virclust_{real}.txt'.format(inp=dir_in, real=real)
                 
                if os.path.isfile(file_in):
                 
                    data_glob = np.loadtxt(file_in)
    
                    if data_glob.ndim >1:
                        data_glob=data_glob[-1]

            
                    
                    #~ data = np.array([times, n_clusters_dbscan_list, n_clusters_dbscan_est_list, n_clusters_dbscan_opt_kvar_list, eps_dbscan_est_list, eps_dbscan_opt_kvar_list,  size_clusters_dbscan_list, size_clusters_dbscan_est_list, size_clusters_dbscan_opt_kvar_list, density_clusters_dbscan_list, density_clusters_dbscan_est_list, density_clusters_dbscan_opt_kvar_list, n_clusters_dbscan_opt_dispersion_list, n_clusters_dbscan_opt_dispersion_CV_list, n_clusters_dbscan_opt_kCV_list, mean_dispersion_dbscan_list, mean_dispersion_dbscan_est_list, mean_dispersion_dbscan_opt_list, mean_dispersion_CV_dbscan_list, mean_dispersion_CV_dbscan_est_list, mean_dispersion_CV_dbscan_opt_list, mean_kCV_dbscan_list, mean_kCV_dbscan_est_list, mean_kCV_dbscan_opt_list, mean_kvar_dbscan_list, mean_kvar_dbscan_est_list, mean_kvar_dbscan_opt_list, inter_clusters_dist_dbscan_opt_kvar_list, inter_clusters_dist_max_dbscan_opt_kvar_list, size_clusters_split, inter_clusters_dist_split, inter_clusters_dist_max_split, density_clusters_split, n_clusters_split, num_vir_clusters_split, num_coords_clusters_split, num_ext_tot, num_split_tot, num_ext_per_cl_tot, num_split_per_cl_tot, num_vir_o_nclust_tot_list, vel_clust_split_alltimes, var_parall_clust_split_alltimes, var_perp_clust_split_alltimes, var_tot_clust_split_alltimes])
                    
                    #~ data_allclust = np.array([vel_clust_split_alltimes_allclust, var_parall_clust_split_alltimes_allclust, var_perp_clust_split_alltimes_allclust, var_tot_clust_split_alltimes_allclust])
                    
                    #~ data_allclust_grads = np.array([times_grad_split_alltimes_allclust, grad_fittest_split_alltimes_allclust, gradx_bulk_split_alltimes_allclust, grady_fittest_split_alltimes_allclust, grady_bulk_split_alltimes_allclust, shortmem_gradx_bulk_split_alltimes_allclust, shortmem_grady_fittest_split_alltimes_allclust, shortmem_grady_bulk_split_alltimes_allclust])
                    
                    
                    #~ data_allclust_diffy = np.array([times_diff_fity_split_alltimes_allclust, diff_fity_split_alltimes_allclust, diff_fitx_split_alltimes_allclust, ang_past_dir_shortmem_split_alltimes_allclust, ang_past_dir_split_alltimes_allclust])
                    
                    
                        
                        
                        
                    #~ file_out='{inp}/space_clustering_fct_time_real_{real}.dat'.format(inp=dir_in, real=real)
                    
                    
            
            
            
            
    
                    #~ data = np.array([frac_1cl, avg_n_clusters_dbscan_opt_kvar_list, avg_size_clusters_split, avg_n_clusters_split, avg_num_vir_clusters_split, avg_num_coords_clusters_split, avg_num_ext_tot, avg_num_split_tot, avg_num_ext_per_cl_tot, avg_num_split_per_cl_tot, avg_num_vir_o_nclust_tot_list, avg_vel_clust_split_alltimes_allclust, avg_var_parall_clust_split_alltimes_allclust, avg_var_perp_clust_split_alltimes_allclust, avg_var_tot_clust_split_alltimes_allclust, avg_inter_clusters_dist_split, max_inter_clusters_dist_split, max_inter_clusters_dist_max_split, var_vel_clust_split_alltimes_allclust, var_var_parall_clust_split_alltimes_allclust, var_var_perp_clust_split_alltimes_allclust, var_var_tot_clust_split_alltimes_allclust, rate_num_ext_tot, rate_num_split_tot, rate_num_ext_per_cl_tot, rate_num_split_per_cl_tot, num_sharp_turns, rate_sharp_turns, avg_shortmem_err_ys, avg_shortmem_err_svsx, frac_shortmem_err_ys, frac_shortmem_err_svsx,    diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, diff_ang_fit_rw_fct_time, err_diff_ang_fit_rw_fct_time, chisq_ang_fit_rw_fct_time, diff_ang_short_fit_rw_fct_time, err_diff_ang_short_fit_rw_fct_time, chisq_ang_short_fit_rw_fct_time ])
    
    
                    frac_1cl                                     = data_glob[0]
                    #~ avg_n_clusters_dbscan_opt_kvar_list       = data_glob[1]
                    avg_size_clusters_split                      = data_glob[2]
                    avg_n_clusters_split                         = data_glob[3]
                    avg_num_vir_clusters_split                   = data_glob[4]
                    avg_num_coords_clusters_split                = data_glob[5]
                    avg_num_ext_tot                              = data_glob[6]
                    avg_num_split_tot                            = data_glob[7]
                    avg_num_ext_per_cl_tot                       = data_glob[8]
                    avg_num_split_per_cl_tot                     = data_glob[9]
                    avg_num_vir_o_nclust_tot_list                = data_glob[10]
                    avg_vel_clust_split_alltimes_allclust        = data_glob[11]
                    avg_var_parall_clust_split_alltimes_allclust = data_glob[12]
                    avg_var_perp_clust_split_alltimes_allclust   = data_glob[13]
                    avg_var_tot_clust_split_alltimes_allclust    = data_glob[14]
                    avg_inter_clusters_dist_split                = data_glob[15]
                    max_inter_clusters_dist_split                = data_glob[16]
                    max_inter_clusters_dist_max_split            = data_glob[17]
                    var_vel_clust_split_alltimes_allclust        = data_glob[18]
                    var_var_parall_clust_split_alltimes_allclust = data_glob[19]
                    var_var_perp_clust_split_alltimes_allclust   = data_glob[20]
                    var_var_tot_clust_split_alltimes_allclust    = data_glob[21]
                    rate_num_ext_tot                             = data_glob[22]
                    rate_num_split_tot                           = data_glob[23]
                    rate_num_ext_per_cl_tot                      = data_glob[24]
                    rate_num_split_per_cl_tot                    = data_glob[25]
                    num_sharp_turns                              = data_glob[26]
                    rate_sharp_turns                             = data_glob[27]
                    avg_shortmem_err_ys                          = data_glob[28]
                    avg_shortmem_err_svsx                        = data_glob[29]
                    frac_shortmem_err_ys                         = data_glob[30]
                    frac_shortmem_err_svsx                       = data_glob[31]
                    diff_maxfit_fit_rw_fct_time                  = data_glob[32]
                    err_diff_maxfit_fit_rw_fct_time              = data_glob[33]
                    chisq_maxfit_fit_rw_fct_time                 = data_glob[34]
                    #~ diff_ang_fit_rw_fct_time                  = data_glob[35]
                    #~ err_diff_ang_fit_rw_fct_time              = data_glob[36]
                    #~ chisq_ang_fit_rw_fct_time                 = data_glob[37]
                    #~ diff_ang_short_fit_rw_fct_time            = data_glob[38]
                    #~ err_diff_ang_short_fit_rw_fct_time        = data_glob[39]
                    #~ chisq_ang_short_fit_rw_fct_time           = data_glob[40]
            
    
                print avg_vel_clust_split_alltimes_allclust
        
    

    
    
                file_in='{inp}/space_clustering_fct_time_real_{real}.dat.npz'.format(inp=dir_in, real=real)
                 
                if os.path.isfile(file_in):
                 
                    z = np.load(file_in)    
                    data_vir_time = z['data']
                    data_vir_time_allclust = z['data_allclust']
                    data_vir_time_grad = z['data_allclust_grads']
                    
                     #~ np.savez_compressed(file_out,data = data, data_allclust = data_allclust, data_allclust_grads = data_allclust_grads, frac_1cl=frac_1cl, times_sharp_turn_split_alltimes_allclust= times_sharp_turn_split_alltimes_allclust)
                   
                    
                    times                                    = data_vir_time[:,0 ] - data_vir_time[0,0 ] # 
                    n_clusters_split                         = data_vir_time[:,33]
                    num_ext_tot                              = data_vir_time[:,36]
                    num_split_tot                            = data_vir_time[:,37]
                    num_ext_per_cl_tot                       = data_vir_time[:,38]
                    num_split_per_cl_tot                     = data_vir_time[:,39]
                    inter_clusters_dist_split                = data_vir_time[:,30]
                    inter_clusters_dist_max_split            = data_vir_time[:,31]
                    
                    size_clusters_split                      = data_vir_time[:,29]
                    num_vir_clusters_split                   = data_vir_time[:,34]
                    num_coords_clusters_split                = data_vir_time[:,35]
                    num_vir_o_nclust_tot_list                = data_vir_time[:,40]
                     
                     
                     
                    vel_clust_split_alltimes_allclust        =  data_vir_time_allclust[0,:]
                    var_parall_clust_split_alltimes_allclust =  data_vir_time_allclust[1,:]
                    var_perp_clust_split_alltimes_allclust   =  data_vir_time_allclust[2,:]
                    var_tot_clust_split_alltimes_allclust    =  data_vir_time_allclust[3,:]
                    
                    grad_fittest_split_alltimes_allclust    = data_vir_time_grad[1,:]
                    shortmem_gradx_bulk_split_alltimes_allclust    = data_vir_time_grad[5,:]
                    
                    
                    #~ print times
                    #~ print grad_fittest_split_alltimes_allclust
                    #~ print times
                    print data_vir_time.shape
                    print data_vir_time_allclust.shape
                    print data_vir_time_allclust[0,:].mean()
                    print data_vir_time_grad.shape
                    
                    
                    
                    thrown_times = times
        
                    file_in='{inp}/thrown_space_clustering_fct_time_real_{real}.dat'.format(inp=dir_io, real=real)
                     
                    if os.path.isfile(file_in):
                        data_vir_time_thrown = np.loadtxt(file_in)
                        
                        thrown_n_clusters_split                = data_vir_time[:,7]
                        thrown_num_ext_tot                     = data_vir_time[:,10]
                        thrown_num_split_tot                   = data_vir_time[:,11]
                        thrown_num_ext_per_cl_tot              = data_vir_time[:,12]
                        thrown_num_split_per_cl_tot            = data_vir_time[:,13]
                        thrown_inter_clusters_dist_split       = data_vir_time[:,5]
                        thrown_inter_clusters_dist_max_split   = data_vir_time[:,6]
                        
            
                        #~ data = np.array([times, n_clusters_dbscan_opt_kvar_list, inter_clusters_dist_dbscan_opt_kvar_list, inter_clusters_dist_max_dbscan_opt_kvar_list, size_clusters_split, inter_clusters_dist_split, inter_clusters_dist_max_split, n_clusters_split, num_vir_clusters_split, num_coords_clusters_split, num_ext_tot, num_split_tot, num_ext_per_cl_tot, num_split_per_cl_tot, num_vir_o_nclust_tot_list, vel_clust_split_alltimes, var_parall_clust_split_alltimes, var_perp_clust_split_alltimes, var_tot_clust_split_alltimes])
            
                
                        n_clusters_split                = np.append(n_clusters_split             , thrown_n_clusters_split              )
                        num_ext_tot                     = np.append(num_ext_tot                  , thrown_num_ext_tot                   )
                        num_split_tot                   = np.append(num_split_tot                , thrown_num_split_tot                 )
                        num_ext_per_cl_tot              = np.append(num_ext_per_cl_tot           , thrown_num_ext_per_cl_tot            )
                        num_split_per_cl_tot            = np.append(num_split_per_cl_tot         , thrown_num_split_per_cl_tot          )
                        inter_clusters_dist_split       = np.append(inter_clusters_dist_split    , thrown_inter_clusters_dist_split     )
                        inter_clusters_dist_max_split   = np.append(inter_clusters_dist_max_split, thrown_inter_clusters_dist_max_split )
                                
                        thrown_times = np.arange(n_clusters_split.shape[0])*100
                        
                        print thrown_times[-1]
                
                

            
                    #~ if avg_n_clusters_split==0: # clustering analysis broke before gthering glob

                    avg_size_clusters_split                       = size_clusters_split                       .mean()  
                    avg_num_vir_clusters_split                    = num_vir_clusters_split                    .mean()  
                    avg_num_coords_clusters_split                 = num_coords_clusters_split                 .mean()  
                    avg_num_vir_o_nclust_tot_list                 = num_vir_o_nclust_tot_list                 .mean()  
                    avg_vel_clust_split_alltimes_allclust         = vel_clust_split_alltimes_allclust         .mean()           
                    avg_var_parall_clust_split_alltimes_allclust  = var_parall_clust_split_alltimes_allclust  .mean()           
                    avg_var_perp_clust_split_alltimes_allclust    = var_perp_clust_split_alltimes_allclust    .mean()           
                    avg_var_tot_clust_split_alltimes_allclust     = var_tot_clust_split_alltimes_allclust     .mean()           
                    var_vel_clust_split_alltimes_allclust         = vel_clust_split_alltimes_allclust         .var()   
                    var_var_parall_clust_split_alltimes_allclust  = var_parall_clust_split_alltimes_allclust  .var()   
                    var_var_perp_clust_split_alltimes_allclust    = var_perp_clust_split_alltimes_allclust    .var()   
                    var_var_tot_clust_split_alltimes_allclust     = var_tot_clust_split_alltimes_allclust     .var()  
                    
                    print vel_clust_split_alltimes_allclust
                    print vel_clust_split_alltimes_allclust.shape
                    print avg_vel_clust_split_alltimes_allclust
                                        
                    #~ if os.path.isfile(file_in) or avg_n_clusters_split==0:
                    
                    frac_1cl=(n_clusters_split[:]==1).sum()/float(n_clusters_split[:].shape[0])
        
                    avg_n_clusters_split              = n_clusters_split                          .mean()  
                    avg_num_ext_tot                   = num_ext_tot                               .mean()  
                    avg_num_split_tot                 = num_split_tot                             .mean()  
                    avg_num_ext_per_cl_tot            = num_ext_per_cl_tot                        .mean()  
                    avg_num_split_per_cl_tot          = num_split_per_cl_tot                      .mean()            
                    avg_inter_clusters_dist_split     = inter_clusters_dist_split                 .mean()
                    max_inter_clusters_dist_split     = inter_clusters_dist_split                  .max()   
                    max_inter_clusters_dist_max_split = inter_clusters_dist_max_split              .max()
                    rate_num_ext_tot                  = np.cumsum(num_ext_tot)         [-1]/float( thrown_times[-1])    
                    rate_num_split_tot                = np.cumsum(num_split_tot)       [-1]/float( thrown_times[-1])      
                    rate_num_ext_per_cl_tot           = np.cumsum(num_ext_per_cl_tot)  [-1]/float( thrown_times[-1])        
                    rate_num_split_per_cl_tot         = np.cumsum(num_split_per_cl_tot)[-1]/float( thrown_times[-1]) 
        
        
                    
             
                    avg_grad_fittest_split_alltimes_allclust =  grad_fittest_split_alltimes_allclust.mean()
                    avg_shortmem_gradx_bulk_split_alltimes_allclust =  shortmem_gradx_bulk_split_alltimes_allclust.mean()
                                
                print avg_vel_clust_split_alltimes_allclust


   #~ data = np.array([pers_l_avg_end2end_fit, err_pers_l_avg_end2end_fit, chisq_avg_end2end_fit, pers_l_avg_end2end_fit_varvel_fixvel, beta_avg_end2end_fit_varvel_fixvel, err_pers_l_avg_end2end_fit_varvel_fixvel, err_beta_avg_end2end_fit_varvel_fixvel, chisq_avg_end2end_fit_varvel_fixvel, diff_avg_angles_subseg_fit_rw_fct_time, err_diff_avg_angles_subseg_fit_rw_fct_time, chisq_avg_angles_subseg_fit_rw_fct_time, pers_l_avg_end2end_fit_constvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel, chisq_avg_end2end_fit_constvel_fct_time_fixvel, pers_l_avg_end2end_fit_varvel_fct_time_fixvel, beta_avg_end2end_fit_varvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel, err_beta_avg_end2end_fit_varvel_fct_time_fixvel, chisq_avg_end2end_fit_varvel_fct_time_fixvel])
    
    
    #~ file_out='{inp}/global_features_pers_l_{real}.txt'.format(inp=dir_out_data, real=real)
       
    
    
                file_in='{inp}/global_features_pers_l_{real}.txt'.format(inp=dir_in, real=real)
                 
                if os.path.isfile(file_in):
                 
                    data_glob = np.loadtxt(file_in)
    
                    if data_glob.ndim >1:
                        data_glob=data_glob[-1]        
        

                    pers_l_avg_end2end_fit                                     = data_glob[0 ]
                    err_pers_l_avg_end2end_fit                                 = data_glob[1 ]
                    chisq_avg_end2end_fit                                      = data_glob[2 ]
                    pers_l_avg_end2end_fit_varvel_fixvel                       = data_glob[3 ]
                    beta_avg_end2end_fit_varvel_fixvel                         = data_glob[4 ]
                    err_pers_l_avg_end2end_fit_varvel_fixvel                   = data_glob[5 ]
                    err_beta_avg_end2end_fit_varvel_fixvel                     = data_glob[6 ]
                    chisq_avg_end2end_fit_varvel_fixvel                        = data_glob[7 ]
                    diff_avg_angles_subseg_fit_rw_fct_time                     = data_glob[8 ]
                    err_diff_avg_angles_subseg_fit_rw_fct_time                 = data_glob[9 ]
                    chisq_avg_angles_subseg_fit_rw_fct_time                    = data_glob[10]
                    pers_l_avg_end2end_fit_constvel_fct_time_fixvel            = data_glob[11]
                    err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel        = data_glob[12]
                    chisq_avg_end2end_fit_constvel_fct_time_fixvel             = data_glob[13]
                    pers_l_avg_end2end_fit_varvel_fct_time_fixvel              = data_glob[14]
                    beta_avg_end2end_fit_varvel_fct_time_fixvel                = data_glob[15]
                    err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel          = data_glob[16]
                    err_beta_avg_end2end_fit_varvel_fct_time_fixvel            = data_glob[17]
                    chisq_avg_end2end_fit_varvel_fct_time_fixvel               = data_glob[18]
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                      
    
    
                file_in='{inp}/global_features_rspeciation_{real}.txt'.format(inp=dir_in, real=real)
                 
                if os.path.isfile(file_in):
                 
                    data_glob = np.loadtxt(file_in)
    
                    if data_glob.ndim >1:
                        data_glob=data_glob[-1]        
        

                    rate_num_split_tot_0d1                                     = data_glob[0 ]
                    rate_num_split_tot_0d5                                     = data_glob[1 ]
                    rate_num_split_tot_1                                     = data_glob[2 ]
                    rate_num_split_tot_10                                     = data_glob[3 ]
                    rate_num_split_per_cl_tot_0d1                                     = data_glob[4 ]
                    rate_num_split_per_cl_tot_0d5                                     = data_glob[5 ]
                    rate_num_split_per_cl_tot_1                                     = data_glob[6 ]
                    rate_num_split_per_cl_tot_10                                     = data_glob[7 ]
 
     
     

    
  
    
    print mem_points, F0, ppl_num, mu, rec_width, time_tot_run, probext, probexpl, numext, numexpl, avg_I, std_I, var_log_I, avg_fit_avg, var_fit_avg, avg_fit_var, var_fit_var, avg_num_IS_coords, avg_num_vir_coords, avg_mean_displ_x, avg_mean_displ_y, avg_var_displ_x, avg_var_displ_y, avg_mean_displ_tot, diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, avg_sig_est, var_sig_est, avg_time_switch, avg_time_fullconv, max_fitn_abserr_switch_max_rel, avg_fitn_abserr_switch_max_rel, frac_1cl, avg_size_clusters_split, avg_n_clusters_split, avg_num_vir_clusters_split, avg_num_coords_clusters_split, avg_num_ext_tot, avg_num_split_tot, avg_num_ext_per_cl_tot, avg_num_split_per_cl_tot, avg_num_vir_o_nclust_tot_list, avg_vel_clust_split_alltimes_allclust, avg_var_parall_clust_split_alltimes_allclust, avg_var_perp_clust_split_alltimes_allclust, avg_var_tot_clust_split_alltimes_allclust, avg_inter_clusters_dist_split, max_inter_clusters_dist_split, max_inter_clusters_dist_max_split, var_vel_clust_split_alltimes_allclust, var_var_parall_clust_split_alltimes_allclust, var_var_perp_clust_split_alltimes_allclust, var_var_tot_clust_split_alltimes_allclust, rate_num_ext_tot, rate_num_split_tot, rate_num_ext_per_cl_tot, rate_num_split_per_cl_tot, num_sharp_turns, rate_sharp_turns, avg_shortmem_err_ys, avg_shortmem_err_svsx, frac_shortmem_err_ys, frac_shortmem_err_svsx, diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, timelast, NUM_VIR, TAU, VEL, S, SIG_SPACE, NOSE_SPACE, R_PERS_L  , avg_grad_fittest_split_alltimes_allclust, avg_shortmem_gradx_bulk_split_alltimes_allclust, avg_S_est, pers_l_avg_end2end_fit, err_pers_l_avg_end2end_fit, chisq_avg_end2end_fit, pers_l_avg_end2end_fit_varvel_fixvel, beta_avg_end2end_fit_varvel_fixvel, err_pers_l_avg_end2end_fit_varvel_fixvel, err_beta_avg_end2end_fit_varvel_fixvel, chisq_avg_end2end_fit_varvel_fixvel, diff_avg_angles_subseg_fit_rw_fct_time, err_diff_avg_angles_subseg_fit_rw_fct_time, chisq_avg_angles_subseg_fit_rw_fct_time, pers_l_avg_end2end_fit_constvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel, chisq_avg_end2end_fit_constvel_fct_time_fixvel, pers_l_avg_end2end_fit_varvel_fct_time_fixvel, beta_avg_end2end_fit_varvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel, err_beta_avg_end2end_fit_varvel_fct_time_fixvel, chisq_avg_end2end_fit_varvel_fct_time_fixvel, rate_num_split_tot_0d1, rate_num_split_tot_0d5, rate_num_split_tot_1, rate_num_split_tot_10, rate_num_split_per_cl_tot_0d1, rate_num_split_per_cl_tot_0d5, rate_num_split_per_cl_tot_1, rate_num_split_per_cl_tot_10
    
    
    
    data_fin = np.array([mem_points, F0, ppl_num, mu, rec_width, time_tot_run, probext, probexpl, numext, numexpl, avg_I, std_I, var_log_I, avg_fit_avg, var_fit_avg, avg_fit_var, var_fit_var, avg_num_IS_coords, avg_num_vir_coords, avg_mean_displ_x, avg_mean_displ_y, avg_var_displ_x, avg_var_displ_y, avg_mean_displ_tot, diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, avg_sig_est, var_sig_est, avg_time_switch, avg_time_fullconv, max_fitn_abserr_switch_max_rel, avg_fitn_abserr_switch_max_rel, frac_1cl, avg_size_clusters_split, avg_n_clusters_split, avg_num_vir_clusters_split, avg_num_coords_clusters_split, avg_num_ext_tot, avg_num_split_tot, avg_num_ext_per_cl_tot, avg_num_split_per_cl_tot, avg_num_vir_o_nclust_tot_list, avg_vel_clust_split_alltimes_allclust, avg_var_parall_clust_split_alltimes_allclust, avg_var_perp_clust_split_alltimes_allclust, avg_var_tot_clust_split_alltimes_allclust, avg_inter_clusters_dist_split, max_inter_clusters_dist_split, max_inter_clusters_dist_max_split, var_vel_clust_split_alltimes_allclust, var_var_parall_clust_split_alltimes_allclust, var_var_perp_clust_split_alltimes_allclust, var_var_tot_clust_split_alltimes_allclust, rate_num_ext_tot, rate_num_split_tot, rate_num_ext_per_cl_tot, rate_num_split_per_cl_tot, num_sharp_turns, rate_sharp_turns, avg_shortmem_err_ys, avg_shortmem_err_svsx, frac_shortmem_err_ys, frac_shortmem_err_svsx, diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, timelast, NUM_VIR, TAU, VEL, S, SIG_SPACE, NOSE_SPACE, R_PERS_L, avg_grad_fittest_split_alltimes_allclust, avg_shortmem_gradx_bulk_split_alltimes_allclust, avg_S_est, pers_l_avg_end2end_fit, err_pers_l_avg_end2end_fit, chisq_avg_end2end_fit, pers_l_avg_end2end_fit_varvel_fixvel, beta_avg_end2end_fit_varvel_fixvel, err_pers_l_avg_end2end_fit_varvel_fixvel, err_beta_avg_end2end_fit_varvel_fixvel, chisq_avg_end2end_fit_varvel_fixvel, diff_avg_angles_subseg_fit_rw_fct_time, err_diff_avg_angles_subseg_fit_rw_fct_time, chisq_avg_angles_subseg_fit_rw_fct_time, pers_l_avg_end2end_fit_constvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel, chisq_avg_end2end_fit_constvel_fct_time_fixvel, pers_l_avg_end2end_fit_varvel_fct_time_fixvel, beta_avg_end2end_fit_varvel_fct_time_fixvel, err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel, err_beta_avg_end2end_fit_varvel_fct_time_fixvel, chisq_avg_end2end_fit_varvel_fct_time_fixvel, rate_num_split_tot_0d1, rate_num_split_tot_0d5, rate_num_split_tot_1, rate_num_split_tot_10, rate_num_split_per_cl_tot_0d1, rate_num_split_per_cl_tot_0d5, rate_num_split_per_cl_tot_1, rate_num_split_per_cl_tot_10])
    
    
    
    
     
    dir_out_data='{inp}/..'.format(inp=dir_io) # directory with input files
     
     
    #file_out='{inp}/summary_analysis_{ppl_num}.dat'.format(inp=dir_out_data, ppl_num=ppl_num)
    file_out=sys.argv[2]
    
    with open(file_out,'a') as f_handle:
        np.savetxt(f_handle, data_fin, fmt='%40.15f', newline=" ")
        f_handle.write("\n")
    
        
     
        
    
