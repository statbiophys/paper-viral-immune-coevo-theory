
"""
COLLECTS SOME AVERAGE STATISTICS TIME SERIES, SUCH AS THE VIRAL POPULATION SIZE

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

import sys
sys.path.append('..')
from lib.mppaper import *
import lib.mpsetup as mpsetup
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from scipy.optimize import fsolve
from scipy.special import factorial
import shutil
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar



def rw_msd_fct_time(x, diff):
    return diff*x



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




def fitness_inv_fct(fitn, M, F_0):
    return  M*(1 -np.exp((fitn - F_0)/M));




def vtau_fromparams(recog_width, mem_points, F_0):
    return (recog_width/(np.exp(F_0/mem_points)-1.))        


dir_io=sys.argv[1] # directory with input files
dir_in_tot='{inp}/realizations'.format(inp=dir_io) # directory with input files
#dir_in = dir_in.replace(".", "d")

dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots
#dir_out_plots = dir_out_plots.replace(".", "d")
  
      
    
    
param_file='{inp}/parameters_backup.dat'.format(inp=dir_io)

params = np.genfromtxt(param_file, dtype="f8,f8,f8,i8,i8,i8,i8, |S10, |S10, |S10, i8, i8, f8", names=['mu','rec_width','jump_size', 'ppl_num', 'maxinfections','save_full_time','n_real','initial_condition','fake_initial_condition','phylo_subsample','mem_points','t_init','F0'])
in_cond=params['initial_condition']
t_init=params['t_init']

print params
print params.shape
n_real=params['n_real']
num_ppl=params['ppl_num']
recog_width=params['rec_width']
mem_points=params['mem_points']
F0=params['F0']

latt_sp=1
vtau=vtau_fromparams(recog_width, mem_points, F0)


print in_cond
#except ValueError:
#    params = np.genfromtxt(param_file, dtype="f8,f8,f8,f8,f8,i8,i8,f8,i8,i8,i8", names=['mu','rec_width','sigma','eps','I','ppl_num','save_cov','f_m','maxinfections','save_full_time','n_real'])
    

print params
print params.shape
n_real=params['n_real']
ppl_num=params['ppl_num']
print n_real

    
thisfigsize = figsize
thisfigsize[1] *= 0.75
      
times_reals                      = []
num_IS_coords_reals                   = []
num_vir_coords_reals                   = []
num_IS_tot_reals                   = []
num_vir_tot_reals                   = []
avg_fitn_reals                   = []
var_fitn_reals                   = []

# get data
for real in np.arange(1,n_real+1):
    
    dir_in='{inp}/realization_{real}'.format(inp=dir_in_tot,real=real) # directory with input files
    dir_out_data=dir_in
    dir_out_plots='{inp}/realization_{real}'.format(inp=dir_out_plots_tot,real=real) # directory with output plots

    file_exploded='{inp}/expl_file.txt'.format(inp=dir_in)
    file_extinct='{inp}/extinct_file.txt'.format(inp=dir_in)
    
#    if os.path.exists(dir_out_plots):
#   shutil.rmtree(dir_out_plots, ignore_errors=False, onerror=None)
       
    if not os.path.exists(dir_out_plots):
        os.makedirs(dir_out_plots) 
        


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
            
            
            time_mod=time[time>10000] # throw first 100 years
            
            sec_thr=min(100000, time_mod[time_mod.size/4])
            
            time_mod=time_mod[time_mod>sec_thr]
            
            print sec_thr, time.size
            
        
            num_vir_tot       = data[:,7]
            
            num_vir_tot_thr=np.mean(num_vir_tot)
            if os.path.isfile(file_exploded):
                ind_times_notexpl=num_vir_tot.size - np.argmax(num_vir_tot[::-1] < num_vir_tot_thr)
            elif os.path.isfile(file_extinct):
                ind_times_notexpl=num_vir_tot.size - np.argmax(num_vir_tot[::-1] > num_vir_tot_thr)
            else:
                ind_times_notexpl=num_vir_tot.size
            
            print ind_times_notexpl, num_vir_tot.size
            
            time_all = time[:ind_times_notexpl]
            
            time_ss=time_all[time_all>sec_thr]
    
            
            #survival=time_all[-1]
            #time_ss=time[(time> t_init/365.) & (time <= survival)]
            #
            #print survival, t_init/365.
            print time_ss.size
            
            if time_ss.size > 3000:
                
                data= data[:ind_times_notexpl,:]
                data= data[time_all>sec_thr,:]
            
                #evo_mean_file<<"# 0 time"<<setw(30)<<" 1 num_x"<<setw(30)<<" 2 num_y "<<setw(30)<<" 3 area"<<setw(30)<<"4 number IS coord"<<setw(30)<<"5 number viruses coord"<<setw(30)<<" 6 tot number IS"<<setw(30)<<"7 tot number viruses"<<setw(30)<<"8 avg viral fitness"<<setw(30)<<"9 var viral fitness" <<setw(30)<<"10 start_x" <<setw(30)<<"11 start_y" <<setw(30)<<"12 avg_vir_x" <<setw(30)<<"13 avg_vir_y" <<endl;
                 
                   
            
                time = data[:,0]
                
                
           # evo_mean_file<<"# 0 time"<<setw(30)<<" 1 num_x"<<setw(30)<<" 2 num_y "<<setw(30)<<" 3 area"<<setw(30)<<"4 number IS coord"<<setw(30)<<"5 number viruses coord"<<setw(30)<<" 6 tot number IS"<<setw(30)<<"7 tot number viruses"<<setw(30)<<"8 avg viral fitness"<<setw(30)<<"9 var viral fitness" <<setw(30)<<"10 start_x" <<setw(30)<<"11 start_y" <<setw(30)<<"12 avg_vir_x" <<setw(30)<<"13 avg_vir_y"<<setw(30)<<"14 maxfit" <<setw(30)<<"15 number IS coord update" <<setw(30)<<"16 tot number IS update"  <<setw(30)<<"17 min vir x"   <<setw(30)<<"18 min vir y"   <<setw(30)<<"19 max vir x"   <<setw(30)<<"20 max vir y"   <<setw(30)<<"21 min IS x "   <<setw(30)<<"22 min IS y "   <<setw(30)<<"23 max IS x "   <<setw(30)<<"24 max IS y "   <<setw(30)<<"25 mean_displ_x_ "   <<setw(30)<<"26 mean_displ_y_ "   <<setw(30)<<"27 mean_displ_tot_ "   <<setw(30)<<"28 var_displ_x_ "   <<setw(30)<<"29 var_displ_y_ "   <<setw(30)<<"30 count_displ_ "   <<setw(30)<<"31 x1_max_fitn "   <<setw(30)<<"32 y1_max_fitn " <<endl;
        

        
                
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
                    
                   
        
                vel_cloud_avg = np.sqrt((vir_x[1:] - vir_x[:-1])**2 + (vir_y[1:] - vir_y[:-1])**2)/(time[1:] - time[:-1])
                
                end2end_coarse_gr_fct_time_list = []
                time_diff_traj_list             = []
                ttot=time[-1] - time[0]
                    
                    
                # OUTDATED OBSERVABLES
                diff_maxfit_fit_rw_fct_time    = 0
                err_diff_maxfit_fit_rw_fct_time   = 0
                chisq_maxfit_fit_rw_fct_time   = 0
                
                
                
                
                var_log_I= np.var(np.log(num_vir_tot/float(ppl_num)))
                
                avg_I=np.mean(num_vir_tot/float(ppl_num))
                std_I=np.std(num_vir_tot/float(ppl_num))
                
                avg_fit_avg=np.mean(avg_fitn)
                var_fit_avg=np.var(avg_fitn)
                
                avg_fit_var=np.mean(var_fitn)
                var_fit_var=np.var(var_fitn)
                
                avg_fit_nose=np.mean(fit_nose)
                var_fit_nose=np.var(fit_nose)
                
                avg_vel=np.mean(vel_cloud_avg[1:])
                var_vel=np.var(vel_cloud_avg[1:])
                
                
                avg_num_IS_coords=np.mean(num_IS_coords)
                avg_num_vir_coords=np.mean(num_vir_coords)
                
                single_mean_displ_x = mean_displ_x[count_displ==1]
                single_mean_displ_y = mean_displ_y[count_displ==1]
                single_mean_displ_tot = np.sqrt(np.sum(mean_displ_x)**2 + np.sum(mean_displ_y)**2)
                
                single_count_displ= np.sum(count_displ==1)
                
                mean_single_mean_displ_x = single_mean_displ_x.mean()
                mean_single_mean_displ_y = single_mean_displ_y.mean()
                var_single_mean_displ_x = single_mean_displ_x.var()
                var_single_mean_displ_y = single_mean_displ_y.var()
                
                new_count_displ    = np.append(count_displ,single_count_displ)
                new_mean_displ_tot = np.append(mean_displ_tot,single_mean_displ_tot)
                new_mean_displ_x   = np.append(mean_displ_x  ,mean_single_mean_displ_x)
                new_mean_displ_y   = np.append(mean_displ_y  ,mean_single_mean_displ_y)
                new_var_displ_x    = np.append(var_displ_x   ,var_single_mean_displ_x )
                new_var_displ_y    = np.append(var_displ_y   ,var_single_mean_displ_y )
                
                
                avg_mean_displ_tot = np.mean(new_mean_displ_tot[new_count_displ>3])
                avg_mean_displ_x   = np.mean(new_mean_displ_x  [new_count_displ>3])
                avg_mean_displ_y   = np.mean(new_mean_displ_y  [new_count_displ>3])
                avg_var_displ_x    = np.mean(new_var_displ_x   [new_count_displ>3])
                avg_var_displ_y    = np.mean(new_var_displ_y   [new_count_displ>3])
                
                avg_sig_est    = np.mean(sig_est )
                var_sig_est    = np.var(sig_est )
                
                #data = np.array([t_p, t_del, t_dup, t_ch, p1, p2, p3, p4, p5, avg_pid_len1, avg_pid_len10, avg_pid_len18, avg_pid_len2_term, avg_pid_len18_term])
                data = np.array([avg_I, std_I, var_log_I, avg_fit_avg, var_fit_avg, avg_fit_var, var_fit_var, avg_fit_nose, var_fit_nose, avg_vel, var_vel, avg_num_IS_coords, avg_num_vir_coords, avg_mean_displ_x, avg_mean_displ_y, avg_var_displ_x, avg_var_displ_y, avg_mean_displ_tot, diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, avg_sig_est, var_sig_est])
                
                
                file_out='{inp}/global_features_{real}.txt'.format(inp=dir_out_data, real=real)
            
                with open(file_out,'a') as f_handle:
                    np.savetxt(f_handle, data, fmt='%15.15f', newline=" ")
                    f_handle.write("\n")
            
        
            
                  
                if real < 10:
                       
                    
                                
                    
                    fig = plt.figure(figsize=thisfigsize)
                    grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.2, wspace=0.4, hspace=0.35)
                    labeled_axes = []
                        
                    ax = plt.Subplot(fig, grid[0, 0])
                    fig.add_subplot(ax)
                    labeled_axes.append(ax)
                    ax.plot(time[1:-1], vel_cloud_avg[1:], linestyle='-', color='g')
                    ax.set_xlabel('time (cy)')
                    ax.set_ylabel('average virus speed')
                    ax.xaxis.labelpad = axis_labelpad
                    ax.yaxis.labelpad = axis_labelpad
                    mpsetup.despine(ax) 
                        
                        
                    #### finish figure ####
                    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                    xycoords=('axes fraction'), fontweight = 'bold')
                    #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
                    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
                    out_file='{out}/vir_avg_speed_{real}.pdf'.format(out=dir_out_plots, real=real)
                    #    print out_file
                    fig.savefig(out_file)
                    #    plt.show()
                    #    print out_file
                        
                        
                        
                        
                        
                        
                        
                    fig = plt.figure(figsize=thisfigsize)
                    grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.2, wspace=0.4, hspace=0.35)
                    labeled_axes = []
                        
                    ax = plt.Subplot(fig, grid[0, 0])
                    fig.add_subplot(ax)
                    labeled_axes.append(ax)
                    ax.plot(time[count_displ>3], mean_displ_x[count_displ>3], linestyle='-', color='g', label='x')
                    ax.plot(time[count_displ>3], mean_displ_y[count_displ>3], linestyle='-', color='r', label='y')
                    ax.set_xlabel('time ')
                    ax.set_ylabel('mean displacement ')
                    ax.xaxis.labelpad = axis_labelpad
                    ax.yaxis.labelpad = axis_labelpad
                    mpsetup.despine(ax) 
                        
                        
                    #### finish figure ####
                    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                    xycoords=('axes fraction'), fontweight = 'bold')
                    #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
                    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
                    out_file='{out}/mean_displ_xy_{real}.pdf'.format(out=dir_out_plots, real=real)
                    #    print out_file
                    fig.savefig(out_file)
                    #    plt.show()
                    #    print out_file
                    
                        
                        
                        
                        
                        
                        
                    fig = plt.figure(figsize=thisfigsize)
                    grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.2,
                             wspace=0.4, hspace=0.35)
                    labeled_axes = []
                        
                    ax = plt.Subplot(fig, grid[0, 0])
                    fig.add_subplot(ax)
                    labeled_axes.append(ax)
                    ax.plot(time, num_vir_coords, linestyle='-', color='g')
                    ax.set_xlabel('time (y)')
                    ax.set_ylabel('virus coordinates')
                    ax.xaxis.labelpad = axis_labelpad
                    ax.yaxis.labelpad = axis_labelpad
                    mpsetup.despine(ax) 
                        
                        
                    #### finish figure ####
                    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                         xycoords=('axes fraction'), fontweight = 'bold')
                    #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
                    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
                    out_file='{out}/num_vir_coords_{real}.pdf'.format(out=dir_out_plots, real=real)
                    #    print out_file
                    fig.savefig(out_file)
                    #    plt.show()
                    #    print out_file
                        
                        
                        
                        
                        
                    ## plot frac_inf    
                        
                    fig = plt.figure(figsize=thisfigsize)
                    grid = gridspec.GridSpec(1, 2, left=0.2, right=0.97, top=0.91, bottom=0.22,
                             wspace=0.4, hspace=0.35)
                    labeled_axes = []
                        
                    ax = plt.Subplot(fig, grid[0, 0])
                    fig.add_subplot(ax)
                    labeled_axes.append(ax)
                    ax.plot(time, num_vir_tot, linestyle='-', color='g')
                    ax.set_xlabel('time (y)')
                    ax.set_ylabel('number viruses')
                    ax.set_yscale('log')
                    ax.xaxis.labelpad = axis_labelpad
                    ax.yaxis.labelpad = axis_labelpad
                    mpsetup.despine(ax) 
                    
                    
                        
                        
                    #### finish figure ####
                    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                        xycoords=('axes fraction'), fontweight = 'bold')
                    #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
                    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
                    out_file='{out}/num_vir_tot_log_{real}.pdf'.format(out=dir_out_plots, real=real)
                    #    print out_file
                    fig.savefig(out_file)
                    #    plt.show()
                    #    print out_file
                    
                    
                    
                    
                        
                        
                    ## plots fig 3A
                        
                    fig = plt.figure(figsize=(thisfigsize[0]/2.,thisfigsize[0]/2.))
                    grid = gridspec.GridSpec(1, 1, left=0.12, right=0.97, top=0.93, bottom=0.13,
                             wspace=0.4, hspace=0.35)
                    labeled_axes = []
                    ax = plt.Subplot(fig, grid[0, 0])
                    fig.add_subplot(ax)
                    labeled_axes.append(ax)
                    ax.plot(vir_x, vir_y, linestyle='-', color='g')#, label='average virus'
                    #ax.plot(IS_x, IS_y, linestyle='-', color='r', label='average IS')
            
                    frame1 = plt.gca()
                    #~ frame1.axes.get_xaxis().set_visible(False)
                    #~ frame1.axes.get_yaxis().set_visible(False)
                    frame1.axes.get_xaxis().set_ticks([])
                    frame1.axes.get_yaxis().set_ticks([])
                    ax.set_xlabel('phenotypic trait 1')
                    ax.set_ylabel('phenotypic trait 2')
                    
                        
                                
                                
                    ylims_diff = ax.get_ylim()[1] - ax.get_ylim()[0]
                    xlims_diff = ax.get_xlim()[1] - ax.get_xlim()[0]
                    
                    print xlims_diff, ylims_diff
                    
                    
                    if ylims_diff > xlims_diff:
                        
                        diff_diff= ylims_diff-xlims_diff
                        
                        maxdiff=ylims_diff
                        
                        #~ xmin=ax.get_xlim()[0] - diff_diff/2.
                        #~ xmax=ax.get_xlim()[1] + diff_diff/2.
                        xmin=ax.get_xlim()[0]  - diff_diff
                        xmax=ax.get_xlim()[1]
                        
                        ax.set_xlim(xmin, xmax)
                        
                    
                    elif xlims_diff > ylims_diff:
                    
                        diff_diff= xlims_diff-ylims_diff
                        
                        #~ ymin=ax.get_ylim()[0] - diff_diff/2.
                        #~ ymax=ax.get_ylim()[1] + diff_diff/2.
                        ymin=ax.get_ylim()[0]
                        ymax=ax.get_ylim()[1]  + diff_diff
                        
                        ax.set_ylim(ymin, ymax)
                        maxdiff=xlims_diff
                        
                        
                    
                    
                    #ax.set_xlim(xmin, xmax)
                    #ax.set_ylim(ymin, ymax)
                    
                    print maxdiff, recog_width
                    
                    if maxdiff>200*recog_width:
                        scale_size=100*recog_width
                        scale_lab=r'$100 r$'
                    elif maxdiff>20*recog_width:
                        scale_size=10*recog_width
                        scale_lab=r'$10 r$'
                    elif maxdiff>2*recog_width:
                        scale_size=recog_width
                        scale_lab=r'$r$'
                    elif maxdiff>recog_width:
                        scale_size=0.5*recog_width
                        scale_lab=r'$0.5 r$'
                    elif maxdiff>0.2*recog_width:
                        scale_size=0.1*recog_width
                        scale_lab=r'$0.1 r$'
                    elif maxdiff>0.1*recog_width:
                        scale_size=0.05*recog_width
                        scale_lab=r'$0.05 r$'
                    
                    print scale_size, scale_lab
                            
                    scalebar = AnchoredSizeBar(ax.transData,
                                   scale_size, scale_lab, 2, 
                                   pad=0.1,
                                   color='black',
                                   frameon=False,
                                   size_vertical=0.,
                                   label_top=True)# ,fontproperties=fontprops
                                   
                    
                    ax.add_artist(scalebar)
                    
                    
                    
                    if maxdiff>200*vtau:
                        scale_size=100*vtau
                        scale_lab=r'$100 v\tau$'
                    elif maxdiff>40*vtau:
                        scale_size=20*vtau
                        scale_lab=r'$20 v\tau$'
                    elif maxdiff>20*vtau:
                        scale_size=10*vtau
                        scale_lab=r'$10 v\tau$'
                    elif maxdiff>2*vtau:
                        scale_size=vtau
                        scale_lab=r'$v\tau$'
                    elif maxdiff>vtau:
                        scale_size=0.5*vtau
                        scale_lab=r'$0.5 v\tau$'
                    elif maxdiff>0.2*vtau:
                        scale_size=0.1*vtau
                        scale_lab=r'$0.1 v\tau$'
                    elif maxdiff>0.1*vtau:
                        scale_size=0.05*vtau
                        scale_lab=r'$0.05 v\tau$'
                    print maxdiff, vtau
                    print scale_size, scale_lab
                            
                    scalebar = AnchoredSizeBar(ax.transData,
                                   scale_size, scale_lab, 3, 
                                   pad=0.1,
                                   color='black',
                                   frameon=False,
                                   size_vertical=0.,
                                   label_top=True)# ,fontproperties=fontprops
                                   
                    
                    ax.add_artist(scalebar)
                    
                    
                    
                    ax.xaxis.labelpad = axis_labelpad
                    ax.yaxis.labelpad = axis_labelpad
                    ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
                          loc='upper right', bbox_to_anchor=(2.05, 1.18))
                    mpsetup.despine(ax) 
                        
                   # ax.set_aspect('equal', 'box')  
                    ax.axis('equal') 
                    #### finish figure ####
                    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                         xycoords=('axes fraction'), fontweight = 'bold')
                    #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
                    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
                    out_file='{out}/trajs_cartoon_{real}.pdf'.format(out=dir_out_plots, real=real)
                    #    print out_file
                    fig.savefig(out_file)
                        
                        
                        
                        
                        
                                    
                    timestep=1000
                 
                    file_in_FFTerr_npz_compr='{inp}/times_errors_benchmark_conv.dat'.format(inp=dir_in)
                    
                    if os.path.isfile(file_in_FFTerr_npz_compr) and sum(1 for line in open(file_in_FFTerr_npz_compr))>2:
        
                        #data_FFTerr = np.loadtxt(file_in_FFTerr_npz_compr)
                        data_FFTerr=[]
                        with open(file_in_FFTerr_npz_compr) as f:
                            lines=f.readlines()
                            for line in lines:
                                if not line.startswith("#"):
                                    line=line.decode('utf-8','ignore').encode("utf-8")
                                    myarray = np.fromstring(line, dtype=float, sep=' ')
                                    data_FFTerr.append(myarray)
                                    #print(myarray)
                            print len(data_FFTerr)
                            
                            data_FFTerr = np.asarray(data_FFTerr)
                            print data_FFTerr.shape
                  
                    
                        
                        time_nfft       = data_FFTerr[:,0]
                        time_farf       = data_FFTerr[:,1]
                        time_farf_nfft       = data_FFTerr[:,2]
                        time_fullconv       = data_FFTerr[:,3]
                        
                        err_fullconv_tot       = data_FFTerr[:,4]
                        
                        #err_nfft_tot_rel       = np.abs(data_FFTerr[:,5] - err_fullconv_tot)/err_fullconv_tot
                        #err_farf_tot_rel       = np.abs(data_FFTerr[:,6] - err_fullconv_tot)/err_fullconv_tot
                        #err_farf_nfft_tot_rel  = np.abs(data_FFTerr[:,7] - err_fullconv_tot)/err_fullconv_tot
                        
                        err_nfft_tot_rel       = (data_FFTerr[:,5] - err_fullconv_tot)/err_fullconv_tot
                        err_farf_tot_rel       = (data_FFTerr[:,6] - err_fullconv_tot)/err_fullconv_tot
                        err_farf_nfft_tot_rel  = (data_FFTerr[:,7] - err_fullconv_tot)/err_fullconv_tot
                        
                        err_nfft_max_rel       = data_FFTerr[:,8] 
                        err_farf_max_rel       = data_FFTerr[:,9] 
                        err_farf_nfft_max_rel  = data_FFTerr[:,10] 
        
                        err_fullconv_tot_avg       = data_FFTerr[:,11]
                        
                        err_nfft_tot_avg_rel       = (data_FFTerr[:,12] - err_fullconv_tot_avg)/err_fullconv_tot_avg
                        err_farf_tot_avg_rel       = (data_FFTerr[:,13] - err_fullconv_tot_avg)/err_fullconv_tot_avg
                        err_farf_nfft_tot_avg_rel  = (data_FFTerr[:,14] - err_fullconv_tot_avg)/err_fullconv_tot_avg
                        
                        
                        conv_err_nfft_max_rel       = data_FFTerr[:,15] 
                        conv_err_farf_max_rel       = data_FFTerr[:,16] 
                        conv_err_farf_nfft_max_rel  = data_FFTerr[:,17] 
                        
                        fitn_avg_abserr_nfft_max_rel       = data_FFTerr[:,18] 
                        fitn_avg_abserr_farf_max_rel       = data_FFTerr[:,19] 
                        fitn_avg_abserr_farf_nfft_max_rel  = data_FFTerr[:,20] 
                        
                        
                        conv_avg_abserr_nfft_max_rel       = data_FFTerr[:,21] 
                        conv_avg_abserr_farf_max_rel       = data_FFTerr[:,22] 
                        conv_avg_abserr_farf_nfft_max_rel  = data_FFTerr[:,23] 
                        
                        
                        conv_resc_avg_abserr_nfft_max_rel       = data_FFTerr[:,24] 
                        conv_resc_avg_abserr_farf_max_rel       = data_FFTerr[:,25] 
                        conv_resc_avg_abserr_farf_nfft_max_rel  = data_FFTerr[:,26] 
                        
                        
                        
                        fitn_fullconv_stdev_bench  = data_FFTerr[:,27] 
                        conv_fullconv_stdev_bench  = data_FFTerr[:,28] 
                        
                        
                        time_switch              = data_FFTerr[:,29]
                        err_switch_tot_rel       = (data_FFTerr[:,30] - err_fullconv_tot)/err_fullconv_tot
                        err_switch_max_rel       = data_FFTerr[:,31] 
                        err_switch_tot_avg_rel   = (data_FFTerr[:,32] - err_fullconv_tot_avg)/err_fullconv_tot_avg
                        conv_err_switch_max_rel  = data_FFTerr[:,33] 
                        fitn_avg_abserr_switch_max_rel  = data_FFTerr[:,34] 
                        conv_avg_abserr_switch_max_rel  = data_FFTerr[:,35] 
                        conv_resc_avg_abserr_switch_max_rel  = data_FFTerr[:,36] 
                        
                        
                        
                        mask=np.abs(fitn_fullconv_stdev_bench)>10**(-14)
                        filt_fitn_avg_abserr_switch = fitn_avg_abserr_switch_max_rel[mask]/fitn_fullconv_stdev_bench[mask]
                                
                        
                        mask=np.abs(conv_fullconv_stdev_bench)>10**(-14)
                        filt_conv_avg_abserr_switch = conv_avg_abserr_switch_max_rel[mask]/conv_fullconv_stdev_bench[mask]
                                
                        max_fitn_abserr_switch_max_rel    = np.amax(filt_fitn_avg_abserr_switch )
                        avg_fitn_abserr_switch_max_rel    = np.mean(filt_fitn_avg_abserr_switch )
                        avg_conv_abserr_switch_max_rel    = np.mean(filt_conv_avg_abserr_switch )
                        avg_conv_err_switch_max_rel    = np.mean(conv_err_switch_max_rel )
                        avg_time_switch    = np.mean(time_switch )
                        avg_time_fullconv    = np.mean(time_fullconv )
                        
                        #data = np.array([t_p, t_del, t_dup, t_ch, p1, p2, p3, p4, p5, avg_pid_len1, avg_pid_len10, avg_pid_len18, avg_pid_len2_term, avg_pid_len18_term])
                        data = np.array([avg_time_switch, avg_time_fullconv, max_fitn_abserr_switch_max_rel, avg_fitn_abserr_switch_max_rel, avg_conv_abserr_switch_max_rel, avg_conv_err_switch_max_rel])
                        
                        
                        file_out='{inp}/global_features_bench_fitupd_{real}.txt'.format(inp=dir_out_data, real=real)
                    
                        with open(file_out,'a') as f_handle:
                            np.savetxt(f_handle, data, fmt='%15.15f', newline=" ")
                            f_handle.write("\n")
                    
        
            
              
                    
                                    
                        fig = plt.figure(figsize=thisfigsize)
                        grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.2, wspace=0.4, hspace=0.35)
                        labeled_axes = []
                        
                        #print diff_fitn_maxIS[np.abs(avg_fitn)>10**(-5)]/np.abs(avg_fitn)[np.abs(avg_fitn)>10**(-5)]
                            
                        ax = plt.Subplot(fig, grid[0, 0])
                        fig.add_subplot(ax)
                        labeled_axes.append(ax)
                        
           
                        mask=np.abs(fitn_fullconv_stdev_bench)>10**(-14)
                        
                        x=np.arange(fitn_fullconv_stdev_bench.size)*timestep
                        x= x[mask]
                        #denom=np.abs(avg_fitn)[-diff_fitn_maxfit.size:]
                        #denom=denom[mask]
                        
                        ax.plot(x, fitn_avg_abserr_nfft_max_rel    [mask]/fitn_fullconv_stdev_bench[mask], linestyle='-', color='b', label='nfft')
                        ax.plot(x, fitn_avg_abserr_farf_nfft_max_rel[mask]/fitn_fullconv_stdev_bench[mask], linestyle='-', color='r', label='far field plus nfft')
                        ax.plot(x, fitn_avg_abserr_switch_max_rel[mask]/fitn_fullconv_stdev_bench[mask], linestyle='-', color='c', label='switch')
                        ax.plot(x, fitn_avg_abserr_farf_max_rel    [mask]/fitn_fullconv_stdev_bench[mask], linestyle='-', color='g', label='far field')
                        ax.set_xlabel('time ')
                        ax.set_ylabel('average fitness error o stddev')
                        ax.set_yscale('log')
                        ax.xaxis.labelpad = axis_labelpad
                        ax.yaxis.labelpad = axis_labelpad
                        mpsetup.despine(ax) 
                        ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
                          loc='upper right', bbox_to_anchor=(1.05, 1.))
        
                            
                        #### finish figure ####
                        labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                        xycoords=('axes fraction'), fontweight = 'bold')
                        #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
                        #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
                        out_file='{out}/fitn_avg_abserr_o_std_benchmark_{real}.pdf'.format(out=dir_out_plots, real=real)
                        #    print out_fconv_ile
                        fig.savefig(out_file)
                        #    plt.show()
                        #    print out_file
                        
                        
                        
                    
                    
                        
                        
                        
                        
                    fig.clf()
                    plt.close('all')
                    
                        
                        
