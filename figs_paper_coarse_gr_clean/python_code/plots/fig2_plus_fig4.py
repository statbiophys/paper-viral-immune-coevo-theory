####
# Figure 3
# needs:
# - data/1d.npz produced by intra_host_neut.py
# - data/2d.npz produced by run2d.py
####
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
#import shutil
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
import matplotlib.tri as tri
import matplotlib.ticker as mticker
from scipy.optimize import root
# ~ import cmocean
# ~ from cmcrameri import cm as cmcram
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit



f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
fmt = mticker.FuncFormatter(g)

#### subfigure A ####
## import data

def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]
    
    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)

def show_colormap(cmap):
    im = np.outer(np.ones(10), np.arange(100))
    fig, ax = plt.subplots(2, figsize=(6, 1.5),
                           subplot_kw=dict(xticks=[], yticks=[]))
    fig.subplots_adjust(hspace=0.1)
    ax[0].imshow(im, cmap=cmap)
    ax[1].imshow(im, cmap=grayify_cmap(cmap))
    
    out_file='{out}/cmap.pdf'.format(out=dir_out_plots_scatter)
    #    print out_file, dpi=300
    fig.savefig(out_file)
    


def get_cmap(N):
    ''' Returns a function that maps each index in 0, 1, ...
        N-1 to a distinct RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    #scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap='plasma')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color




def myExpFunc(x, a, b):
    return a * np.exp(-b*x)
    
    
    

	    

def plot_scatter_colorbar(array_x, array_y, xlabel, ylabel, ax=None, show_cbar=True, fig=None, show_legend=True): # between two observables
    global m_list   
    global f_m_list   
    global pplnum_list
    
    show_cbar_inset=False
    
    if 'fig2' in xlabel:
        show_cbar=False
        show_legend=False
        if "num virs teo" in ylabel:
            show_cbar_inset=True

    #~ markers = np.asarray(["x","+","*"])
    markers = np.asarray(["s","o","v"])
    sizes = np.asarray([7,14])

    array_cycle_ppl=ppl_num
    array_cycle_F0=F0
    array_cycle_M=mem_points
   
    array_z_orig=mu
    array_z=mu
    if 'fig2' in xlabel:
        array_z=D_mod_all/(real_rec_width**2.) 

    #~ array_z=SIG_SPACE
    #~ array_z=S_TOT
    #~ array_z=D_mod_all
    #~ array_z=avg_I
    #array_z=frac_mutants
    zlabel=r'$\mu$ '
    

    
    if array_y.size > array_z.size:
        array_z_orig=all_mu
        array_z=all_mu
        #~ array_z=S_TOT
        array_cycle_ppl=all_ppl_num
        array_cycle_F0 =all_F0
        array_cycle_M =all_mem_points
 
    if ("clust" in xlabel or "lineage"  in xlabel )  and "fails" not in xlabel  and "all " not in xlabel :
 
        #~ array_x=rec_width[~fails_clust.astype(bool)]
        array_z_orig=mu[~fails_clust.astype(bool)]
        array_z=mu[~fails_clust.astype(bool)]
        if 'fig2' in xlabel:
            array_z=D_mod_all[~fails_clust.astype(bool)]/(real_rec_width [~fails_clust.astype(bool)]**2.) 
        #~ array_z=S_TOT[~fails_clust.astype(bool)]
        
        array_cycle_ppl=ppl_num[~fails_clust.astype(bool)]
        array_cycle_F0 = F0[~fails_clust.astype(bool)]
        array_cycle_M = mem_points[~fails_clust.astype(bool)]
        
        
            

    if "pers time" in xlabel and "fails" not in xlabel   :
        #~ array_x=rec_width[~fails_persl.astype(bool)]
        array_z_orig=mu[~fails_persl.astype(bool)]
        array_z=mu[~fails_persl.astype(bool)]
    
        array_cycle_ppl=ppl_num[~fails_persl.astype(bool)]
        array_cycle_F0 = F0[~fails_persl.astype(bool)]
        array_cycle_M = mem_points[~fails_persl.astype(bool)]
        if 'fig2' in xlabel:
            array_z=D_mod_all[~fails_persl.astype(bool)]/(real_rec_width [~fails_persl.astype(bool)]**2.) 


  
    if "altperscolor" in xlabel:
        # ~ array_z=timelast[~fails_persl.astype(bool)]
        # ~ array_z=timescale_system
        array_z=timescale_system/timelast[~fails_persl.astype(bool)].astype(float)
        # ~ array_z=chisq_avg_angles_subseg_fit_rw_fct_time
    if "altperscolor2" in xlabel:
        # ~ array_z=timelast[~fails_persl.astype(bool)]
        # ~ array_z=timescale_system
        # ~ array_z=timescale_system/timelast[~fails_persl.astype(bool)].astype(float)
        array_z=chisq_avg_angles_subseg_fit_rw_fct_time
        
    if "altperscolor3" in xlabel:
        # ~ array_z=timelast[~fails_persl.astype(bool)]
        # ~ array_z=timescale_system
        # ~ array_z=timescale_system/timelast[~fails_persl.astype(bool)].astype(float)
        # ~ array_z=np.abs(pers_l_avg_end2end_fit_varvel_fct_time_fixvel - pers_l_avg_end2end_fit_constvel_fct_time_fixvel)/pers_l_avg_end2end_fit_varvel_fct_time_fixvel
        array_z=np.abs(pers_l_avg_end2end_fit_varvel_fct_time_fixvel - array_y)/array_y
        
        

  
    if "altclustcolor" in xlabel:
        # ~ array_z=timelast[~fails_persl.astype(bool)]
        # ~ array_z=timescale_system
        array_z=timescale_system
        zlabel=r'$T$ '
        # ~ array_z=chisq_avg_angles_subseg_fit_rw_fct_time
    if "altclustcolor2" in xlabel:
        # ~ array_z=timelast[~fails_persl.astype(bool)]
        # ~ array_z=timescale_system
        # ~ array_z=timescale_system/timelast[~fails_persl.astype(bool)].astype(float)
        array_z=S
        zlabel=r'$S$ '
        
    if "altclustcolor3" in xlabel:
        # ~ array_z=timelast[~fails_persl.astype(bool)]
        # ~ array_z=timescale_system
        # ~ array_z=timescale_system/timelast[~fails_persl.astype(bool)].astype(float)
        # ~ array_z=np.abs(pers_l_avg_end2end_fit_varvel_fct_time_fixvel - pers_l_avg_end2end_fit_constvel_fct_time_fixvel)/pers_l_avg_end2end_fit_varvel_fct_time_fixvel
        array_z=SIG_SPACE
        zlabel=r'$\sigma$ '
        
        
    if "altclustcolor4" in xlabel:
        # ~ array_z=timelast[~fails_persl.astype(bool)]
        # ~ array_z=timescale_system
        # ~ array_z=timescale_system/timelast[~fails_persl.astype(bool)].astype(float)
        # ~ array_z=np.abs(pers_l_avg_end2end_fit_varvel_fct_time_fixvel - pers_l_avg_end2end_fit_constvel_fct_time_fixvel)/pers_l_avg_end2end_fit_varvel_fct_time_fixvel
        array_z=avg_I[~fails_clust.astype(bool)]
        zlabel=r'$I$ '
        
        


    if 'fig2' in xlabel:
        array_cycle_ppl=array_cycle_ppl[(array_z_orig >=0.01)]
        array_cycle_F0 = array_cycle_F0[(array_z_orig >=0.01)]
        array_cycle_M = array_cycle_M[(array_z_orig >=0.01)]
        array_x = array_x[(array_z_orig >=0.01)]
        array_y = array_y[(array_z_orig >=0.01)]
        array_z=array_z[(array_z_orig >=0.01)]
        
                    
    mask=np.isfinite(array_x) & np.isfinite(array_y)
    mask=mask & (array_y>0) & (array_x>0)


    if "pers time" in xlabel and "angle" in xlabel:
        chisq_mask=chisq_avg_angles_subseg_fit_rw_fct_time   
        if 'fig2' in xlabel: 
            chisq_mask=chisq_avg_angles_subseg_fit_rw_fct_time[(array_z_orig >=0.01)]
            
    if "pers time" in xlabel and "constvel" in xlabel:
        chisq_mask=chisq_avg_end2end_fit_constvel_fct_time_fixvel  
    if "pers time" in xlabel and "varvel" in xlabel:
        chisq_mask=chisq_avg_end2end_fit_varvel_fct_time_fixvel   
    
    if "pers time" in xlabel or "turn" in xlabel:
        # ~ mask=mask & (chisq_mask < 1000)  & (chisq_mask > 0.4)  
        mask=mask & (array_x>0)  
        # ~ if  "angle" in xlabel:
            # ~ mask=mask  & (array_x<100000)   
        # ~ else:
            # ~ mask=mask  & (array_x<10000)   
        if 'filterpersl' in xlabel:
            mask=mask & (chisq_mask < 3) 
    
    array_z=array_z[mask]
    X=array_x[mask]
    Y=array_y[mask]  
    
    array_cycle_ppl=array_cycle_ppl[mask]
    array_cycle_F0 =array_cycle_F0 [mask]      
    array_cycle_M =array_cycle_M [mask]      
        
        
    array_col_m    =np.unique(array_z)
 


    mincol=np.amin(array_z)
    maxcol=np.amax(array_z) 
    
    if 'fig2' in xlabel:
        mincol= 10.**(-13.)
        maxcol= 2.*10.**(-6.)
 
    
    # ~ if 'fig2' in xlabel:
        # ~ mincol= 0.0005
        # ~ maxcol= 0.01
 
    print array_x
    print array_y
    print array_y[(array_x>0)]
    # ~ print D_mod_all
    print array_z
    print array_col_m
    
    print X
    print Y
    print np.amin(Y)
    print np.amax(Y)
    print np.amin(array_y)
    print np.amax(array_y)
    print np.amin(array_x)
    print np.amax(array_x)
  
    print array_x.size
    print array_y.size
    print array_z.size

    col_list = []
    
    cmap = get_cmap(array_z.size)
    
    for x in range(array_z.size):
	col_list.append(cmap(float(x)))


     
    fm_ind=np.arange(array_z.size)
    print fm_ind
    

    
    
    

 
    savefig=False
    
    
    if ax is None:
        
        savefig=True
    
        fig = plt.figure(figsize=(thisfigsize[0], thisfigsize[0]*4./5))#figsize=thisfigsize
        grid = gridspec.GridSpec(1, 1, left=0.17, right=0.92, top=0.85, bottom=0.17,
                                 wspace=0.4, hspace=0.35)
        labeled_axes = []
        ax = plt.Subplot(fig, grid[0, 0])
        fig.add_subplot(ax)
        labeled_axes.append(ax)
    #color.cycle_cmap(rec_width[0,:].size, cmap='pink', ax=ax)
    
    
    ax.set_prop_cycle( cycler('color', col_list) )
    
     
    sig_ind=np.arange(array_col_m.size)
     
    cmx = plt.cm.get_cmap('plasma')
    #~ c = plt.cm.gist_rainbow(colors)
    #~ c = plt.cm.jey(colors)
    
    s_m = matplotlib.cm.ScalarMappable(cmap=cmx, norm=matplotlib.colors.Normalize(vmin=np.amin(fm_ind), vmax=np.amax(fm_ind)))
    s_m.set_array([])
   
   
    
    if "altperscolor" in xlabel or "altclustcolor" in xlabel:
      
        mincol=np.amin(array_z)
        maxcol=np.amax(array_z)
    
    print array_z
    
    print "color lims"
    print mincol, maxcol
 
 
    for i_m, m_mean in  enumerate(m_list   ):
        for i_f, f_mean in  enumerate(f_m_list   ):
            for i_p, ppl in enumerate(pplnum_list): 
                print i_p, ppl, markers[i_p]
            
                # ~ idxs=np.nonzero((array_cycle_F0==f_mean) & (array_cycle_ppl==ppl))[0]
                idxs=np.nonzero((array_cycle_F0==f_mean) & (array_cycle_ppl==ppl) & (array_cycle_M==m_mean))[0]
                
                array_z_now=array_z[idxs]
                X_now=X[idxs]
                Y_now=Y[idxs]  
                
                #~ if i_f==0:
                    #~ face="'none'"
                #~ else: 
                    #~ face=array_z
                #~ , edgecolor=''  
                
                #~ print i_f,   face      
                
                if show_legend:
                    legend=r'$F0$ {f_ms}, $N_h$ {ppl}'.format(f_ms=f_mean, ppl=fmt(ppl))
                else:
                    legend=''
                
        
            
                #~ 
                #~ sc = ax.scatter(X, Y, c=array_z, vmin=mincol, vmax=maxcol, marker='.', s=30, cmap=cmx, norm=matplotlib.colors.LogNorm())
                # ~ sc = ax.scatter(X_now, Y_now, c=array_z_now, vmin=mincol, vmax=maxcol, marker=markers[i_p], s=10, cmap=cmx, norm=matplotlib.colors.LogNorm(), label=legend)
                sc = ax.scatter(X_now, Y_now, c=array_z_now, vmin=mincol, vmax=maxcol, marker=markers[i_p], s=sizes[i_m], cmap=cmx, norm=matplotlib.colors.LogNorm(), label=legend)
                
                if i_f==1:
                    sc.set_facecolor('none')
     
        show_legend=False
    #~ ax.set_ylim(ymin=np.amin(array_y)/1.1, ymax=np.amax(array_y)*1.1)
    #~ ax.set_xlim(xmin=np.amin(array_x)/1.1, xmax=np.amax(array_x)*1.1)
    

    if show_cbar:
        
        #~ cbar = fig.colorbar(CS) # draw colorbar
        #cbar.set_label('{obs_name}'.format(obs_name=zlabel), rotation=270, labelpad=+8)
        
        if 'fig2' in xlabel:
            cbar = plt.colorbar(sc, ticks=[0.001, 0.01, 0.1]) 
            # ~ cbar = plt.colorbar(sc, ticks=np.unique(array_z)) 
            # ~ cbar.ax.set_yticklabels(['3*10^'])

        else:
            cbar = plt.colorbar(sc) 

        cbar.set_label(zlabel, rotation=270, labelpad=+8)
    

    if "altclustcolor" in xlabel:
        
        cbar = plt.colorbar(sc) 

        cbar.set_label(zlabel, rotation=270, labelpad=+8)
    

    if show_cbar_inset:
        
        
        axins = inset_axes(ax,
                   width="3%",  # width = 5% of parent_bbox width
                   height="50%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(0.08, 0.4, 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )

        # Controlling the placement of the inset axes is basically same as that
        # of the legend.  you may want to play with the borderpad value and
        # the bbox_to_anchor coordinate.

        # ~ fig.colorbar(im, cax=axins, ticks=[1, 2, 3])

        cbar = plt.colorbar(sc, cax=axins, ticks=[10**(-12), 10**(-9), 10**(-6)]) 
        
        cbar.ax.set_title(r'$D$')

        # ~ cbar.set_label(r'$D$', rotation=0, labelpad=+1)#
      
      
      
      
    if "pers time teo" in ylabel:
        ax.set_xlabel(r'numerical persistence time')
        ax.set_ylabel(r'theoretical persistence time')

    elif "var clust perp v3" in ylabel:
        ax.set_xlabel(r'$\sqrt{1.66} \sigma$')
        ax.set_ylabel(r'$\sigma_{\perp}$')

    elif "var clust perp v2" in ylabel:
        ax.set_xlabel(r'$\sqrt{2} \sigma_{\parallel}$')
        ax.set_ylabel(r'$\sigma_{\perp}$')
    elif "vel rescteo" in ylabel:
        ax.set_xlabel(r'numerical $\sigma ^2$')
        ax.set_ylabel(r'numerical $v/s$')
    elif "wave size teo" in ylabel:
        ax.set_xlabel(r'numerical $\sigma$')
        ax.set_ylabel(r'theoretical $\sigma$')
    elif "vel teo" in ylabel:
        ax.set_xlabel(r'numerical $v$')
        ax.set_ylabel(r'theoretical $v$')
    elif "fittest space teo" in ylabel:
        ax.set_xlabel(r'numerical $u_c$')
        ax.set_ylabel(r'theoretical $u_c$')
    elif "sel teo" in ylabel:
        ax.set_xlabel(r'numerical $s$')
        ax.set_ylabel(r'theoretical $s$')
    elif "vtau teo" in ylabel:
        ax.set_xlabel(r'numerical $v\tau$')
        ax.set_ylabel(r'theoretical $v\tau$')
    elif "num virs teo" in ylabel:
        ax.set_xlabel(r'numerical $N$')
        ax.set_ylabel(r'theoretical $N$')
    elif "D X extinction rate" in ylabel:
        ax.set_xlabel(r'numerical $v/N$')
        ax.set_ylabel(r'numerical (extinction rate) $ D$')
    elif "D13 X extinction rate" in ylabel:
        ax.set_xlabel(r'numerical $1/N$')
        ax.set_ylabel(r'numerical (extinction rate) $ D^{1/3}$')
    elif "1oN model" in xlabel:
        ax.set_xlabel(r'numerical $1/N$')
        ax.set_ylabel(r'numerical extinction rate')
    elif "vel resc teo clust" in xlabel:
        ax.set_xlabel(r'theoretical $v/\sqrt{D}$')
        ax.set_ylabel(r'rate of lineage splitting')
    elif "sel teo clust 4" in xlabel:
        ax.set_xlabel(r'theoretical $S$')
        ax.set_ylabel(r'rate of lineage splitting')
    elif "N model clust" in xlabel:
        ax.set_xlabel(r'numerical $N$')
        ax.set_ylabel(r'rate of lineage splitting')
    elif "I model clust 3" in xlabel:
        ax.set_xlabel(r'numerical $I$')
        ax.set_ylabel(r'rate of lineage splitting')
        
    elif "RATE SPEC exp"  in ylabel:
        ax.set_ylabel(r'$\kappa_{\rm split}/((3/2)^{1/2} v^2 / 8 D)$')
        ax.set_xlabel(r'$\mathcal{L}$')
        
        
    else:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
    ax.xaxis.labelpad = axis_labelpad
    ax.yaxis.labelpad = axis_labelpad
    #ax.legend( [sc], "fi", frameon=False, columnspacing=0.5, handletextpad=0.2,
    

    if "single lineage probability"  in ylabel:
        ax.set_xscale('log')
        ax.set_xlim(xmin=np.amin(X)/2., xmax=np.amax(X)*2.)

    elif "RATE SPEC exp"  in ylabel:
                
        
        newX = np.linspace(np.amin(X), np.amax(X), 100)
        popt, pcov = curve_fit(myExpFunc, X, Y)
        
        # ~ ax.loglog(X, X, marker='', linestyle='-', linewidth=0.5, color='grey')
        ax.plot(newX, myExpFunc(newX, *popt), marker='', linestyle='-', linewidth=0.5, color='grey')
        
        print newX
        print myExpFunc(newX, *popt)
        print popt
        
        # ~ sys.exit()
        
        # ~ plt.plot(newX, myExpFunc(newX, *popt), 'r-', 
                 # ~ label="({0:.3f}*x**{1:.3f})".format(*popt))
        
        
        ax.set_yscale('log')
        ax.set_ylim(ymin=np.amin(Y)/2., ymax=np.amax(Y)*2.)

    elif "rescaled pop size" not in xlabel:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.loglog(X, X, marker='', linestyle='-', linewidth=0.5, color='grey')
        
    
        ax.set_ylim(ymin=np.amin(Y)/2., ymax=np.amax(Y)*2.)
        ax.set_xlim(xmin=np.amin(X)/2., xmax=np.amax(X)*2.)
        #~ ax.set_xlim(xmax=0.1)
        #~ ax.set_ylim(ymax=0.1)
    
    elif "extinction rate" in ylabel:

        ax.set_yscale('log')

        ax.set_ylim(ymin=np.amin(Y)/2., ymax=np.amax(Y)*2.)
    
    else:        
        ax.set_xlim(xmin=0, xmax=20)
        ax.set_ylim(ymin=0, ymax=8)
        ax.set_xticks([2,6,10,14,18])
        ax.set_yticks([0,2,4,6,8])
            
        
    #~ custom_lines = [Line2D([0], [0], color=cmap(0.), lw=4),
                    #~ Line2D([0], [0], color=cmap(.5), lw=4),
                    #~ Line2D([0], [0], color=cmap(1.), lw=4)]
    
    #~ fig, ax = plt.subplots()
    #~ lines = ax.plot(data)
    #~ ax.legend( ['Cold', 'Medium', 'Hot'], bbox_to_anchor=(0.7, 1.1))    
     
    #~ plt.gca().legend(('y0','y1'), bbox_to_anchor=(0.7, 1.1))
    
    bbox_to_anchor=(0.95, 1.22)
    if not savefig:
        bbox_to_anchor=(1.5, 1.65)



            
    leg = ax.legend(frameon=False, ncol=2, columnspacing=0.5, handletextpad=0.2,
	  loc='upper right', bbox_to_anchor=bbox_to_anchor, prop={'size': 7})
    mpsetup.despine(ax) 
    
    #leg = ax.get_legend()
    for handler in leg.legendHandles:
        handler.set_color('black') 
    for handler in leg.legendHandles[-len(pplnum_list):]:
        handler.set_facecolor('none')
    
    #### finish figure ####
    labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
	     xycoords=('axes fraction'), fontweight = 'bold')
    #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
    #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
    
    xstring=xlabel.replace(" ", "_")
    ystring=ylabel.replace(" ", "_")
    title_string='{xstring}-{ystring}'.format(xstring=xstring, ystring=ystring)

    
    print title_string
    
    out_file='{out}/{title_string}.pdf'.format(out=dir_out_plots_scatter, title_string=title_string)
    #    print out_file, dpi=300
    
    if "RATE SPEC exp"  in ylabel:
        
        out_file='{out}/fig4.pdf'.format(out=dir_o_4)
        
    if savefig:

        fig.savefig(out_file)
        fig.clf()
        plt.close('all')
    
	  
      









    



def plot_scatter_allparams_4panels(array_x1, xlabel1, array_y1, ylabel1, array_x2, xlabel2, array_y2, ylabel2, array_x3, xlabel3, array_y3, ylabel3, array_x4, xlabel4, array_y4, ylabel4, title):


    if 'fig2' in xlabel1:
    
        width_ratios=[1, 1]
        
        fig_4pan = plt.figure(figsize=(thisfigsize[0]*4./3.,thisfigsize[0]*4./3))
        grid = gridspec.GridSpec(2, 2, left=0.13, right=0.96, top=0.96, bottom=0.1,
             wspace=0.35, hspace=0.3, width_ratios=width_ratios)


    else:
            
        width_ratios=[1, 5./4]
        
        fig_4pan = plt.figure(figsize=(thisfigsize[0]*4./3.,thisfigsize[1]*2.*4./3))
        grid = gridspec.GridSpec(2, 2, left=0.13, right=0.92, top=0.83, bottom=0.13,
             wspace=0.5, hspace=0.55, width_ratios=width_ratios)
   
    labeled_axes = []
    
   
    ax = plt.Subplot(fig_4pan, grid[0, 0]) 
    fig_4pan.add_subplot(ax)
    labeled_axes.append(ax)
    
 
 
    plot_scatter_colorbar(array_x1, array_y1, xlabel1, ylabel1, ax, show_cbar=False, fig=fig_4pan, show_legend=True)


   
    ax = plt.Subplot(fig_4pan, grid[0, 1]) 
    fig_4pan.add_subplot(ax)
    labeled_axes.append(ax)
 
 
    plot_scatter_colorbar(array_x2, array_y2, xlabel2, ylabel2, ax, show_cbar=True, fig=fig_4pan, show_legend=False)


       
   
    ax = plt.Subplot(fig_4pan, grid[1, 0]) 
    fig_4pan.add_subplot(ax)
    labeled_axes.append(ax)
 
 
    plot_scatter_colorbar(array_x3, array_y3, xlabel3, ylabel3, ax, show_cbar=False, fig=fig_4pan, show_legend=False)


   
    ax = plt.Subplot(fig_4pan, grid[1, 1]) 
    fig_4pan.add_subplot(ax)
    labeled_axes.append(ax)
 
 
    plot_scatter_colorbar(array_x4, array_y4, xlabel4, ylabel4, ax, show_cbar=True, fig=fig_4pan, show_legend=False)






    labeldict = dict(labelstyle=r'%s', fontsize='large',
		     xycoords=('axes fraction'), fontweight = 'bold')
    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.3,  0.97), **labeldict)
    mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3,  0.97), **labeldict)
    mpsetup.label_axes([labeled_axes[2]], labels='C', xy=(-0.3,  0.97), **labeldict)
    mpsetup.label_axes([labeled_axes[3]], labels='D', xy=(-0.3,  0.97), **labeldict)


    print title
    
    out_file='{out}/{title_string}.pdf'.format(out=dir_out_plots_scatter, title_string=title)
    if 'fig2' in title:
        out_file='{out}/{title_string}.pdf'.format(out=dir_out_plots_scatter, title_string=title)
    #    print out_file, dpi=300
    fig_4pan.savefig(out_file)
    fig_4pan.clf()
    plt.close('all')
    
	  




	  
	 

def vel_fromparams(recog_width, diff_const, mem_points, F_0, vir_number):
    return (diff_const**(2/3.))*((mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)**(1/3.))*((24*np.log((diff_const**(1/3.)) * ((mem_points*(np.exp(F_0/mem_points) -1.)/recog_width)**(2/3.)) * vir_number))**(1/3.)) 
    

def sig_fromparams(recog_width, diff_const, mem_points, F_0, vir_number):
    return (diff_const**(1/3.))*((mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)**(-1/3.))*((24*np.log((diff_const**(1/3.)) * ((mem_points*(np.exp(F_0/mem_points) -1.)/recog_width)**(2/3.)) * vir_number))**(1/6.)) 
    

def nose_fromparams(recog_width, diff_const, mem_points, F_0, vir_number):
    return (diff_const**(1/3.))*((mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)**(-1/3.))*((24*np.log((diff_const**(1/3.)) * ((mem_points*(np.exp(F_0/mem_points) -1.)/recog_width)**(2/3.)) * vir_number))**(2/3.))/4. 
    

def sel_fromparams(recog_width, mem_points, F_0):
    return (mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)
    

def vtauor_fromparams(recog_width, mem_points, F_0):
    return (1./(np.exp(F_0/mem_points)-1.))
    

def vtau_fromparams(recog_width, mem_points, F_0):
    return (recog_width/(np.exp(F_0/mem_points)-1.))
    

def tau_fromparams(mem_points, host_number, vir_number):
    return (mem_points*host_number)/vir_number.astype(float)
    
        
    

def R_persl_fromparams(recog_width, mem_points, F_0):
     
    # ~ return mem_points*recog_width/(sel_fromparams(recog_width, mem_points, F_0)*recog_width + mem_points) 
    return recog_width* np.exp(-F_0/mem_points)  
    
    



def Trascendental_N_iterative_fill(N, recog_width, diff_const, mem_points, F_0, ppl_num, count, res, i_n ):
    
    #~ return 0
    
    if np.isfinite(N):
        
        count+=1
        
        N_next=ppl_num*(diff_const**(2/3.))*((mem_points*(np.exp(F_0/mem_points)-1.)/recog_width)**(4/3.))*((24*np.log((diff_const**(1/3.)) * ((mem_points*(np.exp(F_0/mem_points) -1.)/recog_width)**(2/3.)) * N))**(1/3.))
        
        print count, N, N_next
        
        if np.abs( N - N_next) < 100:
            print"found root"
            res[i_n]  =  N_next
        
        else:
            Trascendental_N_iterative_fill(N_next, recog_width, diff_const, mem_points, F_0, ppl_num, count, res, i_n)
 

        


def Trascendental_N_iterative_vec(N_vec, recog_width_vec, diff_const_vec, mem_points_vec, F_0_vec, ppl_num_vec, count):
    
    #~ return 0
  
    res=np.zeros_like(N_vec.ravel()) - 1
    
    
    for i_n, N in enumerate(N_vec.ravel()):
        
        recog_width = recog_width_vec.ravel()[i_n] 
        diff_const          = diff_const_vec         .ravel()[i_n] 
        mem_points     = mem_points_vec    .ravel()[i_n] 
        F_0             = F_0_vec            .ravel()[i_n] 
        ppl_num        = ppl_num_vec       .ravel()[i_n]
    
        
        if np.isfinite(N):
            
            Trascendental_N_iterative_fill(N, recog_width, diff_const, mem_points, F_0, ppl_num, count, res, i_n )
            
        print i_n, res[i_n]
        print
        
    return res
   
   
    

        



#rom=['I','II','III','IV']
rom=['i','ii','iii','iv']

thisfigsize = figsize
thisfigsize[1] *= 0.75

title="pheno_space_ph_diag"
title_2pan="var_vs_vel"


# ~ jump_size=2.
# ~ cycles_per_year=100.
cycles_per_year=1. # time units
 
    
dir_io_orig='../../fig2' # directory with input files
dir_o_4='../../fig4' # directory with input files
    
dir_out_plots_tot='{inp}/'.format(inp=dir_io_orig)

if not os.path.exists(dir_o_4):
    os.makedirs(dir_o_4)	


data_arr=np.empty(0)

for dir_io_tot  in ['{inp}/from_zuzia_coarse_grained'.format(inp=dir_io_orig), '{inp}/from_vision_coarse_grained'.format(inp=dir_io_orig)   ] :# directory with input files
 
    
    if 'zuzia' in dir_io_tot:
        m_list=[1.,5.]
        f_m_list=[1.,3.]
        #~ f_m_list=[1.]
        pplnum_list=[100000000,10000000000, 1000000000000]
        #~ pplnum_list=[1000000000000]
        width_ratios=[1, 1, 1, 6./5]
        left_sp=0.1
    
    if 'vision' in dir_io_tot:
        m_list=[1., 5.]
        f_m_list=[1., 3.]
        pplnum_list=[100000000,10000000000, 1000000000000]
        #~ pplnum_list=[1000000000000]
        width_ratios=[1, 1, 1, 6./5]
        left_sp=0.1
        
    numcol=len(f_m_list)
    
        
    for i_f, f_mean in enumerate(f_m_list):
        for i_p, ppl in enumerate(pplnum_list):
    
            f_m_string='{f_m:.{prec}f}'.format(f_m=f_mean, prec=1)
            
            f_m_string = f_m_string.replace(".", "d")[:-1]
            
            dir_io='{inp}/summary_analysis_{f_m_string}_{ppl}'.format(inp=dir_io_tot, f_m_string=f_m_string, ppl=ppl)

            in_file= '{inp}/summary_analysis.dat'.format(inp=dir_io)
            
            #in_file=sys.argv[2]
    
    
    
            data_arr_spec = np.loadtxt(in_file)
            
            if data_arr.size==0:
                data_arr=data_arr_spec.copy()
            else:
                data_arr=np.concatenate((data_arr, data_arr_spec), axis=0)
            
            print data_arr_spec.shape
            print data_arr.shape
                
                
            data_out=data_arr_spec[data_arr_spec[:,3] >= 0.01,:]
            
            file_out= '{inp}/summary_analysis.dat'.format(inp=dir_io)
            
            with open(file_out,'w') as f_handle:
                np.savetxt(f_handle, data_out, fmt='%15.15f')
                
                

print data_arr.shape



#~ data_arr = data_arr[(avg_frac_inf_reals_mean>0) & (cluster_counts>0) & (np.isfinite(coalescent)),:]


print data_arr[0,:]
print data_arr[:,0]


print data_arr[0,0]
print data_arr[0,2]
print data_arr[0,3]

print data_arr[1,2]
print data_arr[1,3]

ppl_num=int(data_arr[0,0])
ppl_num_param=int(data_arr[0,0])


dir_out_plots_scatter='{inp}/'.format(inp=dir_out_plots_tot)


    
    
# ~ show_colormap('magma')
    
# eliminate points that terminate early too often, and that don't have "early explosions", meaning eplosion prob not too high or valid sim time not too low


#data_arr=data_arr[(end_bef_init_prob_all<0.8) & ((explosion_prob_all<0.7) | (valid_simulation_time_all>60)),:]

fails                                       = data_arr[:,6]
fails           = np.where(fails ==-1, 1, 0)

print fails


data_arr_all=data_arr.copy()
data_arr=data_arr[(fails==0),:]

print data_arr.shape



all_real_rec_width                                     = data_arr_all[:,4] 

all_mem_points                                    = data_arr_all[:,0]
all_F0                                            = data_arr_all[:,1]
all_ppl_num                                       = data_arr_all[:,2]
all_mu                                            = data_arr_all[:,3]
all_time_tot_run                                  = data_arr_all[:,5]/cycles_per_year





real_rec_width                                     = data_arr[:,4] 


mem_points                                    = data_arr[:,0]
F0                                            = data_arr[:,1]
ppl_num                                       = data_arr[:,2]
mu                                            = data_arr[:,3]
mu                                            = np.where(mu ==0, 10**-6, mu)

rec_width                                     = 1.
time_tot_run                                  = data_arr[:,5]/cycles_per_year
probext                                       = data_arr[:,6]
probexpl                                      = data_arr[:,7]
numext                                        = data_arr[:,8]
numexpl                                       = data_arr[:,9]
avg_I                                         = data_arr[:,10]
std_I                                         = data_arr[:,11]
var_log_I                                     = data_arr[:,12]
avg_fit_avg                                   = data_arr[:,13]*(cycles_per_year)
var_fit_avg                                   = data_arr[:,14]*(cycles_per_year**2)
avg_fit_var                                   = data_arr[:,15]*(cycles_per_year**2)
var_fit_var                                   = data_arr[:,16]*(cycles_per_year**4)
avg_num_IS_coords                             = data_arr[:,17]
avg_num_vir_coords                            = data_arr[:,18]
avg_mean_displ_x                              = data_arr[:,19]# keep simulation units
avg_mean_displ_y                              = data_arr[:,20]
avg_var_displ_x                               = data_arr[:,21]
avg_var_displ_y                               = data_arr[:,22]
avg_mean_displ_tot                            = data_arr[:,23]
#~ diff_maxfit_fit_rw_fct_time                   = data_arr[:,24]
#~ err_diff_maxfit_fit_rw_fct_time               = data_arr[:,25]
#~ chisq_maxfit_fit_rw_fct_time                  = data_arr[:,26]
avg_sig_est                                   = data_arr[:,27]/(real_rec_width)
var_sig_est                                   = data_arr[:,28]/(real_rec_width**2)

avg_time_switch                               = data_arr[:,29]# keep sim units
avg_time_fullconv                             = data_arr[:,30]
max_fitn_abserr_switch_max_rel                = data_arr[:,31]
avg_fitn_abserr_switch_max_rel                = data_arr[:,32]


#~ frac_1cl                                      = data_arr[:,33]
#~ avg_size_clusters_split                       = data_arr[:,34]
avg_n_clusters_split                          = data_arr[:,35]
#~ avg_num_vir_clusters_split                    = data_arr[:,36]
#~ avg_num_coords_clusters_split                 = data_arr[:,37]
#~ avg_num_ext_tot                               = data_arr[:,38]
#~ avg_num_split_tot                             = data_arr[:,39]
#~ avg_num_ext_per_cl_tot                        = data_arr[:,40]
#~ avg_num_split_per_cl_tot                      = data_arr[:,41]
#~ avg_num_vir_o_nclust_tot_list                 = data_arr[:,42]
#~ avg_vel_clust_split_alltimes_allclust         = data_arr[:,43]
#~ avg_var_parall_clust_split_alltimes_allclust  = data_arr[:,44]
#~ avg_var_perp_clust_split_alltimes_allclust    = data_arr[:,45]
#~ avg_var_tot_clust_split_alltimes_allclust     = data_arr[:,46]
#~ avg_inter_clusters_dist_split                 = data_arr[:,47]
#~ max_inter_clusters_dist_split                 = data_arr[:,48]
#~ max_inter_clusters_dist_max_split             = data_arr[:,49]
#~ var_vel_clust_split_alltimes_allclust         = data_arr[:,50]
#~ var_var_parall_clust_split_alltimes_allclust  = data_arr[:,51]
#~ var_var_perp_clust_split_alltimes_allclust    = data_arr[:,52]
#~ var_var_tot_clust_split_alltimes_allclust     = data_arr[:,53]
#~ rate_num_ext_tot                              = data_arr[:,54]
#~ rate_num_split_tot                            = data_arr[:,55]
#~ rate_num_ext_per_cl_tot                       = data_arr[:,56]
#~ rate_num_split_per_cl_tot                     = data_arr[:,57]
#~ num_sharp_turns                               = data_arr[:,58]
#~ rate_sharp_turns                              = data_arr[:,59]
#~ avg_shortmem_err_ys                           = data_arr[:,60]
#~ avg_shortmem_err_svsx                         = data_arr[:,61]
#~ frac_shortmem_err_ys                          = data_arr[:,62]
#~ frac_shortmem_err_svsx                        = data_arr[:,63]
#~ diff_maxfit_fit_rw_fct_time                   = data_arr[:,64]
#~ err_diff_maxfit_fit_rw_fct_time               = data_arr[:,65]
#~ chisq_maxfit_fit_rw_fct_time                  = data_arr[:,66]
timelast                                      = data_arr[:,67]
#~ NUM_VIR                                       = data_arr[:,68]
#~ TAU                                           = data_arr[:,69]
#~ VEL                                           = data_arr[:,70]
#~ S                                             = data_arr[:,71]
#~ SIG_SPACE                                     = data_arr[:,72]
#~ NOSE_SPACE                                    = data_arr[:,73]
#~ R_PERS_L                                      = data_arr[:,74]
avg_S_est                                     = data_arr[:,77]*(real_rec_width*cycles_per_year)

#~ pers_l_avg_end2end_fit                               = data_arr[:,78]
#~ err_pers_l_avg_end2end_fit                           = data_arr[:,79]
chisq_avg_end2end_fit                                = data_arr[:,80]
#~ pers_l_avg_end2end_fit_varvel_fixvel                 = data_arr[:,81]
#~ beta_avg_end2end_fit_varvel_fixvel                   = data_arr[:,82]
#~ err_pers_l_avg_end2end_fit_varvel_fixvel             = data_arr[:,83]
#~ err_beta_avg_end2end_fit_varvel_fixvel               = data_arr[:,84]
#~ chisq_avg_end2end_fit_varvel_fixvel                  = data_arr[:,85]
#~ diff_avg_angles_subseg_fit_rw_fct_time               = data_arr[:,86]
#~ err_diff_avg_angles_subseg_fit_rw_fct_time           = data_arr[:,87]
#~ chisq_avg_angles_subseg_fit_rw_fct_time              = data_arr[:,88]
#~ pers_l_avg_end2end_fit_constvel_fct_time_fixvel      = data_arr[:,89]
#~ err_pers_l_avg_end2end_fit_constvel_fct_time_fixvel  = data_arr[:,90]
#~ chisq_avg_end2end_fit_constvel_fct_time_fixvel       = data_arr[:,91]
#~ pers_l_avg_end2end_fit_varvel_fct_time_fixvel        = data_arr[:,92]
#~ beta_avg_end2end_fit_varvel_fct_time_fixvel          = data_arr[:,93]
#~ err_pers_l_avg_end2end_fit_varvel_fct_time_fixvel    = data_arr[:,94]
#~ err_beta_avg_end2end_fit_varvel_fct_time_fixvel      = data_arr[:,95]
#~ chisq_avg_end2end_fit_varvel_fct_time_fixvel         = data_arr[:,96]










D=mu*(real_rec_width**2.)
D_mod=mu*(avg_var_displ_x/2.) # simulations units!
# ~ D_mod=D # simulations units!



count=0
#~ res_N= np.asarray([Trascendental_N_iterative(n, real_rec_width.ravel()[i_n], D_mod.ravel()[i_n], mem_points.ravel()[i_n], F0.ravel()[i_n], ppl_num.ravel()[i_n], count) for i_n, n in enumerate((avg_I*ppl_num).ravel())])
res_N= Trascendental_N_iterative_vec(avg_I*ppl_num, real_rec_width, D_mod, mem_points, F0, ppl_num, count)

#~ res_N=np.reshape(res_N, (-1, ppl_num.shape[0]))



VEL          = vel_fromparams(real_rec_width, D_mod, mem_points, F0, res_N)/(real_rec_width/cycles_per_year)
VEL_Nteo_all          = vel_fromparams(real_rec_width, D_mod, mem_points, F0, res_N)/(real_rec_width/cycles_per_year)
VEL_Nteo          = vel_fromparams(real_rec_width, D_mod, mem_points, F0, res_N)/(real_rec_width/cycles_per_year)
#~ VEL          = vel_fromparams(real_rec_width, D_mod, mem_points, F0, res_N)/(real_rec_width/cycles_per_year)
SIG_SPACE    = sig_fromparams(real_rec_width, D_mod, mem_points, F0, res_N)/(real_rec_width)
NOSE_SPACE   = nose_fromparams(real_rec_width, D_mod, mem_points, F0, res_N)/(real_rec_width)
S            = sel_fromparams(real_rec_width, mem_points, F0)*(real_rec_width*cycles_per_year)  # units 1/(space*time)
R_PERS_L     = R_persl_fromparams(real_rec_width, mem_points, F0)/(real_rec_width) # units space

VTAUOR     = vtauor_fromparams(real_rec_width, mem_points, F0) # units space
TAUOR        = (tau_fromparams(mem_points, ppl_num, res_N)/cycles_per_year) # units space
TAU        = (tau_fromparams(mem_points, ppl_num, res_N)/cycles_per_year) # units space
TAU_Nteo        = (tau_fromparams(mem_points, ppl_num, res_N)/cycles_per_year) # units space
VTAU    = vtau_fromparams(real_rec_width, mem_points, F0)/(real_rec_width) # units space

avg_S_est_TOT=avg_S_est.copy()
S_TOT=S.copy()



D_mod_all=D_mod.copy() # now right units
D_mod=D_mod/(real_rec_width**2./cycles_per_year) # now right units



#~ Trasc_N= Trascendental_N_params(real_rec_width, D_mod, mem_points, F0, ppl_num)

#~ mask=np.isfinite(np.log((diff_const**(1/3.)) * ((mem_points*(np.exp(F_0/mem_points) -1.)/recog_width)**(2/3.)) * N))
#~ N=N[mask]
#~ N=N[mask]
#~ N=N[mask]
        

#~ res_N = root(Trasc_N, 0.5).x[0] 

print res_N
print res_N.shape, ppl_num.shape



#~ print VEL         
#~ print SIG_SPACE   
#~ print NOSE_SPACE  
#~ print S           
#~ print R_PERS_L    
#~ print (D_mod**(1/3.)) * ((mem_points*(np.exp(F0/mem_points) -1.)/real_rec_width)**(2/3.)) *  avg_I*ppl_num    
#~ print  avg_I*ppl_num    

#~ sys.exit()

thisfigsize = figsize
thisfigsize[1] *= 0.75


print time_tot_run




fails_clust           = np.where(avg_n_clusters_split ==0, 1, 0)


data_arr_clust=data_arr[avg_n_clusters_split >0,:]

frac_1cl                                      = data_arr_clust[:,33]
avg_size_clusters_split                       = data_arr_clust[:,34]
avg_n_clusters_split                          = data_arr_clust[:,35]
avg_num_vir_clusters_split                    = data_arr_clust[:,36]
avg_num_coords_clusters_split                 = data_arr_clust[:,37]
#~ avg_num_ext_tot                               = data_arr_clust[:,38]
#~ avg_num_split_tot                             = data_arr_clust[:,39]
#~ avg_num_ext_per_cl_tot                        = data_arr_clust[:,40]
#~ avg_num_split_per_cl_tot                      = data_arr_clust[:,41]
avg_num_vir_o_nclust_tot_list                 = data_arr_clust[:,42]
avg_vel_clust_split_alltimes_allclust         = data_arr_clust[:,43]/(real_rec_width [~fails_clust.astype(bool)]/cycles_per_year)
avg_var_parall_clust_split_alltimes_allclust  = data_arr_clust[:,44]/(real_rec_width [~fails_clust.astype(bool)]**2)
avg_var_perp_clust_split_alltimes_allclust    = data_arr_clust[:,45]/(real_rec_width [~fails_clust.astype(bool)]**2)
avg_var_tot_clust_split_alltimes_allclust     = data_arr_clust[:,46]/(real_rec_width [~fails_clust.astype(bool)]**2)
avg_inter_clusters_dist_split                 = data_arr_clust[:,47]/(real_rec_width [~fails_clust.astype(bool)])
max_inter_clusters_dist_split                 = data_arr_clust[:,48]/(real_rec_width [~fails_clust.astype(bool)])
max_inter_clusters_dist_max_split             = data_arr_clust[:,49]/(real_rec_width [~fails_clust.astype(bool)])
var_vel_clust_split_alltimes_allclust         = data_arr_clust[:,50]/((real_rec_width [~fails_clust.astype(bool)]/cycles_per_year)**2)
var_var_parall_clust_split_alltimes_allclust  = data_arr_clust[:,51]/(real_rec_width [~fails_clust.astype(bool)]**4)
var_var_perp_clust_split_alltimes_allclust    = data_arr_clust[:,52]/(real_rec_width [~fails_clust.astype(bool)]**4)
var_var_tot_clust_split_alltimes_allclust     = data_arr_clust[:,53]/(real_rec_width [~fails_clust.astype(bool)]**4)
rate_num_ext_tot                              = data_arr_clust[:,54]*cycles_per_year
rate_num_split_tot                            = data_arr_clust[:,55]*cycles_per_year
rate_num_ext_per_cl_tot                       = data_arr_clust[:,56]*cycles_per_year
rate_num_split_per_cl_tot                     = data_arr_clust[:,57]*cycles_per_year
#~ num_sharp_turns                            = data_arr_clust[:,58]
#~ rate_sharp_turns                           = data_arr_clust[:,59]
#~ avg_shortmem_err_ys                        = data_arr_clust[:,60]
#~ avg_shortmem_err_svsx                      = data_arr_clust[:,61]
#~ frac_shortmem_err_ys                       = data_arr_clust[:,62]
#~ frac_shortmem_err_svsx                     = data_arr_clust[:,63]
#~ diff_maxfit_fit_rw_fct_time                = data_arr_clust[:,64]
#~ err_diff_maxfit_fit_rw_fct_time            = data_arr_clust[:,65]
#~ chisq_maxfit_fit_rw_fct_time               = data_arr_clust[:,66]


avg_grad_fittest_split_alltimes_allclust      = data_arr_clust[:,75]*(real_rec_width [~fails_clust.astype(bool)]*cycles_per_year)
avg_shortmem_gradx_bulk_split_alltimes_allclust      = data_arr_clust[:,76]*(real_rec_width [~fails_clust.astype(bool)]*cycles_per_year)
#~ VEL                                           = data_arr_clust[:,70]
#~ S                                             = data_arr_clust[:,71]
#~ SIG_SPACE                                     = data_arr_clust[:,72]
#~ NOSE_SPACE                                    = data_arr_clust[:,73]
avg_S_est                                     = data_arr_clust[:,77]*(real_rec_width [~fails_clust.astype(bool)]*cycles_per_year)






rate_num_split_tot_0d1            = data_arr_clust[:,97]*cycles_per_year
rate_num_split_tot_0d5            = data_arr_clust[:,98]*cycles_per_year
rate_num_split_tot_1              =data_arr_clust[:,99]*cycles_per_year
rate_num_split_tot_10             =data_arr_clust[:,100]*cycles_per_year
rate_num_split_per_cl_tot_0d1     =data_arr_clust[:,101]*cycles_per_year
rate_num_split_per_cl_tot_0d5     =data_arr_clust[:,102]*cycles_per_year
rate_num_split_per_cl_tot_1       =data_arr_clust[:,103]*cycles_per_year
rate_num_split_per_cl_tot_10      =data_arr_clust[:,104]*cycles_per_year








#~ D_mod=D_mod[~fails_clust.astype(bool)]
#~ D=D[~fails_clust.astype(bool)]

D_mod          = D_mod         [~fails_clust.astype(bool)]
VEL          = VEL         [~fails_clust.astype(bool)]
VEL_Nteo          = VEL_Nteo_all         [~fails_clust.astype(bool)]
SIG_SPACE    = SIG_SPACE   [~fails_clust.astype(bool)]
# ~ NOSE_SPACE   = NOSE_SPACE  [~fails_clust.astype(bool)]
S            = S           [~fails_clust.astype(bool)]
VTAUOR            = VTAUOR           [~fails_clust.astype(bool)]
TAUOR               = TAUOR             [~fails_clust.astype(bool)]
TAU                = TAU              [~fails_clust.astype(bool)]
TAU_Nteo                = TAU_Nteo              [~fails_clust.astype(bool)]
VTAU               = VTAU             [~fails_clust.astype(bool)]
#~ R_PERS_L     = R_PERS_L    [~fails_clust.astype(bool)]

timescale_system = R_PERS_L[~fails_clust.astype(bool)]/VEL_Nteo_all[~fails_clust.astype(bool)]



plot_scatter_allparams_4panels(avg_I*ppl_num, "num virs model fig2", res_N, "num virs teo", avg_vel_clust_split_alltimes_allclust, "vel clust model fig2", VEL, "vel teo", np.sqrt(avg_var_parall_clust_split_alltimes_allclust),  "wave size clust model fig2", SIG_SPACE, "wave size teo", np.sqrt(1.66)*np.sqrt(avg_var_parall_clust_split_alltimes_allclust), "var clust paral fig2", np.sqrt(avg_var_perp_clust_split_alltimes_allclust), "var clust perp v3", "fig2")
# ~ sys.exit()


RATE_SPEC_exp=((2.**9) * (S**3.) * (D_mod**2.)* (0.1**2)* (R_PERS_L[~fails_clust.astype(bool)]**2.)/(9* (VEL_Nteo**5.)))**(1./4)

plot_scatter_colorbar(RATE_SPEC_exp,rate_num_split_per_cl_tot_0d1/(np.sqrt(3./8)*((VEL_Nteo**2.)/(4*D_mod))), "rate of speciation per clust 0d1 resc fig2",  "RATE SPEC exp")#*( avg_I [~fails_clust.astype(bool)]*ppl_num [~fails_clust.astype(bool)])**(0.)


data_out = np.array([mu[~fails_clust.astype(bool)], ppl_num[~fails_clust.astype(bool)], F0[~fails_clust.astype(bool)], mem_points[~fails_clust.astype(bool)], real_rec_width[~fails_clust.astype(bool)],  D_mod, S, res_N[~fails_clust.astype(bool)], timescale_system, TAU_Nteo, SIG_SPACE, NOSE_SPACE[~fails_clust.astype(bool)],  VEL_Nteo , VTAU, avg_I[~fails_clust.astype(bool)], rate_num_split_per_cl_tot ])

data_out=data_out.T

print data_out.shape

file_out='{inp}/splitting_rate.txt'.format(inp=dir_io_orig)


with open(file_out,'w') as f_handle:
    f_handle.write("# 1 mu  2 N_h   3 F0    4 M     5 r     6 D (units t_i/r**2)    7 S (units r/t_i)   8 Theoretical N     9 T     10 theoretical tau      11 theoretical sigma    12 theoretical u_c      13 theoretical v    14 theoretical vtau     15 numerical I      16 numerical lineage splitting rate")
    f_handle.write("\n")
    
with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_out, fmt='%15.15f')#, newline=" "
    # ~ f_handle.write("\n")
       


data_out = np.array([mu[~fails_clust.astype(bool)], ppl_num[~fails_clust.astype(bool)], F0[~fails_clust.astype(bool)], mem_points[~fails_clust.astype(bool)], real_rec_width[~fails_clust.astype(bool)],  D_mod, S, res_N[~fails_clust.astype(bool)], timescale_system, TAU_Nteo, SIG_SPACE, NOSE_SPACE[~fails_clust.astype(bool)],  VEL_Nteo , VTAU, avg_I[~fails_clust.astype(bool)], rate_num_split_tot_0d1, rate_num_split_per_cl_tot_0d1,  rate_num_split_tot_1, rate_num_split_per_cl_tot_1, numexpl[~fails_clust.astype(bool)]/time_tot_run[~fails_clust.astype(bool)] ])

data_out=data_out.T              
print data_out.shape      

data_out=data_out[data_out[:,0] >= 0.01,:]
print data_out.shape

file_out='{inp}/speciation_rate.txt'.format(inp=dir_io_orig)

with open(file_out,'w') as f_handle:
    f_handle.write("# 1 mu  2 N_h   3 F0    4 M     5 r     6 D (units t_i/r**2)    7 S (units r/t_i)   8 Theoretical N     9 T     10 theoretical tau      11 theoretical sigma    12 theoretical u_c      13 theoretical v    14 theoretical vtau     15 numerical I      16 numerical rate speciation threshold 0.1r   17   numerical rate speciation per lineage threshold 0.1r     18   numerical rate speciation threshold r     19   numerical rate speciation per lineage threshold r      20   explosion rate ")
    f_handle.write("\n")
    
with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_out, fmt='%15.15f')#, newline=" "
    # ~ f_handle.write("\n")
       
       
