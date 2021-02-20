

"""
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
#~ matplotlib.use('Agg')
#~ matplotlib.use('webagg')
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
import shutil
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import norm

col_dict=colors.cnames

#### subfigure A ####
## import data




def clusters_dist_vars(X, labels, clust, clust_dir): # returns three floats: the variance of distances from the centroid , parrallel to the direction of motion, perpendicular, and total
    #OLD
	print "cluster rotation"


	points_in_clust=X[labels==clust,:]
	
	points_in_clust=points_in_clust - np.mean(points_in_clust, axis=0)[None,:] # regularize
	
	print points_in_clust.shape
	print np.mean(points_in_clust, axis=0)


	clust_perpdir=np.array([clust_dir[1], -clust_dir[0]])
	
	Q = np.array([clust_dir, clust_perpdir]).T
	print 'Transformed to original system with \n Q={}'.format(Q)
	print 'Orthogonality check \n {}'.format(Q.dot(Q.T))
	print 'Orthogonality check \n {}'.format(np.dot(Q.T, X[0,:]) == np.linalg.solve(Q, X[0,:]))
	
	print np.dot(Q.T, points_in_clust[0,:])
	print np.linalg.solve(Q, points_in_clust[0,:])
	
	points_in_clust_newbasis=np.dot(Q.T, points_in_clust.T).T
	
	print "transformed direction of motion ", np.dot(Q.T, clust_dir)
	print "transformed direction perpendicular to motion ", np.dot(Q.T, clust_perpdir)
	
	print points_in_clust_newbasis.shape, np.mean(points_in_clust_newbasis, axis=0)
	
	# v_ = np.linalg.solve(Q, v)
	# print('The vector in the new coordinates \n v_={}'.format(v_))
	

	#dists_sq=np.linalg.norm(points_in_clust_newbasis - np.mean(points_in_clust_newbasis, axis=0)[None,:], axis=1)**2
	dists_sq=np.linalg.norm(points_in_clust_newbasis , axis=1)**2  
	print clust, points_in_clust.shape, dists_sq.shape
	var_dist=np.mean(dists_sq)
	print var_dist
	
	dists_sq_x=points_in_clust_newbasis[:,0]**2  
	print clust, dists_sq_x.shape
	var_dist_x=np.mean(dists_sq_x) # direction of motion
	print var_dist_x
	
	dists_sq_y=points_in_clust_newbasis[:,1]**2  
	print clust, dists_sq_y.shape
	var_dist_y=np.mean(dists_sq_y) # direction of motion
	print var_dist_y
	
	
	return [var_dist_x, var_dist_y, var_dist]



def clusters_dist_vars_all(X, clust_dir): # returns three floats: the variance of distances from the centroid , parrallel to the direction of motion, perpendicular, and total
    
	print "cluster rotation"


	points_in_clust=X[:,:2]
    
    
	print np.mean(points_in_clust, axis=0)
	print np.average(points_in_clust, weights=X[:,3], axis=0)
	
	#points_in_clust=points_in_clust - np.mean(points_in_clust, axis=0)[None,:] # regularize
	points_in_clust=points_in_clust - np.average(points_in_clust, weights=X[:,3], axis=0)[None,:] # regularize

	print points_in_clust.shape
	print np.average(points_in_clust, weights=X[:,3], axis=0)


	clust_perpdir=np.array([clust_dir[1], -clust_dir[0]])
	
	Q = np.squeeze(np.array([clust_dir, clust_perpdir])).T
	print 'Transformed to original system with \n Q={}'.format(Q)
	print 'Orthogonality check \n {}'.format(Q.dot(Q.T))
	print 'Orthogonality check \n {}'.format(np.dot(Q.T, X[0,:2]) == np.linalg.solve(Q, X[0,:2]))
	
	print np.dot(Q.T, points_in_clust[0,:])
	print np.linalg.solve(Q, points_in_clust[0,:])
	
	points_in_clust_newbasis=np.dot(Q.T, points_in_clust.T).T
	
	print "transformed direction of motion ", np.dot(Q.T, clust_dir)
	print "transformed direction perpendicular to motion ", np.dot(Q.T, clust_perpdir)
	
	print points_in_clust_newbasis.shape, np.mean(points_in_clust_newbasis, axis=0)
	
	# v_ = np.linalg.solve(Q, v)
	# print('The vector in the new coordinates \n v_={}'.format(v_))
	

	#dists_sq=np.linalg.norm(points_in_clust_newbasis - np.mean(points_in_clust_newbasis, axis=0)[None,:], axis=1)**2
	dists_sq=np.linalg.norm(points_in_clust_newbasis , axis=1)**2  
	print points_in_clust.shape, dists_sq.shape
	#var_dist=np.mean(dists_sq)
	var_dist=np.average(dists_sq, weights=X[:,3])
	print var_dist
	
	dists_sq_x=points_in_clust_newbasis[:,0]**2  
	print dists_sq_x.shape
	var_dist_x=np.average(dists_sq_x, weights=X[:,3])
	print var_dist_x
	
	dists_sq_y=points_in_clust_newbasis[:,1]**2  
	print dists_sq_y.shape
	var_dist_y=np.average(dists_sq_y, weights=X[:,3]) # direction of motion
	print var_dist_y
	
	
	return [var_dist_x, var_dist_y, var_dist]



dir_io='../../fig1_bcd' # directory with input files

   
def get_cmap(N):
    ''' Returns a function that maps each index in 0, 1, ...
        N-1 to a distinct RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    #scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap='jet')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

    

class MidpointNormalize_symlog(matplotlib.colors.SymLogNorm):

    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False, linthresh=None, linscale=None):
        self.vcenter = vcenter
        colors.SymLogNorm.__init__(self, vmin, vmax, clip, linthresh, linscale)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))



def vtau_fromparams(recog_width, mem_points, F_0):
    return (recog_width/(np.exp(F_0/mem_points)-1.))        



rec_width=2000
vtau=vtau_fromparams(rec_width, 1, 3.)

latt_sp=1.




thisfigsize = figsize
thisfigsize[1] *= 0.75


    
    

file_in_frames_npz_compr = "{inp}/IS_update_antigenic_space_time_49999.dat".format(inp=dir_io)
filename, _ = os.path.splitext(file_in_frames_npz_compr)
time=os.path.basename(filename).split('_')[-1]
time=int(time)

data_space = np.loadtxt(file_in_frames_npz_compr)

#print data_space
print data_space.shape
print data_space.ndim


if data_space.ndim == 1:     
    vir_x   = data_space[0] # 
    vir_y   = data_space[1] # 
    num_vir = data_space[3] # 
else:    
    vir_x = data_space[:,0] # 
    vir_y = data_space[:,1] # 
    num_vir = data_space[:,3]

fitn = data_space[:,4] # 
num_IS = data_space[:,2] # 
num_IS_upd = data_space[:,5] # 


print num_vir.shape

    
xmin_zoom=np.amin(vir_x)
xmax_zoom=np.amax(vir_x)
ymin_zoom=np.amin(vir_y)
ymax_zoom=np.amax(vir_y)

updmin=np.amin(num_IS_upd)
updmax=np.amax(num_IS_upd)




nx_zoom          = int((xmax_zoom - xmin_zoom)/(latt_sp)) +1
ny_zoom          = int((ymax_zoom - ymin_zoom)/(latt_sp)) +1  
x_zoom           = np.linspace(xmin_zoom, xmax_zoom, nx_zoom)
y_zoom           = np.linspace(ymin_zoom, ymax_zoom, ny_zoom)
xv_zoom, yv_zoom = np.meshgrid(x_zoom, y_zoom)    
    

Z_zoom = np.zeros(xv_zoom.shape, dtype=float) #+ 0.1
upd_zoom = np.zeros(xv_zoom.shape, dtype=float)  #- 1000
vir_zoom = np.zeros(xv_zoom.shape, dtype=float)  #- 1000


print Z_zoom.shape


for i_x, x_f in enumerate(vir_x):
    
    y_f   = vir_y[i_x]
    num_f = num_IS[i_x]
    upd_f = num_IS_upd[i_x]
    num_v = num_vir[i_x]
    
    
    
    idx_x_zoom = int((x_f - xmin_zoom)/latt_sp)
    idx_y_zoom = int((y_f - ymin_zoom)/latt_sp)
    
    

    
    
    #Z[idx_y, idx_x] = num_f
    Z_zoom[idx_y_zoom, idx_x_zoom] = num_f
    
    upd_zoom[idx_y_zoom, idx_x_zoom] += - upd_f
    vir_zoom[idx_y_zoom, idx_x_zoom] = num_v


max_maxcol=max([np.abs(updmin), np.abs(updmax)])
min_maxcol=min([np.abs(updmin), np.abs(updmax)])/1.

nummax=np.amax(num_IS)
numvirmax=np.amax(num_vir)
            
            
normalize_vir = matplotlib.colors.LogNorm(vmin=1, vmax=numvirmax)
normalize_is = matplotlib.colors.LogNorm(vmin=1, vmax=nummax)
#~ normalize_upd = matplotlib.colors.SymLogNorm(linthresh=1., linscale=1.,  vmin=updmin, vmax=updmax)
normalize_upd = matplotlib.colors.SymLogNorm(linthresh=1., linscale=1.,  vmin=-min_maxcol, vmax=min_maxcol)
#~ normalize_upd = matplotlib.colors.DivergingNorm(vmin=-500., vcenter=0, vmax=4000)
#~ normalize_upd = MidpointNormalize_symlog(linthresh=1., linscale=1., vmin=updmin, vcenter=0, vmax=updmax)
#~ normalize_upd = MidpointNormalize(vmin=updmin, vcenter=0, vmax=updmax)



fig_ins = plt.figure(figsize=(thisfigsize[0]*4./3, thisfigsize[0]))
grid_ins = gridspec.GridSpec(2, 2, left=0.05, right=0.93, top=0.95, bottom=0.07,
     wspace=0.35, hspace=0.4)
labeled_axes = []
ax_ins = plt.Subplot(fig_ins, grid_ins[1, 0])#, rasterized=True
fig_ins.add_subplot(ax_ins)
labeled_axes.append(ax_ins)
  
     
#ax_ins.plot(vir_x, vir_y, linestyle='', marker='x', color='r', markersize=0.5, label='viruses', rasterized=True)
#heatmap_fitn = ax_ins.pcolormesh(xv_zoom, yv_zoom, fitn_zoom,  cmap='Greys') # faster than pcolor , alpha=.9
    
#~ heatmap = ax_ins.pcolormesh(xv_zoom, yv_zoom, upd_zoom, cmap='RdBu', norm=normalize_upd) # faster than pcolor
ax_ins.set_rasterization_zorder(1)
heatmap = ax_ins.pcolormesh(xv_zoom, yv_zoom, Z_zoom, cmap='Blues', norm=normalize_is,zorder=0) # faster than pcolor, rasterized=True
heatmap_vir = ax_ins.pcolormesh(xv_zoom, yv_zoom, vir_zoom, cmap='Reds', norm=normalize_vir,zorder=0) # faster than pcolor
#heatmap = ax_ins.imshow(vir_dens2d, interpolation='none', rasterized=True, extent=[xmin_tot, xmax_tot, ymin_tot, ymax_tot]) # faster than pcolor

cbar = plt.colorbar(heatmap)
cbar.set_label(r'$N_h h(\mathbf{x})$', rotation=270, labelpad=+9)
#~ plt.colorbar(heatmap_vir)

#~ s_m = matplotlib.cm.ScalarMappable(cmap='RdBu', norm=normalize_upd)
#~ s_m = matplotlib.cm.ScalarMappable(cmap='RdBu', norm=matplotlib.colors.SymLogNorm(linthresh=1., linscale=1.,  vmin=0, vmax=min_maxcol))
#~ s_m.set_array([])
#~ plt.colorbar(s_m,  ticks=[10.,updmax/2.])

#    ax_ins.plot(IS_x, IS_y, linestyle='-', color='r', label='average IS')
ax_ins.set_xlabel('phenotypic trait 1')
ax_ins.set_ylabel('phenotypic trait 2')
#~ ax_ins.set_xlim(xmin_zoom, xmax_zoom)
#~ ax_ins.set_ylim(ymin_zoom, ymax_zoom)

# shitty hack
ax_ins.set_xlim(1500, xmax_zoom + 10)
ax_ins.set_ylim(ymin_zoom, -1880)


frame1 = plt.gca()

frame1.axes.get_xaxis().set_ticks([])
frame1.axes.get_yaxis().set_ticks([])

ylims_diff = ax_ins.get_ylim()[1] - ax_ins.get_ylim()[0]
xlims_diff = ax_ins.get_xlim()[1] - ax_ins.get_xlim()[0]

print xlims_diff, ylims_diff


if ylims_diff > xlims_diff:
    
    diff_diff= ylims_diff-xlims_diff
    
    maxdiff=ylims_diff
    
    xmin=ax_ins.get_xlim()[0] - diff_diff/2.
    xmax=ax_ins.get_xlim()[1] + diff_diff/2.
    #xmin=ax_ins.get_xlim()[0]  - diff_diff
    #xmax=ax_ins.get_xlim()[1]
    
    ax_ins.set_xlim(xmin, xmax)
    

elif xlims_diff > ylims_diff:

    diff_diff= xlims_diff-ylims_diff
    
    ymin=ax_ins.get_ylim()[0] - diff_diff/2.
    ymax=ax_ins.get_ylim()[1] + diff_diff/2.
    #~ ymin=ax_ins.get_ylim()[0]
    #~ ymax=ax_ins.get_ylim()[1]  + diff_diff
    
    ax_ins.set_ylim(ymin, ymax)
    maxdiff=xlims_diff
    



#ax_ins.set_xlim(xmin, xmax)
#ax_ins.set_ylim(ymin, ymax)

print maxdiff, rec_width

if maxdiff>2*rec_width:
    scale_size=rec_width
    scale_lab=r'$r$'
elif maxdiff>rec_width:
    scale_size=0.5*rec_width
    scale_lab=r'$0.5 r$'
elif maxdiff>0.2*rec_width:
    scale_size=0.1*rec_width
    scale_lab=r'$0.1 r$'
elif maxdiff>0.1*rec_width:
    scale_size=0.05*rec_width
    scale_lab=r'$0.05 r$'

print scale_size, scale_lab
        
scalebar = AnchoredSizeBar(ax_ins.transData,
               scale_size, scale_lab, 2, 
               pad=0.1,
               color='black',
               frameon=False,
               size_vertical=0.,
               label_top=True)# ,fontproperties=fontprops
               

ax_ins.add_artist(scalebar)


if maxdiff>2*vtau:
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
        
scalebar = AnchoredSizeBar(ax_ins.transData,
               scale_size, scale_lab, 3, 
               pad=0.1,
               color='black',
               frameon=False,
               size_vertical=0.,
               label_top=True)# ,fontproperties=fontprops
               

ax_ins.add_artist(scalebar)


ax_ins.xaxis.labelpad = axis_labelpad
ax_ins.yaxis.labelpad = axis_labelpad
ax_ins.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
  loc='upper right', bbox_to_anchor=(2.5, 1.2))
mpsetup.despine(ax_ins) 




#~ file_in_frames_npz_compr = "{inp}/antigenic_space_time_49999.dat".format(inp=dir_io)

xs = []
ys = []
fs = []
nums = []
times = []

avg_xs = []
avg_ys = []

#centroids_dict={}


for file_in_frames_npz_compr in glob.glob("{inp}/antigenic_space_time_*.dat".format(inp=dir_io)):
    filename, _ = os.path.splitext(file_in_frames_npz_compr)
    time=os.path.basename(filename).split('_')[-1]
    time=int(time)
    data_space = np.loadtxt(file_in_frames_npz_compr)
    #     outstream<<"# 1 x"<<setw(30)<<" 2 y "<<setw(30)<<" 3 number IS"<<setw(30)<<"4 number viruses"<<setw(30)<<"5 viral fitness" <<endl;
    #~ print data_space.ndim
            
    if data_space.ndim >1 and data_space[data_space[:,3]>0,:].ndim >1:

        data_space = data_space[data_space[:,3]>0,:]
        
        #print data_space
        #~ print data_space.shape
        
        if data_space.ndim == 1:     
            vir_x   = data_space[0] # 
            vir_y   = data_space[1] # 
            num_vir = data_space[3] # 
        else:    
            #dim_reshape=min(13, data_space.shape[1])
            #data_space=data_space[:,:dim_reshape]
            vir_x = data_space[:,0] # 
            vir_y = data_space[:,1] # 
            num_vir = data_space[:,3]

        fitn = data_space[:,4] # 

        #vir_y=vir_y[vir_x<31]
        #vir_x=vir_x[vir_x<31]
        
        #vir_x =vir_x[num_vir>0]
        #vir_y =vir_y[num_vir>0]	
        xs.extend([np.amin(vir_x),np.amax(vir_x)])
        ys.extend([np.amin(vir_y),np.amax(vir_y)])
        fs.extend([np.amin(fitn),np.amax(fitn)])
        nums.extend([np.amax(num_vir)])
        
        times.append(time)
        
        vir_x_marginal= np.unique(vir_x)
        num_vir_x_marginal= np.asarray([ np.sum(num_vir[x==vir_x])  for x in vir_x_marginal])
        
        avg_x=np.average(vir_x_marginal, weights=num_vir_x_marginal) # lazy way to compute distribution average
        avg_xs.append(avg_x)
        
        
        vir_y_marginal= np.unique(vir_y)
        num_vir_y_marginal= np.asarray([ np.sum(num_vir[y==vir_y])  for y in vir_y_marginal])
        
        avg_y=np.average(vir_y_marginal, weights=num_vir_y_marginal) # lazy way to compute distribution average
        avg_ys.append(avg_y)
        
    
#times.sort()

times, avg_xs, avg_ys= zip(*sorted(zip(times, avg_xs, avg_ys)))    
times =np.asarray(times)	
avg_xs=np.asarray(avg_xs)	
avg_ys=np.asarray(avg_ys)	

print times
print times.shape


fmin=min(fs)
fmax=max(fs)

N_lev=20
fit_cont_lev= np.linspace(fmin, fmax, N_lev)[1:-1]

ind0 = np.searchsorted(fit_cont_lev, 0., side='left')

if fit_cont_lev[ind0] !=0:
    fit_cont_lev=np.insert(fit_cont_lev, ind0, 0.)

num_max=max(nums)


#~ normalize = matplotlib.colors.LogNorm(vmin=0.1, vmax=num_max)
    
for file_in_frames_npz_compr in glob.glob("{inp}/antigenic_space_time_*.dat".format(inp=dir_io)):
    filename, _ = os.path.splitext(file_in_frames_npz_compr)
    time=os.path.basename(filename).split('_')[-1]
    time=int(time)
    #print(frame_number)
    
    #print(file_in_viruses)
    
    data_space = np.loadtxt(file_in_frames_npz_compr)
    #     outstream<<"# 1 x"<<setw(30)<<" 2 y "<<setw(30)<<" 3 number IS"<<setw(30)<<"4 number viruses"<<setw(30)<<"5 viral fitness" <<endl;
   
    #~ print data_space.ndim
   
    if data_space.ndim >1 and data_space[data_space[:,3]>0,:].ndim >1:
       
        data_space = data_space[data_space[:,3]>0,:]
        
        
        #print data_space
        #~ print data_space.shape
        
        if data_space.ndim == 1:     
            vir_x = data_space[0] # 
            vir_y = data_space[1] # 
            num_vir = data_space[3] # 
        else:    
            #dim_reshape=min(13, data_space.shape[1])
            #data_space=data_space[:,:dim_reshape]
            vir_x = data_space[:,0] # 
            vir_y = data_space[:,1] # 
            num_vir = data_space[:,3]
        
        fitn = data_space[:,4] # 

        #vir_x   =vir_x[num_vir>0]
        #vir_y   =vir_y[num_vir>0]	
        #num_vir =num_vir[num_vir>0]	
        
        #~ print num_vir.shape
        
            
        xmin_zoom=np.amin(vir_x)
        xmax_zoom=np.amax(vir_x)
        ymin_zoom=np.amin(vir_y)
        ymax_zoom=np.amax(vir_y)
        
        if xmax_zoom - xmin_zoom > ymax_zoom - ymin_zoom:
             ymax_zoom += ((xmax_zoom - xmin_zoom) - (ymax_zoom - ymin_zoom))/2.
             ymin_zoom -= ((xmax_zoom - xmin_zoom) - (ymax_zoom - ymin_zoom))/2.
        else:
             xmax_zoom += ((ymax_zoom - ymin_zoom) - (xmax_zoom - xmin_zoom) )/2.
             xmin_zoom -= ((ymax_zoom - ymin_zoom) - (xmax_zoom - xmin_zoom) )/2.
        
    
        
        nx_zoom          = int((xmax_zoom - xmin_zoom)/(latt_sp)) +1
        ny_zoom          = int((ymax_zoom - ymin_zoom)/(latt_sp)) +1  
        x_zoom           = np.linspace(xmin_zoom, xmax_zoom, nx_zoom)
        y_zoom           = np.linspace(ymin_zoom, ymax_zoom, ny_zoom)
        xv_zoom, yv_zoom = np.meshgrid(x_zoom, y_zoom)    
            

        Z_zoom = np.zeros(xv_zoom.shape, dtype=float) 
        fitn_zoom = np.zeros(xv_zoom.shape, dtype=float) - 1000 #+ 0.1
        #fitn_zoom = ma.masked_all(xv_zoom.shape, dtype=float)  #+ 0.1 - 1000
        
        #~ print Z_zoom.shape
        
        
        for i_x, x_f in enumerate(vir_x):
            
            y_f   = vir_y[i_x]
            num_f = num_vir[i_x]
            fit_f = fitn[i_x]
            
            
            idx_x_zoom = int((x_f - xmin_zoom)/latt_sp)
            idx_y_zoom = int((y_f - ymin_zoom)/latt_sp)
            
            #print x_f, xmax , xmin, idx_x
            #print y_f, ymax, ymin , idx_y
            #print nx, ny

            
            
            Z_zoom[idx_y_zoom, idx_x_zoom] += num_f
            if fit_f!=0: #temporary workaraund
                if fitn_zoom[idx_y_zoom, idx_x_zoom] == -1000:
                    fitn_zoom[idx_y_zoom, idx_x_zoom] = 0.
                    
                fitn_zoom[idx_y_zoom, idx_x_zoom] += fit_f
                
        fitn_zoom=np.ma.masked_where(fitn_zoom <= -1000, fitn_zoom)


        idx_time= np.nonzero(times == time)[0]
        
        #~ print idx_time, time
        
        if idx_time==0:
            centr_son_x=avg_xs[idx_time+1]
            centr_son_y=avg_ys[idx_time+1]

            centr_prec_x=avg_xs[idx_time]
            centr_prec_y=avg_ys[idx_time]
        
        else:
            centr_son_x=avg_xs[idx_time]
            centr_son_y=avg_ys[idx_time]

            centr_prec_x=avg_xs[idx_time-1]
            centr_prec_y=avg_ys[idx_time-1]
            
            
        centr_son_x =np.squeeze(centr_son_x )
        centr_son_y =np.squeeze(centr_son_y )    
        centr_prec_x=np.squeeze(centr_prec_x)
        centr_prec_y=np.squeeze(centr_prec_y)  
        
        
          
        #~ print centr_son_x, centr_son_y
        #~ print centr_prec_x, centr_prec_y

        clust_dir=np.array([centr_son_x - centr_prec_x , centr_son_y - centr_prec_y])/np.sqrt((centr_prec_x - centr_son_x)**2 + (centr_prec_y - centr_son_y)**2)

        
        vars_dist = clusters_dist_vars_all(data_space, clust_dir)
        
        var_parall= vars_dist[0]
        var_perp  = vars_dist[1]
        var_tot   = vars_dist[2]
        
        #var_parall_cloud.append(var_parall)
        #var_perp_cloud.append(var_perp    )
        #var_tot_cloud.append(var_tot      )
                    
        
        #~ var_parall_cloud[idx_time] = var_parall
        #~ var_perp_cloud  [idx_time] = var_perp  
        #~ var_tot_cloud   [idx_time] = var_tot   
                 
                    
                    
        # put top 5 % coords in black
          
        
        
        k=int(fitn.size*2./100)
        if k==0:
            k=1
        ind = np.argpartition(fitn, -k)[-k:] # average the 6 biggest distances to baricenter of core samples
        top_vir_x=vir_x[ind]
        top_vir_y=vir_y[ind]

                
                    
                    


# Create an inset outside the axes
#~ axins = inset_axes(ax_ins, width="100%", height="100%",
                   #~ bbox_to_anchor=(1.2, 0., 1., 1.),
                   #~ bbox_transform=ax_ins.transAxes, loc=2, borderpad=0)
#~ axins.tick_params(left=False, right=False, labelleft=False, labelright=False)
 
axins = plt.Subplot(fig_ins, grid_ins[1, 1])#, rasterized=True
fig_ins.add_subplot(axins)
labeled_axes.append(axins)
  



        
numvirmax=np.amax(num_vir)
            
            
normalize_vir = matplotlib.colors.LogNorm(vmin=1, vmax=numvirmax)
    
heatmap2 = axins.pcolormesh(xv_zoom, yv_zoom, Z_zoom, cmap='Reds', norm=normalize_vir) # faster than pcolor
#heatmap = axins.imshow(vir_dens2d, interpolation='none', rasterized=True, extent=[xmin_tot, xmax_tot, ymin_tot, ymax_tot]) # faster than pcolor

        
cbar = plt.colorbar(heatmap2, ax = axins)
cbar.set_label(r'$n(\mathbf{x})$', rotation=270, labelpad=+9)

fit_cont_lev= np.linspace(fmin, fmax, N_lev)[1:-1]

indfmin = np.searchsorted(fit_cont_lev, np.amin(fitn), side='right')

indfmax = np.searchsorted(fit_cont_lev, np.amax(fitn), side='left')

fit_cont_lev_frame=fit_cont_lev[indfmin:indfmax]


CS = plt.contour(xv_zoom, yv_zoom, fitn_zoom, fit_cont_lev_frame, linewidths=.5, colors='k' ) # negative contours will be dashed by default
            
axins.plot(top_vir_x, top_vir_y + latt_sp*0.5, linestyle='', marker='.', color='k', markersize=2.7)#, label='top viruses',rasterized=True

#    axins.plot(IS_x, IS_y, linestyle='-', color='r', label='average IS')
axins.set_xlabel('phenotypic trait 1')
axins.set_ylabel('phenotypic trait 2')
axins.set_xlim(xmin_zoom -1, xmax_zoom +2)
axins.set_ylim(ymin_zoom, ymax_zoom)

axins.arrow(centr_son_x, centr_son_y, 0.2*(xmax_zoom - xmin_zoom)*clust_dir[0], 0.2*(xmax_zoom - xmin_zoom)*clust_dir[1], head_width=0.02*(xmax_zoom - xmin_zoom), head_length=0.04*(xmax_zoom - xmin_zoom), length_includes_head=True, fc='k', ec='k')
#axins.annotate("", xy=(0.5, 0.5), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))

#axins.axis('equal') 





# shitty hack
#~ axins.set_xlim(1500, xmax_zoom + 10)
#~ axins.set_ylim(ymin_zoom, -1880)


ylims_diff = axins.get_ylim()[1] - axins.get_ylim()[0]
xlims_diff = axins.get_xlim()[1] - axins.get_xlim()[0]

print xlims_diff, ylims_diff


if ylims_diff > xlims_diff:
    
    diff_diff= ylims_diff-xlims_diff
    
    maxdiff=ylims_diff
    
    #~ xmin=axins.get_xlim()[0] - diff_diff/2.
    #~ xmax=axins.get_xlim()[1] + diff_diff/2.
    xmin=axins.get_xlim()[0]  - diff_diff
    xmax=axins.get_xlim()[1]
    
    axins.set_xlim(xmin, xmax)
    

elif xlims_diff > ylims_diff:

    diff_diff= xlims_diff-ylims_diff
    
    #~ ymin=axins.get_ylim()[0] - diff_diff/2.
    #~ ymax=axins.get_ylim()[1] + diff_diff/2.
    ymin=axins.get_ylim()[0] - diff_diff*0.6
    ymax=axins.get_ylim()[1] + diff_diff*0.4
    #~ ymin=axins.get_ylim()[0]
    #~ ymax=axins.get_ylim()[1]  + diff_diff
    
    axins.set_ylim(ymin, ymax)
    maxdiff=xlims_diff
    



#axins.set_xlim(xmin, xmax)
#axins.set_ylim(ymin, ymax)

print maxdiff, rec_width

if maxdiff>2*rec_width:
    scale_size=rec_width
    scale_lab=r'$r$'
elif maxdiff>rec_width:
    scale_size=0.5*rec_width
    scale_lab=r'$0.5 r$'
elif maxdiff>0.2*rec_width:
    scale_size=0.1*rec_width
    scale_lab=r'$0.1 r$'
elif maxdiff>0.1*rec_width:
    scale_size=0.05*rec_width
    scale_lab=r'$0.05 r$'
elif maxdiff>0.02*rec_width:
    scale_size=0.01*rec_width
    scale_lab=r'$10^{-2} r$'
elif maxdiff>0.01*rec_width:
    scale_size=0.005*rec_width
    #~ scale_lab=r'$5*10^{-3} r$'
    scale_lab=r'$0.005 r$'
elif maxdiff>0.002*rec_width:
    scale_size=0.001*rec_width
    scale_lab=r'$10^{-3} r$'
elif maxdiff>0.001*rec_width:
    scale_size=0.0005*rec_width
    scale_lab=r'$5 10^{-4} r$'

elif maxdiff>0.0002*rec_width:
    scale_size=0.0001*rec_width
    scale_lab=r'$10^{-4} r$'

print scale_size, scale_lab
        
scalebar = AnchoredSizeBar(axins.transData,
               scale_size, scale_lab, 1, 
               pad=0.1,
               color='black',
               frameon=False,
               size_vertical=0.,
               label_top=True)# ,fontproperties=fontprops
               

axins.add_artist(scalebar)


if maxdiff>2*vtau:
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
        
scalebar = AnchoredSizeBar(axins.transData,
               scale_size, scale_lab, 4, 
               pad=0.1,
               color='black',
               frameon=False,
               size_vertical=0.,
               label_top=True)# ,fontproperties=fontprops
               

axins.add_artist(scalebar)


#~ axins.xaxis.labelpad = axis_labelpad
#~ axins.yaxis.labelpad = axis_labelpad
#~ axins.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
  #~ loc='upper right', bbox_to_anchor=(2.5, 1.2))
#~ mpsetup.despine(axins) 

 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax_ins, axins, loc1=2, loc2=3, fc="none", ec="0.2")#
 

frame1 = fig_ins.gca()

#~ frame1.axes.get_xaxis().set_ticks([])
#~ frame1.axes.get_yaxis().set_ticks([])

frame1.axes.get_xaxis().set_visible(False)
frame1.axes.get_yaxis().set_visible(False)

#### panel D ####


ax = plt.Subplot(fig_ins, grid_ins[0, 1])#, rasterized=True
fig_ins.add_subplot(ax)
labeled_axes.append(ax)
  
  

file_in_frames_npz_compr = "{inp}/antigenic_space_time_49999.dat".format(inp=dir_io)

filename, _ = os.path.splitext(file_in_frames_npz_compr)
time=os.path.basename(filename).split('_')[-1]
time=int(time)
#print(frame_number)

#print(file_in_viruses)

data_space = np.loadtxt(file_in_frames_npz_compr)
#     outstream<<"# 1 x"<<setw(30)<<" 2 y "<<setw(30)<<" 3 number IS"<<setw(30)<<"4 number viruses"<<setw(30)<<"5 viral fitness" <<endl;

#~ print data_space.ndim

if data_space.ndim >1 and data_space[data_space[:,3]>0,:].ndim >1:
   
    data_space = data_space[data_space[:,3]>0,:]
    
    
    
    fitn = data_space[:,4] # 


# Create an inset outside the axes
#~ axins = inset_axes(ax_ins, width="100%", height="100%",
                   #~ bbox_to_anchor=(1.2, 0., 1., 1.),
                   #~ bbox_transform=ax_ins.transAxes, loc=2, borderpad=0)
#~ axins.tick_params(left=False, right=False, labelleft=False, labelright=False)
 
 
 
 
#    axins.plot(IS_x, IS_y, linestyle='-', color='r', label='average IS')
ax.set_xlabel('fitness')
ax.set_ylabel('probability density')

window=1

avg=np.mean(fitn)#/33.

stddev=np.std(fitn)

n1, bins1 = np.histogram(fitn, 50, density=False) # , normed=True

n_smooth1=np.convolve(n1, np.ones((window,))/window, mode='same')

center1 = (bins1[:-1] + bins1[1:]) / 2

indmax=np.argmax(n_smooth1)

ax.plot(center1, n_smooth1, c='r',  linestyle='-', label='raw data')

ax.arrow(center1[indmax], n_smooth1[indmax], 1.*stddev, 0., head_width=4*stddev, head_length=0.2*stddev, length_includes_head=True, fc='k', ec='k')
#axins.annotate("", xy=(0.5, 0.5), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))


mean,std=norm.fit(fitn)
y = norm.pdf(center1, mean, std)
ax.plot(center1, y, c='g',  linestyle='-', label='gaussian fit')
ax.set_xlim(-0.11, 0.11)
ax.set_ylim(ymin=0, ymax=17)





ax = plt.Subplot(fig_ins, grid_ins[0, 0])#, rasterized=True
fig_ins.add_subplot(ax)
labeled_axes.append(ax)
  
  
frame1 = fig_ins.gca()

frame1.axes.get_xaxis().set_visible(False)
frame1.axes.get_yaxis().set_visible(False)

# ~ ax.set_visible(False)

  
#### finish figure ####




labeldict = dict(labelstyle=r'%s', fontsize='large',
         xycoords=('axes fraction'), fontweight = 'bold')
mpsetup.label_axes([labeled_axes[0]], labels='B', xy=(-0.15,  0.97), **labeldict)
mpsetup.label_axes([labeled_axes[1]], labels='C', xy=(-0.15,  0.97), **labeldict)
mpsetup.label_axes([labeled_axes[2]], labels='D', xy=(-0.12,  0.97), **labeldict)
mpsetup.label_axes([labeled_axes[3]], labels='A', xy=(-0.12,  0.97), **labeldict)

out_file='{out}/fig1all.pdf'.format(out=dir_io)
#    print out_file
#~ fig_ins.savefig(out_file, dpi=3000)
fig_ins.savefig(out_file, dpi=1000)
fig_ins.clf()
plt.close('all')
   
