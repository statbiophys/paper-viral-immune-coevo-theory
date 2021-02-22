
"""
CLUSTERS AND TRACKS VIRAL LINEAGES FROM TIME SNAPSHOTS. SAVE A FILE WITH LINEAGES TRACKING IN OUTPUT, TO BE USED BY DOWNSTREAM SCRIPTS

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
from scipy.interpolate import splrep, splev, griddata, SmoothBivariateSpline 
from scipy import stats
from sklearn.neighbors import KDTree
from scipy.optimize import curve_fit
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


def get_cmap(N):
    ''' Returns a function that maps each index in 0, 1, ...
        N-1 to a distinct RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    #scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap='jet')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color




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

def rw_msd_fct_time(x, diff):
    return diff*x




def mean_dispersion(X, labels): # returns a float: the mean dispersion averaged over all the clusters
    labels_unique=set(labels)
    labels_unique.discard(-1)
    clust_dispersions=[]
    for clust in labels_unique:
        points_in_clust=X[labels==clust,:]
        squared_dists=np.linalg.norm(points_in_clust - np.mean(points_in_clust, axis=0)[None,:], axis=1)**2
        #print clust, points_in_clust.shape, squared_dists.shape
        mean_squared_dist=np.mean(squared_dists)
        #print mean_squared_dist
        clust_dispersions.append(mean_squared_dist)
        
    if len(labels_unique)>0:
        return [sum(clust_dispersions)/float(len(labels_unique)), sum(clust_dispersions)]
    else:
        return [-1., sum(clust_dispersions)]




def clusters_size(X, labels): # returns a float: the mean dispersion averaged over all the clusters
    labels_unique=set(labels)
    labels_unique.discard(-1)
    clust_dispersions=[]
    densities=[]
    for clust in labels_unique:
        points_in_clust=X[labels==clust,:]
        num=points_in_clust.shape[0]
        dists=np.linalg.norm(points_in_clust - np.mean(points_in_clust, axis=0)[None,:], axis=1)
        #print clust, points_in_clust.shape, squared_dists.shape
        mean_dist=np.mean(dists)
        #print mean_squared_dist
        clust_dispersions.append(np.pi*(mean_dist**2))
	if mean_dist>0.1 and num>=2:
	    densities.append(num/(np.pi*(mean_dist**2)))
	else:
	    densities.append(0.)
	
    if len(labels_unique)>0:
        return [sum(clust_dispersions)/float(len(labels_unique)), sum(densities)/float(len(labels_unique))]
    else:
        return [-1., -1.]





def clusters_dist_vars_all(X, clust_dir): # returns three floats: the variance of distances from the centroid , parrallel to the direction of motion, perpendicular, and total
    
	#~ print "cluster rotation"


	points_in_clust=X[:,:2]
    
    
	#~ print np.mean(points_in_clust, axis=0)
	#~ print np.average(points_in_clust, weights=X[:,3], axis=0)
	
	#points_in_clust=points_in_clust - np.mean(points_in_clust, axis=0)[None,:] # regularize
	points_in_clust=points_in_clust - np.average(points_in_clust, weights=X[:,3], axis=0)[None,:] # regularize

	#~ print points_in_clust.shape
	#~ print np.average(points_in_clust, weights=X[:,3], axis=0)


	clust_perpdir=np.array([clust_dir[1], -clust_dir[0]])
	
	Q = np.squeeze(np.array([clust_dir, clust_perpdir])).T
	#~ print 'Transformed to original system with \n Q={}'.format(Q)
	#~ print 'Orthogonality check \n {}'.format(Q.dot(Q.T))
	#~ print 'Orthogonality check \n {}'.format(np.dot(Q.T, X[0,:2]) == np.linalg.solve(Q, X[0,:2]))
	
	#~ print np.dot(Q.T, points_in_clust[0,:])
	#~ print np.linalg.solve(Q, points_in_clust[0,:])
	
	points_in_clust_newbasis=np.dot(Q.T, points_in_clust.T).T
	
	#~ print "transformed direction of motion ", np.dot(Q.T, clust_dir)
	#~ print "transformed direction perpendicular to motion ", np.dot(Q.T, clust_perpdir)
	
	#~ print points_in_clust_newbasis.shape, np.mean(points_in_clust_newbasis, axis=0)
	
	# v_ = np.linalg.solve(Q, v)
	# print('The vector in the new coordinates \n v_={}'.format(v_))
	

	#dists_sq=np.linalg.norm(points_in_clust_newbasis - np.mean(points_in_clust_newbasis, axis=0)[None,:], axis=1)**2
	dists_sq=np.linalg.norm(points_in_clust_newbasis , axis=1)**2  
	#~ print points_in_clust.shape, dists_sq.shape
	#var_dist=np.mean(dists_sq)
	var_dist=np.average(dists_sq, weights=X[:,3])
	#~ print var_dist
	
	dists_sq_x=points_in_clust_newbasis[:,0]**2  
	#~ print dists_sq_x.shape
	var_dist_x=np.average(dists_sq_x, weights=X[:,3])
	#~ print var_dist_x
	
	dists_sq_y=points_in_clust_newbasis[:,1]**2  
	#~ print dists_sq_y.shape
	var_dist_y=np.average(dists_sq_y, weights=X[:,3]) # direction of motion
	#~ print var_dist_y
	
	
	return [var_dist_x, var_dist_y, var_dist]



def clusters_dist_vars(X, labels, clust, clust_dir, num_vir): # returns three floats: the variance of distances from the centroid , parrallel to the direction of motion, perpendicular, and total
    
    #~ print "cluster rotation"
    
    
    points_in_clust=X[labels==clust,:]
    num_in_clust=num_vir[labels==clust]
    
    #~ print np.average(points_in_clust, weights=num_in_clust, axis=0)
    
    points_in_clust=points_in_clust - np.average(points_in_clust, weights=num_in_clust, axis=0)[None,:] # regularize
    
    #~ print points_in_clust.shape
    #~ print np.mean(points_in_clust, axis=0)
    #~ print np.average(points_in_clust, weights=num_in_clust, axis=0)
    
    
    
    clust_perpdir=np.array([clust_dir[1], -clust_dir[0]])
    
    Q = np.array([clust_dir, clust_perpdir]).T
    #~ print 'Transformed to original system with \n Q={}'.format(Q)
    #~ print 'Orthogonality check \n {}'.format(Q.dot(Q.T))
    #~ print 'Orthogonality check \n {}'.format(np.dot(Q.T, X[0,:]) == np.linalg.solve(Q, X[0,:]))
    
    #~ print np.dot(Q.T, points_in_clust[0,:])
    #~ print np.linalg.solve(Q, points_in_clust[0,:])
    
    points_in_clust_newbasis=np.dot(Q.T, points_in_clust.T).T
    
    #~ print "transformed direction of motion ", np.dot(Q.T, clust_dir)
    #~ print "transformed direction perpendicular to motion ", np.dot(Q.T, clust_perpdir)
    
    #~ print points_in_clust_newbasis.shape, np.mean(points_in_clust_newbasis, axis=0)
    
    # v_ = np.linalg.solve(Q, v)
    # print('The vector in the new coordinates \n v_={}'.format(v_))
    
    
    #dists_sq=np.linalg.norm(points_in_clust_newbasis - np.mean(points_in_clust_newbasis, axis=0)[None,:], axis=1)**2
    dists_sq=np.linalg.norm(points_in_clust_newbasis , axis=1)**2  
    #~ print clust, points_in_clust.shape, dists_sq.shape
    #var_dist=np.mean(dists_sq)
    var_dist=np.average(dists_sq, weights=num_in_clust)
    #~ print var_dist
    
    dists_sq_x=points_in_clust_newbasis[:,0]**2  
    #~ print clust, dists_sq_x.shape
    var_dist_x=np.average(dists_sq_x, weights=num_in_clust)
    #~ print var_dist_x
    
    dists_sq_y=points_in_clust_newbasis[:,1]**2  
    #~ print clust, dists_sq_y.shape
    var_dist_y=np.average(dists_sq_y, weights=num_in_clust) # direction of motion
    #~ print var_dist_y
    
    
    return [var_dist_x, var_dist_y, var_dist]
    
    
    
    
    
def inter_clusters_distance(X, labels): # returns a float: the mean dispersion averaged over all the clusters
    labels_unique=set(labels)
    labels_unique.discard(-1)
    centroids=[]
    for clust in labels_unique:
        points_in_clust=X[labels==clust,:]
        num=points_in_clust.shape[0]
        centroid=np.mean(points_in_clust, axis=0)
    
        #print mean_squared_dist
        centroids.append(centroid)
            
        
    #print len(centroids), len(labels_unique)
    
    if len(centroids)<2:
        inter_cl=0.
        inter_cl_max=0.
    
    else:
        centroids=np.asarray(centroids)
        #print centroids.shape
        dists=pdist(centroids)
        #print dists.shape 
        
        inter_cl=dists.mean()
        inter_cl_max=np.amax(dists)
        
    return [inter_cl, inter_cl_max] 
    




def recluster_split(X, labels, core_samples_mask, num_vir): # resets the clustering labels so that splits are identified only when I see clusters splitted in the frame. 
    #Returns new labels and the corresponding centroids, and all basic cluster characteristics: number ppl inside, size
    labels_unique=set(labels)
    labels_unique.discard(-1)
    labels_clust=[]
    centroids=[]
    max_centroid_dists=[]
    max_centroid_dists=[]
    nums =[]
    sizes=[]
    clusts=[]
    coords=[]
    
    #print labels[core_samples_mask]
    #~ print "reclustering"
    
    #print labels.shape
    #print core_samples_mask.shape
    #~ print len(labels_unique)
    
    for clust in labels_unique:
        #print clust
        #print core_samples_mask[labels==clust].sum()
        #print core_samples_mask[:20]
        #print (labels==clust)[:20]
        #print (core_samples_mask & (labels==clust)).sum()
        #print (core_samples_mask[:20] & (labels==clust)[:20])
        
        points_in_clust=X[labels==clust,:]
        coord=points_in_clust.shape[0]
        num=num_vir[labels==clust]
        
        #centroid=np.mean(points_in_clust, axis=0)
        centroid=np.average(points_in_clust, weights=num, axis=0)
        
        core_points_in_clust=X[(labels==clust) & core_samples_mask,:]
    
        core_dists=core_points_in_clust - np.average(points_in_clust, weights=num, axis=0)[None,:]
        #
        #k=num
        #if num>10:
        #    k=10
        
        dists=np.linalg.norm(points_in_clust - np.average(points_in_clust, weights=num, axis=0)[None,:], axis=1)
        #print clust, points_in_clust.shape, squared_dists.shape
        #mean_dist=np.mean(dists)
        mean_dist=np.average(dists, weights=num)
        #print mean_squared_dist
        size=np.pi*(mean_dist**2)
        
    
        modulo_dists=np.linalg.norm(core_dists, axis=1)
        #ind = np.argpartition(modulo_dists, -k)[-k:] # average the 6 biggest distances to baricenter of core samples
        if modulo_dists.size<2:
            modulo_dists=dists
        
        maxdist=np.amax(modulo_dists)
        #maxdist=np.mean(modulo_dists[ind])
        
        #~ print maxdist
        #~ print centroid
        #print core_points_in_clust[np.argmax(modulo_dists),:]
        #print num, np.argmax(modulo_dists)
    
    
        #print mean_squared_dist
        centroids.append(centroid)
        max_centroid_dists.append(maxdist)
        labels_clust.append(clust)
        nums.append(np.sum(num) )
        sizes.append(size)
        clusts.append(clust)
        coords.append(coord)
            
            
    #~ print len(centroids), len(labels_unique)
    
    if len(centroids)<2:
        labels_split=labels
        centroids_split=centroids
        nums_split=nums 
        sizes_split=sizes
        clust_labels_split=clusts
        nums_coords_split=coords
        return [labels_split, centroids_split, nums_split, sizes_split, clust_labels_split, nums_coords_split]
	
    else:
        
        #merge=np.zeros((len(labels_clust),len(labels_clust)), dtype=bool) # all false by default, do not merge unless...
        tuples_to_merge=[] # list of tuple of indexes to merge.
        cm=0
        for (lab1_idx, lab2_idx) in itertools.combinations(range(len(labels_clust)), r=2): # check which original clusters are to merge
            
            lab1     =labels_clust[lab1_idx]
            centroid1=centroids[lab1_idx]
            max_centroid_dist1=max_centroid_dists[lab1_idx]
            
            lab2     =labels_clust[lab2_idx]
            centroid2=centroids[lab2_idx]
            max_centroid_dist2=max_centroid_dists[lab2_idx]
            
            #print centroid1
            #print centroid2
            #
            #
            #print max_centroid_dist1, max_centroid_dist2
            
            #max_centroid_dist_both=2*max([max_centroid_dist1, max_centroid_dist2])
            max_centroid_dist_both=max_centroid_dist1 + max_centroid_dist2
            
            dist_centroids=np.linalg.norm(centroid1 - centroid2)
            
            #merge[lab1_idx, lab2_idx]=  merge[lab2_idx, lab1_idx]= (dist_centroids<2*max_centroid_dist_both) # merge condition as square matrix
            
            #print dist_centroids, max_centroid_dist_both
            
            if (dist_centroids<max_centroid_dist_both + max([max_centroid_dist1, max_centroid_dist2])/2.):
                tuples_to_merge.append((lab1_idx, lab2_idx))
            
            
        #~ print "tuples to merge"
        #~ print len(tuples_to_merge)
        #~ print tuples_to_merge
        
        
        list_merges_tot=[]
            
        if len(tuples_to_merge)==0:
            labels_split=labels
            centroids_split=centroids
            nums_split=nums 
            sizes_split=sizes
            clust_labels_split=clusts
            nums_coords_split=coords    
            return [labels_split, centroids_split, nums_split, sizes_split, clust_labels_split, nums_coords_split]
    
        
        elif len(tuples_to_merge)==1:
            (lab1_idx, lab2_idx) = tuples_to_merge[0]
            lab1     =labels_clust[lab1_idx]
            lab2     =labels_clust[lab2_idx]
    
            list_merges_tot.append(set([lab1, lab2]))
            
        else:
            
            #for lab_idx, lab in enumerate(labels_clust): # build a list of clusters to merge all together
            
            for tup_idx1 in range(len(tuples_to_merge)):
                (lab1_idx, lab2_idx) = tuples_to_merge[tup_idx1]
        
                lab1     =labels_clust[lab1_idx]
                lab2     =labels_clust[lab2_idx]
                
                
                set_loop=set([lab1, lab2])
                
                if tup_idx1 < len(tuples_to_merge) -1:
                
                    for tup_idx2 in range(tup_idx1 +1,len(tuples_to_merge)):
                        
                        (lab1_idx2, lab2_idx2) = tuples_to_merge[tup_idx2]
                
                        lab1_sec     =labels_clust[lab1_idx2]
                        lab2_sec     =labels_clust[lab2_idx2]
                        
                        set_tmp=set([lab1_sec, lab2_sec])
                        
                        if len(set_loop.intersection(set_tmp))>0:
                            set_loop |= set_tmp
                    
                    
        
                if len(list_merges_tot)==0:
                    list_merges_tot.append(set_loop)
                else:
                    idxs_intersect=[i for i in range(len(list_merges_tot)) if set_loop.intersection(list_merges_tot[i])>0]
                    if len(idxs_intersect)==0:
                        list_merges_tot.append(set_loop)
                    elif len(idxs_intersect)==1:
                        list_merges_tot[idxs_intersect[0]] |= set_loop
                    else:
                        list_merges_tot[idxs_intersect[0]] |= set_loop
                        for idx_intersect in idxs_intersect[1:]:
                            list_merges_tot[idxs_intersect[0]] |= list_merges_tot[idx_intersect]
                            
                        list_merges_tot[:]=[list_merges_tot[i] for i in range(len(list_merges_tot)) if i not in idxs_intersect[1:]]
                
            
        labels_split=labels.copy()
        
        #~ print len(list_merges_tot)
        #~ print list_merges_tot
        
        for i_merge, set_merge in enumerate(list_merges_tot):
            
            if len(list_merges_tot)>1:
                list_merges_tot_others=list_merges_tot.copy()
                list_merges_tot_others.pop(i_merge)
                list_merges_tot_others=[list_merges_tot_others[0].update(merge) for merge in list_merges_tot_others]
                list_merges_tot_others=list_merges_tot_others[0]
                if set_merge.intersection(list_merges_tot_others)>0:
                    print "ERROR WITH SET WHEN MERGING CLUSTERS"
                    sys.exit()
            
            to_be_called=min(set_merge)
            idxs_change_label=[i for i in range(labels.size) if labels[i] in set_merge]
            
            labels_split[idxs_change_label] = to_be_called
            
        labels_unique=set(labels_split)
        labels_unique.discard(-1)
        centroids_split=[]
        nums_split=[]
        sizes_split=[]
        clust_labels_split=[]
        nums_coords_split=[]
        
        for clust in labels_unique:
            points_in_clust=X[labels_split==clust,:]
            num=num_vir[labels_split==clust]
            coord=points_in_clust.shape[0]

            
            #centroid=np.mean(points_in_clust, axis=0)
            centroid=np.average(points_in_clust, weights=num, axis=0)
    
        
            #dists=np.linalg.norm(points_in_clust - np.mean(points_in_clust, axis=0)[None,:], axis=1)
            dists=np.linalg.norm(points_in_clust - np.average(points_in_clust, weights=num, axis=0)[None,:], axis=1)
            #print clust, points_in_clust.shape, squared_dists.shape
            mean_dist=np.average(dists, weights=num)
            #print mean_squared_dist
            size=np.pi*(mean_dist**2)
    
    
            centroids_split.append(centroid)
            nums_split.append(np.sum(num))
            sizes_split.append(size)
            clust_labels_split.append(clust)
            nums_coords_split.append(coord)
                
    #~ print "reclustered"
		

  
    return [labels_split, centroids_split, nums_split, sizes_split, clust_labels_split, nums_coords_split]




def mean_dispersion_renorm(X, labels, core_samples_mask): # returns a float: the mean dispersion averaged over all the clusters
    labels_unique=set(labels)
    labels_unique.discard(-1)
    clust_dispersions=[]
    for clust in labels_unique:
        points_in_clust=X[labels==clust,:]
        core_points_in_clust=X[(labels==clust) & core_samples_mask,:]
        dists=points_in_clust - np.mean(points_in_clust, axis=0)[None,:]
        core_dists=core_points_in_clust - np.mean(points_in_clust, axis=0)[None,:]
        squared_dists=np.linalg.norm(dists, axis=1)**2
        
        
        num=core_points_in_clust.shape[0]
        
        k=num
        if num>10:
            k=10
        
        modulo_dists=np.linalg.norm(core_dists, axis=1)
        ind = np.argpartition(modulo_dists, -k)[-k:] # average the 6 biggest distances to baricenter of core samples
        #maxdist=np.amax(np.linalg.norm(core_dists, axis=1))
        maxdist=np.mean(modulo_dists[ind])
        
        #print clust, points_in_clust.shape, squared_dists.shape
        #renorm_squared_dist=np.sum(squared_dists)/(maxdist**2)
        renorm_squared_dist=np.mean(squared_dists)/(maxdist**2)
        #print mean_squared_dist
        clust_dispersions.append(renorm_squared_dist)
        
    #return sum(clust_dispersions)/float(labels[labels!=-1].size)
	
    if len(labels_unique)>0:
        return sum(clust_dispersions)/float(len(labels_unique))
    else:
        return -1.



def mean_dispersion_renorm_other(X, labels): # returns a float: the mean dispersion averaged over all the clusters
    labels_unique=set(labels)
    labels_unique.discard(-1)
    clust_dispersions=[]
    for clust in labels_unique:
        points_in_clust=X[labels==clust,:]
        dists=points_in_clust - np.mean(points_in_clust, axis=0)[None,:]
        squared_dists=np.linalg.norm(dists, axis=1)**2
        meandist=np.mean(np.linalg.norm(dists, axis=1))
        #print clust, points_in_clust.shape, squared_dists.shape
        #renorm_squared_dist=np.sum(squared_dists)/(meandist**2)
        renorm_squared_dist=np.mean(squared_dists)/(meandist**2) - 1.
        #print mean_squared_dist
        clust_dispersions.append(renorm_squared_dist)
        
    #return sum(clust_dispersions)/float(labels[labels!=-1].size)
	
    if len(labels_unique)>0:
        return sum(clust_dispersions)/float(len(labels_unique))
    else:
        return -1.




def kth_dist_distr(X, labels): # returns a float: the CV of pdf of the 10th nearest neighbor within cluster
    labels_unique=set(labels)
    labels_unique.discard(-1)
    clust_dispersions=[]
    variances=[]
    for clust in labels_unique:
        points_in_clust=X[labels==clust,:]
        num=points_in_clust.shape[0]
        
        k=num
        if num>10:
            k=10
        
        nbrs = NearestNeighbors(n_neighbors=k).fit(points_in_clust)
        distances, indices = nbrs.kneighbors(points_in_clust)
        #print distances.shape
        kth_dist=distances[:,-1]
        kth_dist_CV=np.std(kth_dist)/np.mean(kth_dist)
    
        clust_dispersions.append(kth_dist_CV)
        variances.append(np.std(kth_dist))
        
    #return sum(clust_dispersions)/float(labels[labels!=-1].size)
	
    if len(labels_unique)>0:
        return [sum(clust_dispersions)/float(len(labels_unique)), sum(variances)/float(len(labels_unique)), sum(variances)]
    else:
        return [-1., -1., sum(variances)]




def banfeld(X, labels): # returns a float: the mean dispersion averaged over all the clusters
    labels_unique=set(labels)
    labels_unique.discard(-1)
    clust_dispersions=[]
    for clust in labels_unique:
        points_in_clust=X[labels==clust,:]
        squared_dists=np.linalg.norm(points_in_clust - np.mean(points_in_clust, axis=0)[None,:], axis=1)**2
        #print clust, points_in_clust.shape, squared_dists.shape
        num=points_in_clust.shape[0]
        if num>0:
            banfeld=num*np.log(np.mean(squared_dists))
        else:
            banfeld=0
        #print mean_squared_dist
        clust_dispersions.append(banfeld)
	
    return sum(clust_dispersions)








def vtau_fromparams(recog_width, mem_points, F_0):
    return (recog_width/(np.exp(F_0/mem_points)-1.))        


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
mem_points=params['mem_points']
F0=params['F0']

latt_sp=1
vtau=vtau_fromparams(recog_width, mem_points, F0)


#if n_real >10:
#    n_real=10

print n_real   

thisfigsize = figsize
thisfigsize[1] *= 0.75

# get data

for real in np.arange(1,n_real+1):
    dir_in='{inp}/realization_{real}'.format(inp=dir_in_tot,real=real) # directory with input files
    dir_out_plots='{inp}/realization_{real}/clust'.format(inp=dir_out_plots_tot,real=real) # directory with output plots
    dir_out_frames_subs='{inp}/frames_subs'.format(inp=dir_out_plots) # directory with output plots
    dir_out_frames_zoom='{inp}/frames_zoom'.format(inp=dir_out_plots) # directory with output plots
    #dir_in_frames='{inp}/frames_npz_compr'.format(inp=dir_in,real=real) # directory with input files
    dir_in_frames='{inp}/frames/'.format(inp=dir_in) # directory with input files
    

    file_exploded='{inp}/expl_file.txt'.format(inp=dir_in)
    file_extinct='{inp}/extinct_file.txt'.format(inp=dir_in)
    
    
    
    
    file_in_avg_npz_compr='{inp}/evo_mean_stats_real_{real}.dat'.format(inp=dir_in, real=real)
    
    data=[]
    time_ss=[]
    if os.path.isfile(file_in_avg_npz_compr):
        
        
        
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

        
        survival=np.amax(time_ss)
        #time_ss=time[(time> t_init/365.) & (time <= survival)]
        #
        #print survival, t_init/365.
        print time_ss.size
        
        #if time_ss.size > 30000:
        #    
        #    data= data[:ind_times_notexpl,:]
        #    data= data[time_all>sec_thr,:]
            
                     
#    if os.path.exists(dir_in_frames):
    
    if os.path.exists(dir_in_frames) and data.ndim > 1 and time_ss.size > 3000:
    
        if not os.path.exists(dir_out_plots):
            os.makedirs(dir_out_plots)	
            
        if os.path.exists(dir_out_frames_subs):
            shutil.rmtree(dir_out_frames_subs, ignore_errors=False, onerror=None)
           
        if not os.path.exists(dir_out_frames_subs):
            os.makedirs(dir_out_frames_subs)	
            
           
       
       
       
       
            
        xs = []
        ys = []
        fs = []
        nums = []
        times = []
        
        avg_xs = []
        avg_ys = []
        
        #centroids_dict={}
    
        ######################################## LOAD DATA
        
        for file_in_frames_npz_compr in glob.glob("{inp}/antigenic_space_time_*.dat".format(inp=dir_in_frames, real=real)):
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
                
                if int(time)>sec_thr and int(time)<=survival and data_space.ndim > 1: 
            
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
        
        #~ print times
        print times.shape
        
        xmin=min(xs)
        xmax=max(xs)
        ymin=min(ys)
        ymax=max(ys)
        
        fmin=min(fs)
        fmax=max(fs)
        
        N_lev=20
        fit_cont_lev= np.linspace(fmin, fmax, N_lev)[1:-1]
        
        ind0 = np.searchsorted(fit_cont_lev, 0., side='left')
        
        if fit_cont_lev[ind0] !=0:
            fit_cont_lev=np.insert(fit_cont_lev, ind0, 0.)
        
        num_max=max(nums)
        
        
        normalize = matplotlib.colors.LogNorm(vmin=0.1, vmax=num_max)
            
        #if real<=10:
         
            
        var_parall_cloud = np.zeros((times.shape))
        var_perp_cloud   = np.zeros((times.shape))
        var_tot_cloud    = np.zeros((times.shape))
            
        xs = [] # to print in all times
        ys = []
        times_printed= []
        
        #ind_times=np.arange(0, 10)* (times_stat.size-1)/10
        #print ind_times
        #print ind_times.shape
        #times_print=times_stat[ind_times] # every 5 years
        ##times_print=times_stat[::5] # every 5 years
        
        
     
        
        file_in_trw='{inp}/IC_trav_wave.txt'.format(inp=dir_in)
        
        if os.path.isfile(file_in_trw):
        
            data_trw = np.loadtxt(file_in_trw)
            
            sigma=data_trw[4]
            tau=data_trw[1]
        else:
            sigma=3
            tau=1000
        
        c=0
        
        data_viruses_list         =[]
        labels_clusters_opt_list  =[]
        core_samples_mask_opt_list=[]
        
        min_considered_list              =[]
        
        radial_distr_list_binstat              =[]
        radial_distr_matr_list_binstat              =[]
        radius_list_binstat                    =[]
        
        
        
        n_clusters_dbscan_list       =[]
        n_clusters_dbscan_est_list       =[]
        n_clusters_dbscan_opt_kvar_list     =[]
        eps_dbscan_est_list       =[]
        eps_dbscan_opt_kvar_list     =[]
        size_clusters_dbscan_list       =[]
        size_clusters_dbscan_est_list       =[]
        size_clusters_dbscan_opt_kvar_list     =[]
        inter_clusters_dist_dbscan_opt_kvar_list     =[]
        inter_clusters_dist_max_dbscan_opt_kvar_list     =[]
        density_clusters_dbscan_list       =[]
        density_clusters_dbscan_est_list       =[]
        density_clusters_dbscan_opt_kvar_list     =[]
        n_clusters_dbscan_opt_dispersion_list     =[]
        n_clusters_dbscan_opt_dispersion_CV_list     =[]
        n_clusters_dbscan_opt_kCV_list     =[]
        mean_dispersion_dbscan_list       =[]
        mean_dispersion_dbscan_est_list       =[]
        mean_dispersion_dbscan_opt_list     =[]
        mean_dispersion_CV_dbscan_list       =[]
        mean_dispersion_CV_dbscan_est_list       =[]
        mean_dispersion_CV_dbscan_opt_list     =[]
        
        mean_kCV_dbscan_list       =[]
        mean_kCV_dbscan_est_list       =[]
        mean_kCV_dbscan_opt_list     =[]
        mean_kvar_dbscan_list       =[]
        mean_kvar_dbscan_est_list       =[]
        mean_kvar_dbscan_opt_list     =[]
        
        vir_dists_pairwise_mins_max=0
    
           
        times_prep=times
        
        times=[]
        num_vir_tot_list=[]
        fitn_tot_list=[]
        
        
        ################ PERFORM DBSCAN CLUSTERING AND CHOOSE BEST EPSILON BASED ON CLUSTERING SCORE


        for file_in_frames_npz_compr in glob.glob("{inp}/antigenic_space_time_*.dat".format(inp=dir_in_frames, real=real)):
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
                
                
                if int(time)>sec_thr and int(time)<=survival and data_space.ndim > 1:
        
        
                    times.append(time)
            
            
            
                    if len(times)%1 == 0 :
                        #~ print "collecting"
                        gc.collect() 
                        
                        
                        
                        
                    vir_xy = data_space[:,0:2] # 
                    vir_x = data_space[:,0] # 
                    vir_y = data_space[:,1] # 
                    num_vir = data_space[:,3]
                    
                    num_vir_tot_list.append(num_vir)
                    #~ print num_vir.shape
                
                    fitn = data_space[:,4] # 
                    fitn_tot_list.append(fitn)
    #
 
 
            
                    idx_time= np.nonzero(times_prep == time)[0]
                    
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
                                
                    
                    var_parall_cloud[idx_time] = var_parall
                    var_perp_cloud  [idx_time] = var_perp  
                    var_tot_cloud   [idx_time] = var_tot   


                    
                    
                    
                    X=vir_xy
                    #print X.shape
                    
                    mun_cl_coord=10
                    
                    if X.shape[0]<=mun_cl_coord: # NOT ENOUGH TIPES TO CLUSTER
                        
                        labels = np.asarray([0]*X.shape[0])
                        labels_clusters_opt = labels
                        core_samples_mask = np.zeros_like(labels_clusters_opt, dtype=bool)
                        core_samples_mask[:] = True
                        core_samples_mask_opt=core_samples_mask
                        
                        size_clusters = clusters_size(X[labels!=-1], labels[labels!=-1])[0]
                        inter_clusters_dist = inter_clusters_distance(X[labels!=-1], labels[labels!=-1])[0]
                        inter_clusters_dist_max = inter_clusters_distance(X[labels!=-1], labels[labels!=-1])[1]
                        density_clusters = clusters_size(X[labels!=-1], labels[labels!=-1])[1]
    
                        mean_dispersion_avg=mean_dispersion(X, labels)[0]
                        kth_dist_distr_avg=kth_dist_distr(X, labels)[0]
                        kth_dist_distr_var_avg=kth_dist_distr(X, labels)[1]
                        mean_dispersion_renorm_other_avg=mean_dispersion_renorm_other(X, labels)
                                
                       
                        n_clusters_dbscan                   = 1
                        n_clusters_dbscan_est               = 1
                        n_clusters_dbscan_opt_kvar          = 1
                        eps_est                             = 1
                        eps_dbscan_opt_kvar                 = 1
                        size_clusters_dbscan                = size_clusters
                        size_clusters_dbscan_est            = size_clusters
                        size_clusters_dbscan_opt_kvar       = size_clusters
                        inter_clusters_dist_dbscan_opt_kvar = inter_clusters_dist
                        inter_clusters_dist_max_dbscan_opt_kvar = inter_clusters_dist_max
                        density_clusters_dbscan             = density_clusters
                        density_clusters_dbscan_est         = density_clusters
                        density_clusters_dbscan_opt_kvar    = density_clusters
                        n_clusters_dbscan_opt_dispersion    = 1
                        n_clusters_dbscan_opt_dispersion_CV = 1
                        n_clusters_dbscan_opt_kCV           = 1
                        mean_dispersion_dbscan              = mean_dispersion_avg
                        mean_dispersion_dbscan_est          = mean_dispersion_avg
                        mean_dispersion_dbscan_opt          = mean_dispersion_avg
                        mean_dispersion_CV_dbscan           = mean_dispersion_renorm_other_avg
                        mean_dispersion_CV_dbscan_opt       = mean_dispersion_renorm_other_avg
                        mean_dispersion_CV_dbscan_est       = mean_dispersion_renorm_other_avg
                        mean_kCV_dbscan_est                 = kth_dist_distr_avg
                        mean_kCV_dbscan_opt                 = kth_dist_distr_avg
                        mean_kCV_dbscan                     = kth_dist_distr_avg
                        mean_kvar_dbscan_est                = kth_dist_distr_var_avg
                        mean_kvar_dbscan_opt                = kth_dist_distr_var_avg
                        mean_kvar_dbscan                    = kth_dist_distr_var_avg
                        
                        
                        
                        
                    else: # CLUSTER!
                        
                        
                        #~ print "DBSCAN start"
                           
                    
                        if sigma<2:
                            sigma=2
                        
                        #~ print sigma
                        
                        db = DBSCAN(eps=2*sigma, min_samples=mun_cl_coord).fit(X)
                        #db = DBSCAN(eps=140, min_samples=10, metric=dist_nodes_tree(t_full_ext)).fit(alive_viruses[:,None])
                        #core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                        #core_samples_mask[db.core_sample_indices_] = True
                        labels = db.labels_
                        
                        # Number of clusters in labels, ignoring noise if present.
                        n_clusters_dbscan = len(set(labels)) - (1 if -1 in labels else 0)
                        size_clusters_dbscan = clusters_size(X[labels!=-1], labels[labels!=-1])[0]
                        density_clusters_dbscan = clusters_size(X[labels!=-1], labels[labels!=-1])[1]
                        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                        core_samples_mask[db.core_sample_indices_] = True
                        
                        #print "DBSCAN end"
                        
                        
                        mean_dispersion_dbscan=mean_dispersion(X, labels)[0]
                        mean_dispersion_renorm_dbscan=mean_dispersion_renorm(X, labels, core_samples_mask)
                        mean_dispersion_CV_dbscan=mean_dispersion_renorm_other(X, labels)
                        mean_kCV_dbscan=kth_dist_distr(X, labels)[0]
                        mean_kvar_dbscan=kth_dist_distr(X, labels)[1]
                        mean_banfeld_dbscan=banfeld(X, labels)
                        
                        
                        try: # SELECT BEST EPSILON, I TRIED A BUNCH OF SCORES, BUT ULTIMATELY I USE THE VARIANCE OF DISTANCE TO THE 10TH NEIGHBOR
                    
                            range_eps = np.arange(4, 1.9, -0.5) # sigma multiples
                            range_eps=range_eps*sigma
                            
                            mean_dispersion_tot=np.zeros(range_eps.size)
                            mean_dispersion_sum_tot=np.zeros(range_eps.size)
                            mean_dispersion_renorm_tot=np.zeros(range_eps.size)
                            kth_dist_distr_tot=np.zeros(range_eps.size)
                            kth_dist_distr_var_tot=np.zeros(range_eps.size)
                            kth_dist_distr_var_sum_tot=np.zeros(range_eps.size)
                            mean_dispersion_renorm_other_tot=np.zeros(range_eps.size)
                            mean_banfeld_tot=np.zeros(range_eps.size)
                            #mean_banfeld_unif_tot=np.zeros(range_eps.size)
                            n_clusters_tot=np.zeros(range_eps.size)
                            size_clusters_tot=np.zeros(range_eps.size)
                            inter_clusters_dist_tot=np.zeros(range_eps.size)
                            inter_clusters_dist_max_tot=np.zeros(range_eps.size)
                            density_clusters_tot=np.zeros(range_eps.size)
                            
                            
                            
                            labels_clusters_opt  =np.zeros((range_eps.size, X.shape[0]), dtype=int)
                            core_samples_mask_opt=np.zeros((range_eps.size, X.shape[0]), dtype=bool)
                            
                            #~ print "DBSCAN opt" # faster than kmeans
                            for j, eps in enumerate(range_eps):
                            
                                # Initialize the clusterer with n_clusters value and a random generator
                                # seed of 10 for reproducibility.
                                db = DBSCAN(eps=eps, min_samples=mun_cl_coord).fit(X)
                                labels = db.labels_
                                
                                core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                                core_samples_mask[db.core_sample_indices_] = True
                                
                                # clusters
                                n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
                                size_clusters = clusters_size(X[labels!=-1], labels[labels!=-1])[0]
                                inter_clusters_dist = inter_clusters_distance(X[labels!=-1], labels[labels!=-1])[0]
                                inter_clusters_dist_max = inter_clusters_distance(X[labels!=-1], labels[labels!=-1])[1]
                                density_clusters = clusters_size(X[labels!=-1], labels[labels!=-1])[1]
                    
                                #if n_clusters>1:
                                #    print n_clusters, eps, eps/sigma
                    
                    
                                mean_dispersion_avg=mean_dispersion(X, labels)[0]
                                mean_dispersion_sum_avg=mean_dispersion(X, labels)[1]
                                mean_dispersion_renorm_avg=mean_dispersion_renorm(X, labels,core_samples_mask)
                                kth_dist_distr_avg=kth_dist_distr(X, labels)[0]
                                kth_dist_distr_var_avg=kth_dist_distr(X, labels)[1]
                                kth_dist_distr_var_sum_avg=kth_dist_distr(X, labels)[2]
                                mean_dispersion_renorm_other_avg=mean_dispersion_renorm_other(X, labels)
                                mean_banfeld_avg=banfeld(X, labels)
                                #mean_banfeld_unif_avg=banfeld_unif(X, labels)
                                
                                compar=max(labels.size/5., 10)
                                
                                if labels[labels==-1].size < compar:
                        
                                    n_clusters_tot[j]                   =n_clusters
                                    size_clusters_tot[j]                =size_clusters
                                    inter_clusters_dist_tot[j]          =inter_clusters_dist
                                    inter_clusters_dist_max_tot[j]          =inter_clusters_dist_max
                                    density_clusters_tot[j]             =density_clusters
                                    mean_dispersion_tot[j]              =mean_dispersion_avg
                                    mean_dispersion_sum_tot[j]          =mean_dispersion_sum_avg
                                    mean_dispersion_renorm_tot[j]       =mean_dispersion_renorm_avg
                                    kth_dist_distr_tot[j]               =kth_dist_distr_avg
                                    kth_dist_distr_var_tot[j]           =kth_dist_distr_var_avg
                                    kth_dist_distr_var_sum_tot[j]       =kth_dist_distr_var_sum_avg
                                    mean_dispersion_renorm_other_tot[j] =mean_dispersion_renorm_other_avg
                                    mean_banfeld_tot[j]                 =mean_banfeld_avg
                                    
                                    labels_clusters_opt[j,:]=labels
                                    core_samples_mask_opt[j,:]=core_samples_mask
                                
                                else:
                                    print eps, labels.size, labels[labels==-1].size
                                
                                
                                
                            
                            
                            range_eps=range_eps[n_clusters_tot>0]
                            
                            labels_clusters_opt  =  labels_clusters_opt[n_clusters_tot>0,:]
                            core_samples_mask_opt=core_samples_mask_opt[n_clusters_tot>0,:]
                            
                            mean_dispersion_tot=mean_dispersion_tot[n_clusters_tot>0]
                            mean_dispersion_sum_tot=mean_dispersion_sum_tot[n_clusters_tot>0]
                            mean_dispersion_renorm_tot=mean_dispersion_renorm_tot[n_clusters_tot>0]
                            kth_dist_distr_tot=kth_dist_distr_tot[n_clusters_tot>0]
                            kth_dist_distr_var_tot=kth_dist_distr_var_tot[n_clusters_tot>0]
                            kth_dist_distr_var_sum_tot=kth_dist_distr_var_sum_tot[n_clusters_tot>0]
                            mean_dispersion_renorm_other_tot=mean_dispersion_renorm_other_tot[n_clusters_tot>0]
                            mean_banfeld_tot=mean_banfeld_tot[n_clusters_tot>0]
                            #mean_banfeld_unif_tot=mean_banfeld_unif_tot[n_clusters_tot>0]
                            size_clusters_tot     =size_clusters_tot[n_clusters_tot>0]
                            inter_clusters_dist_tot     =inter_clusters_dist_tot[n_clusters_tot>0]
                            inter_clusters_dist_max_tot =inter_clusters_dist_max_tot[n_clusters_tot>0]
                            density_clusters_tot     =density_clusters_tot[n_clusters_tot>0]
                            n_clusters_tot     =n_clusters_tot[n_clusters_tot>0]
                    
                            
                            disp_acceleration=mean_dispersion_tot[1:] - mean_dispersion_tot[:-1]#slope
                            disp_acceleration=disp_acceleration[1:] - disp_acceleration[:-1] #acceleration
                            
                            indmax_dispersion=np.argmax(disp_acceleration) + 1
                            indmax_dispersion_CV=np.argmin(mean_dispersion_renorm_other_tot)
                            indmax_kCV=np.argmin(kth_dist_distr_tot) 
                            indmax_kvar=np.argmin(kth_dist_distr_var_tot) 
                            
                            mean_dispersion_dbscan_opt=mean_dispersion_tot[indmax_dispersion]
                            mean_dispersion_CV_dbscan_opt=mean_dispersion_renorm_other_tot[indmax_dispersion_CV]
                            #mean_banfeld_dbscan_opt=mean_banfeld_tot[indmax_banfeld]
                            mean_kCV_dbscan_opt=kth_dist_distr_tot[indmax_kCV]
                            mean_kvar_dbscan_opt=kth_dist_distr_var_tot[indmax_kvar]
                            n_clusters_dbscan_opt_dispersion=n_clusters_tot[indmax_dispersion]
                            n_clusters_dbscan_opt_dispersion_CV=n_clusters_tot[indmax_dispersion_CV]
                            #n_clusters_dbscan_opt_banfeld=n_clusters_tot[indmax_banfeld]
                            n_clusters_dbscan_opt_kCV=n_clusters_tot[indmax_kCV]
                            n_clusters_dbscan_opt_kvar=n_clusters_tot[indmax_kvar]
                            eps_dbscan_opt_kvar=range_eps[indmax_kvar]
                            size_clusters_dbscan_opt_kvar=size_clusters_tot[indmax_kvar]
                            inter_clusters_dist_dbscan_opt_kvar=inter_clusters_dist_tot[indmax_kvar]
                            inter_clusters_dist_max_dbscan_opt_kvar=inter_clusters_dist_max_tot[indmax_kvar]
                            density_clusters_dbscan_opt_kvar=density_clusters_tot[indmax_kvar]
                            #~ print "DBSCAN opt end, eps ", range_eps[indmax_dispersion], range_eps[indmax_kCV], range_eps[indmax_kvar], real, time
                            #~ print n_clusters_tot
                            #~ print mean_dispersion_tot
                            #~ print kth_dist_distr_var_tot
                            
                            
                            
                            
                
                
                            labels_clusters_opt=labels_clusters_opt[indmax_kvar,:]
                            core_samples_mask_opt=core_samples_mask_opt[indmax_kvar,:] 
                            
                            #~ print labels_clusters_opt.shape, core_samples_mask_opt.shape,  X.shape
                            
                            
                            
                
                    
                    
                            
                            numX=X.shape[0]
                            kX=numX
                            if numX>10:
                                kX=10
                        
                            
                            
                            nbrs = NearestNeighbors(n_neighbors=kX).fit(X)
                            distances, indices = nbrs.kneighbors(X)
                            #print distances.shape
                            kth_dist=distances[:,-1]
                            kth_dist[::-1].sort()
                            window=3
                            kth_dist=np.convolve(kth_dist, np.ones((window,))/window, mode='same')
                            #print kth_dist
                            kth_dist_acceleration=kth_dist[5:] - kth_dist[:-5]#slope
                            kth_dist_acceleration=kth_dist_acceleration[1:] - kth_dist_acceleration[:-1] #acceleration
                            
                    
                            indmax_kth_dist=np.argmax(kth_dist_acceleration) + 5
                            eps_est=kth_dist[indmax_kth_dist]
                            if eps_est<2*sigma:
                                eps_est=2*sigma
                    
                            
                        
                        except ValueError:
                            print("ERROR discarding noise in DBSCAN")
                            
                            mean_dispersion_dbscan_opt=0
                            mean_dispersion_CV_dbscan_opt=0
                            mean_kCV_dbscan_opt=0
                            mean_kvar_dbscan_opt=0
                            n_clusters_dbscan_opt_dispersion=0
                            n_clusters_dbscan_opt_dispersion_CV=0
                            n_clusters_dbscan_opt_kCV=0
                            n_clusters_dbscan_opt_kvar=0
                            eps_dbscan_opt_kvar=0
                            size_clusters_dbscan_opt_kvar=0
                            inter_clusters_dist_dbscan_opt_kvar=0
                            inter_clusters_dist_max_dbscan_opt_kvar=0
                            density_clusters_dbscan_opt_kvar=0
                            eps_est=2*sigma
                            
                            
                            labels_clusters_opt  =np.zeros((X.shape[0]), dtype=int)		
                            core_samples_mask_opt=np.zeros((X.shape[0]), dtype=bool)		
                        
                        #~ print eps_est
            
                        # ALSO TRY THE "ELBOW ESTIMATION" FOR EPSILON
                        
                        db = DBSCAN(eps=eps_est, min_samples=mun_cl_coord).fit(X)
                        #db = DBSCAN(eps=140, min_samples=10, metric=dist_nodes_tree(t_full_ext)).fit(alive_viruses[:,None])
                        #core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                        #core_samples_mask[db.core_sample_indices_] = True
                        labels = db.labels_
                        
                        # Number of clusters in labels, ignoring noise if present.
                        n_clusters_dbscan_est = len(set(labels)) - (1 if -1 in labels else 0)
                        size_clusters_dbscan_est = clusters_size(X[labels!=-1], labels[labels!=-1])[0]
                        density_clusters_dbscan_est = clusters_size(X[labels!=-1], labels[labels!=-1])[1]
                        
                        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                        core_samples_mask[db.core_sample_indices_] = True
                        
                        #print "DBSCAN end"
                        
                        
                        mean_dispersion_dbscan_est=mean_dispersion(X, labels)[0]
                        mean_dispersion_renorm_dbscan_est=mean_dispersion_renorm(X, labels, core_samples_mask)
                        mean_dispersion_CV_dbscan_est=mean_dispersion_renorm_other(X, labels)
                        #print mean_dispersion_CV_dbscan_est
                        mean_kCV_dbscan_est=kth_dist_distr(X, labels)[0]
                        mean_kvar_dbscan_est=kth_dist_distr(X, labels)[1]
                        mean_banfeld_dbscan_est=banfeld(X, labels)
                        
                        
                        
                    
                    #print "clusters ", n_clusters_dbscan, n_clusters_dbscan_bigeps, n_clusters_hier, n_clusters_hier_opt
                    
                    n_clusters_dbscan_list.append(n_clusters_dbscan)
                    n_clusters_dbscan_est_list.append(n_clusters_dbscan_est)
                    n_clusters_dbscan_opt_kvar_list.append(n_clusters_dbscan_opt_kvar)
                    eps_dbscan_est_list.append(eps_est)
                    eps_dbscan_opt_kvar_list.append(eps_dbscan_opt_kvar)
                    size_clusters_dbscan_list.append(size_clusters_dbscan)
                    size_clusters_dbscan_est_list.append(size_clusters_dbscan_est)
                    size_clusters_dbscan_opt_kvar_list.append(size_clusters_dbscan_opt_kvar)
                    inter_clusters_dist_dbscan_opt_kvar_list.append(inter_clusters_dist_dbscan_opt_kvar)
                    inter_clusters_dist_max_dbscan_opt_kvar_list.append(inter_clusters_dist_max_dbscan_opt_kvar)
                    density_clusters_dbscan_list.append(density_clusters_dbscan)
                    density_clusters_dbscan_est_list.append(density_clusters_dbscan_est)
                    density_clusters_dbscan_opt_kvar_list.append(density_clusters_dbscan_opt_kvar)
                    n_clusters_dbscan_opt_dispersion_list.append(n_clusters_dbscan_opt_dispersion)
                    n_clusters_dbscan_opt_dispersion_CV_list.append(n_clusters_dbscan_opt_dispersion_CV)
                    n_clusters_dbscan_opt_kCV_list.append(n_clusters_dbscan_opt_kCV)
                    mean_dispersion_dbscan_list.append(mean_dispersion_dbscan)
                    mean_dispersion_dbscan_est_list.append(mean_dispersion_dbscan_est)
                    mean_dispersion_dbscan_opt_list.append(mean_dispersion_dbscan_opt)
                    mean_dispersion_CV_dbscan_list.append(mean_dispersion_CV_dbscan)
                    mean_dispersion_CV_dbscan_opt_list.append(mean_dispersion_CV_dbscan_opt)
                    mean_dispersion_CV_dbscan_est_list.append(mean_dispersion_CV_dbscan_est)
                    mean_kCV_dbscan_est_list.append(mean_kCV_dbscan_est)
                    mean_kCV_dbscan_opt_list.append(mean_kCV_dbscan_opt)
                    mean_kCV_dbscan_list.append(mean_kCV_dbscan)
                    mean_kvar_dbscan_est_list.append(mean_kvar_dbscan_est)
                    mean_kvar_dbscan_opt_list.append(mean_kvar_dbscan_opt)
                    mean_kvar_dbscan_list.append(mean_kvar_dbscan)
                    
                    
                    
                    
                    
                    
                    data_viruses_list.append(X)
                    labels_clusters_opt_list.append(labels_clusters_opt  )
                    core_samples_mask_opt_list.append(core_samples_mask_opt)
                    
                    
                    gc.collect()
                    
        
        print len(times), len(n_clusters_dbscan_opt_dispersion_CV_list), len(n_clusters_dbscan_opt_kCV_list), len(n_clusters_dbscan_opt_kvar_list), len(mean_dispersion_CV_dbscan_list), len(mean_dispersion_CV_dbscan_est_list), len(mean_dispersion_CV_dbscan_opt_list), len(mean_kCV_dbscan_list), len(mean_kCV_dbscan_est_list), len(mean_kCV_dbscan_opt_list), len(mean_kvar_dbscan_list), len(mean_kvar_dbscan_est_list), len(mean_kvar_dbscan_opt_list)
        
        print len(data_viruses_list), len(labels_clusters_opt_list), len(core_samples_mask_opt_list)
        
        
        times, n_clusters_dbscan_list, n_clusters_dbscan_est_list, n_clusters_dbscan_opt_kvar_list, eps_dbscan_est_list, eps_dbscan_opt_kvar_list, size_clusters_dbscan_list, size_clusters_dbscan_est_list, size_clusters_dbscan_opt_kvar_list, inter_clusters_dist_dbscan_opt_kvar_list, inter_clusters_dist_max_dbscan_opt_kvar_list, density_clusters_dbscan_list, density_clusters_dbscan_est_list, density_clusters_dbscan_opt_kvar_list, n_clusters_dbscan_opt_dispersion_list, n_clusters_dbscan_opt_dispersion_CV_list, n_clusters_dbscan_opt_kCV_list, mean_dispersion_dbscan_list, mean_dispersion_dbscan_est_list, mean_dispersion_dbscan_opt_list, mean_dispersion_CV_dbscan_list, mean_dispersion_CV_dbscan_est_list, mean_dispersion_CV_dbscan_opt_list, mean_kCV_dbscan_list, mean_kCV_dbscan_est_list, mean_kCV_dbscan_opt_list, mean_kvar_dbscan_list, mean_kvar_dbscan_est_list, mean_kvar_dbscan_opt_list, data_viruses_list, labels_clusters_opt_list, core_samples_mask_opt_list, num_vir_tot_list, fitn_tot_list= zip(*sorted(zip(times, n_clusters_dbscan_list, n_clusters_dbscan_est_list, n_clusters_dbscan_opt_kvar_list, eps_dbscan_est_list, eps_dbscan_opt_kvar_list, size_clusters_dbscan_list, size_clusters_dbscan_est_list, size_clusters_dbscan_opt_kvar_list, inter_clusters_dist_dbscan_opt_kvar_list, inter_clusters_dist_max_dbscan_opt_kvar_list, density_clusters_dbscan_list, density_clusters_dbscan_est_list, density_clusters_dbscan_opt_kvar_list, n_clusters_dbscan_opt_dispersion_list, n_clusters_dbscan_opt_dispersion_CV_list, n_clusters_dbscan_opt_kCV_list, mean_dispersion_dbscan_list, mean_dispersion_dbscan_est_list, mean_dispersion_dbscan_opt_list, mean_dispersion_CV_dbscan_list, mean_dispersion_CV_dbscan_est_list, mean_dispersion_CV_dbscan_opt_list, mean_kCV_dbscan_list, mean_kCV_dbscan_est_list, mean_kCV_dbscan_opt_list, mean_kvar_dbscan_list, mean_kvar_dbscan_est_list, mean_kvar_dbscan_opt_list, data_viruses_list, labels_clusters_opt_list, core_samples_mask_opt_list, num_vir_tot_list, fitn_tot_list)))
        
        
        
        
        times= np.asarray(times)
        #num_vir_tot_list= np.asarray(num_vir_tot_list)
        
        
        if times.size > 1:
            
            n_clusters_dbscan_list       =np.asarray(n_clusters_dbscan_list       )
            n_clusters_dbscan_est_list       =np.asarray(n_clusters_dbscan_est_list       )
            n_clusters_dbscan_opt_kvar_list     =np.asarray(n_clusters_dbscan_opt_kvar_list     )
            eps_dbscan_est_list       =np.asarray(eps_dbscan_est_list       )
            eps_dbscan_opt_kvar_list     =np.asarray(eps_dbscan_opt_kvar_list     )
            size_clusters_dbscan_list       =np.asarray(size_clusters_dbscan_list       )
            size_clusters_dbscan_est_list       =np.asarray(size_clusters_dbscan_est_list       )
            size_clusters_dbscan_opt_kvar_list     =np.asarray(size_clusters_dbscan_opt_kvar_list     )
            inter_clusters_dist_dbscan_opt_kvar_list     =np.asarray(inter_clusters_dist_dbscan_opt_kvar_list     )
            inter_clusters_dist_max_dbscan_opt_kvar_list =np.asarray(inter_clusters_dist_max_dbscan_opt_kvar_list     )
            density_clusters_dbscan_list       =np.asarray(density_clusters_dbscan_list       )
            density_clusters_dbscan_est_list       =np.asarray(density_clusters_dbscan_est_list       )
            density_clusters_dbscan_opt_kvar_list     =np.asarray(density_clusters_dbscan_opt_kvar_list     )
            n_clusters_dbscan_opt_dispersion_list     =np.asarray(n_clusters_dbscan_opt_dispersion_list     )
            n_clusters_dbscan_opt_dispersion_CV_list     =np.asarray(n_clusters_dbscan_opt_dispersion_CV_list     )
            n_clusters_dbscan_opt_kCV_list     =np.asarray(n_clusters_dbscan_opt_kCV_list     )
            mean_dispersion_dbscan_list       =np.asarray(mean_dispersion_dbscan_list       )
            mean_dispersion_dbscan_est_list       =np.asarray(mean_dispersion_dbscan_est_list       )
            mean_dispersion_dbscan_opt_list     =np.asarray(mean_dispersion_dbscan_opt_list     )
            
            mean_dispersion_CV_dbscan_list       =np.asarray(mean_dispersion_CV_dbscan_list       )
            mean_dispersion_CV_dbscan_est_list       =np.asarray(mean_dispersion_CV_dbscan_est_list       )
            mean_dispersion_CV_dbscan_opt_list     =np.asarray(mean_dispersion_CV_dbscan_opt_list     )
            
            mean_kCV_dbscan_list       =np.asarray(mean_kCV_dbscan_list       )
            mean_kCV_dbscan_est_list       =np.asarray(mean_kCV_dbscan_est_list       )
            mean_kCV_dbscan_opt_list     =np.asarray(mean_kCV_dbscan_opt_list     )
            
            mean_kvar_dbscan_list       =np.asarray(mean_kvar_dbscan_list       )
            mean_kvar_dbscan_est_list       =np.asarray(mean_kvar_dbscan_est_list       )
            mean_kvar_dbscan_opt_list     =np.asarray(mean_kvar_dbscan_opt_list     )
                       
            
            #CLUSTER TRACKING! USINH OPT KVAR
            print "CLUSTER TRACKING! USING OPT KVAR"
    
            
            
            vel_clust_split_alltimes        =[]
            var_parall_clust_split_alltimes        =[]
            var_perp_clust_split_alltimes        =[]
            var_tot_clust_split_alltimes        =[]
            
            times_diff_fity_split_alltimes_allclust               =[]
            diff_fity_split_alltimes_allclust                     =[]
            diff_fitx_split_alltimes_allclust                     =[]
            ang_past_dir_shortmem_split_alltimes_allclust         =[]
            ang_past_dir_split_alltimes_allclust                  =[]
            
            times_sharp_turn_split_alltimes_allclust        =[]
            times_grad_split_alltimes_allclust        =[]
            grad_fittest_split_alltimes_allclust        =[]
            shortmem_gradx_bulk_split_alltimes_allclust        =[]
            shortmem_grady_fittest_split_alltimes_allclust        =[]
            shortmem_grady_bulk_split_alltimes_allclust        =[]
            gradx_bulk_split_alltimes_allclust        =[]
            grady_fittest_split_alltimes_allclust        =[]
            grady_bulk_split_alltimes_allclust        =[]
            
            vel_clust_split_alltimes_allclust        =[]
            var_parall_clust_split_alltimes_allclust        =[]
            var_perp_clust_split_alltimes_allclust        =[]
            var_tot_clust_split_alltimes_allclust        =[]
            
            size_clusters_split        =[]
            inter_clusters_dist_split  =[]
            inter_clusters_dist_max_split  =[]
            density_clusters_split     =[]
            n_clusters_split           =[]
            num_vir_clusters_split     =[]
            num_coords_clusters_split  =[]
            
            data_track_tot     =[]
            
            num_ext_tot  =np.zeros(times.size)
            num_split_tot=np.zeros(times.size)
            
            num_ext_per_cl_tot  =np.zeros(times.size)
            num_split_per_cl_tot=np.zeros(times.size)
            
            last_centr_ID=0
            
            xs_all_centr_all_time =[] # list of the lists of all xs a centroid had for all times, in ID order
            ys_all_centr_all_time =[] # list of the lists of all ys a centroid had for all times, in ID order
            ts_all_centr_all_time =[] # list of the lists of all ys a centroid had for all times, in ID order
    
            maxfit_wrt_bulk_newbasis_cumulpath_x=0
            maxfit_wrt_bulk_newbasis_cumulpath_y=0
            
            
            #####  REFINE CLUSTERING BASED ON CLUSTERS DISTANCE WITH RESPECT TO CLUSTER SIZES, THEN TRACK THE LINEAGES
    
    
            for i, time in enumerate(times):
                
                
                X=data_viruses_list[i]
                num_vir=num_vir_tot_list[i]
                fitn=fitn_tot_list[i]
	
                gradtot=np.array([0])
                
                
                
                
                labels_clusters_opt = labels_clusters_opt_list[i]
                
                core_samples_mask_opt = core_samples_mask_opt_list[i]
                
                #~ print i, time, times.shape
                #~ print X.shape, labels_clusters_opt.shape, core_samples_mask_opt.shape
                
                
                
                 
                recluster_split_data=recluster_split(X, labels_clusters_opt, core_samples_mask_opt, num_vir) #####  REFINE CLUSTERING BASED ON CLUSTERS DISTANCE WITH RESPECT TO CLUSTER SIZES
                
                labels_split=recluster_split_data[0]
                centroids_split=np.asarray(recluster_split_data[1])
                nums_split =np.asarray(recluster_split_data[2])
                sizes_split=np.asarray(recluster_split_data[3])
                clust_labels_split=np.asarray(recluster_split_data[4])
                nums_coords_split=np.asarray(recluster_split_data[5])
                
                vel_clust_split=[]
            
                var_parall_clust_split = []
                var_perp_clust_split   = []
                var_tot_clust_split    = []
        
             
                #~ print labels_split.shape, centroids_split.shape, nums_split.shape, sizes_split.shape
                
                IDs_split           =np.zeros(nums_split.size) # ID of the centroid
                t_birth_split       =np.zeros(nums_split.size) # time of appearence of the centroid
                IDs_fath_split      =np.zeros(nums_split.size) # ID of the centroid that originated him
                num_vir_fath_split  =np.zeros(nums_split.size) # 
                num_coords_fath_split  =np.zeros(nums_split.size)
                sizes_fath_split    =np.zeros(nums_split.size) # 
                centroids_fath_split=np.zeros(centroids_split.shape) # 
                
        
                num_ext=0
                num_split=0
                
                    
                if i==0:
                    for i_c, centr in enumerate(centroids_split):
                        
                        #~ print i_c, centr
                        #~ print centr.shape
                        
                        last_centr_ID+=1
                        IDs_split[i_c]=last_centr_ID
                        t_birth_split[i_c]=time
                        
                        # fake father data
                        IDs_fath_split[i_c]=-1
                        num_vir_fath_split[i_c]=-1
                        num_coords_fath_split[i_c]=-1
                        sizes_fath_split[i_c]=-1
                        centroids_fath_split[i_c]=np.array([0.,0.])
                        
                        state=-1
                        
                        data_new = np.array([time, state, last_centr_ID, t_birth_split[i_c], centr[0], centr[1], nums_split[i_c], sizes_split[i_c], -1, 0., 0., -1, -1, 0, -1, 0., 0., 0., nums_coords_split[i_c], -1, 0, 0])
                        
                        data_track_tot.append(data_new)
                        
                        ts_all_centr_all_time.append([time])
                        xs_all_centr_all_time.append([centr[0]])
                        ys_all_centr_all_time.append([centr[1]])
                        vel_clust_split.append(0.)
                        
                        var_parall_clust_split.append(0.)
                        var_perp_clust_split.append(0.  )
                        var_tot_clust_split.append(0.   )
            
            
                        vel_clust_split_alltimes_allclust.append(0.)
                        
                        var_parall_clust_split_alltimes_allclust.append(0.)
                        var_perp_clust_split_alltimes_allclust.append(0.  )
                        var_tot_clust_split_alltimes_allclust.append(0.   )
            
                
                    
                else:
                    idxs_closest_prec=np.zeros(nums_split.size) # 
                    for i_c, centr in enumerate(centroids_split): #find closest precedent centroid
                        
                        dists_prec=np.linalg.norm(centroids_split_prec - centr[None,:], axis=1)
                        
                        idx_closest_prec=np.argmin(dists_prec)
                        #if not idx_closest_prec.isscalar():
                        #    idx_closest_prec=idx_closest_prec[0]
                        idxs_closest_prec[i_c]=idx_closest_prec
            
                        
        
                    
                    for i_c_prec, centr_prec in enumerate(centroids_split_prec): 
                        idx_sons=np.nonzero(idxs_closest_prec==i_c_prec)[0]
                        
                        #~ print idx_sons
                        
                        if idx_sons.shape[0]==0: # extinction
                            
                            ID_ext          =IDs_split_prec[i_c_prec]
                            num_vir_ext     =num_vir_split_prec[i_c_prec]
                            num_coords_ext     =num_coords_split_prec[i_c_prec]
                            size_ext        =sizes_split_prec[i_c_prec]
                                    
                            ID_fath_ext     =IDs_fath_split_prec[i_c_prec]
                            num_vir_fath_ext=num_vir_fath_split_prec[i_c_prec]
                            num_coords_fath_ext=num_coords_fath_split_prec[i_c_prec]
                            size_fath_ext   =sizes_fath_split_prec[i_c_prec]
                            centr_fath_ext  =centroids_fath_split_prec[i_c_prec]
                            
                            t_birth_ext   =  t_birth_split_prec[i_c_prec]
                            
                            state=0
            
                            data_ext = np.array([time, state, ID_ext, t_birth_ext, centr_prec[0], centr_prec[1], num_vir_ext, size_ext, ID_fath_ext, centr_fath_ext[0], centr_fath_ext[1], num_vir_fath_ext, size_fath_ext, 0, -1, 0., 0., 0., num_coords_ext, num_coords_fath_ext, 0, 0])
                            
                            data_track_tot.append(data_ext)
                            
                            num_ext+=1
                            
                            
                            
                        elif idx_sons.shape[0]==1: # cluster goes on
                            
                            idx_son=idx_sons[0]
            
                            ID_son=IDs_split_prec[i_c_prec] # just tracks, so ID stays the same
                            t_birth_son=t_birth_split_prec[i_c_prec] # just tracks, so t birth stays the same
                            
                            num_coords_son=nums_coords_split[idx_son]
                            num_vir_son=nums_split[idx_son]
                            size_son=sizes_split[idx_son]
                            centr_son=centroids_split[idx_son]
                            label_son=clust_labels_split[idx_son]
                            
                            ID_fath_son=IDs_fath_split_prec[i_c_prec] # father data is the same as before
                            num_coords_fath_son=num_coords_fath_split_prec[i_c_prec]
                            num_vir_fath_son=num_vir_fath_split_prec[i_c_prec]
                            size_fath_son =sizes_fath_split_prec[i_c_prec]
                            centr_fath_son=centroids_fath_split_prec[i_c_prec]
                            
                            state=1
                            
                            vel_clust = np.sqrt((centr_prec[0] - centr_son[0])**2 + (centr_prec[1] - centr_son[1])**2)/(time - times[i-1])
                            vel_clust_split.append(vel_clust)
                            vel_clust_split_alltimes_allclust.append(vel_clust)
                            
                            
                            clust_dir=np.array([centr_son[0] - centr_prec[0] , centr_son[1] - centr_prec[1]])/np.sqrt((centr_prec[0] - centr_son[0])**2 + (centr_prec[1] - centr_son[1])**2)
            
                            clust_label_prec =clust_labels_split_prec[i_c_prec]
                            
                            #~ print "tracking check "
                            #~ print clust_label_prec, num_vir_split_prec[i_c_prec], (labels_split_prec==clust_label_prec).sum()
                            
                            
                            
                            
                            vars_dist = clusters_dist_vars(data_viruses_list[i-1], labels_split_prec, clust_label_prec, clust_dir, num_vir_tot_list[i-1])
                            
                            var_parall= vars_dist[0]
                            var_perp  = vars_dist[1]
                            var_tot   = vars_dist[2]
                            
                            var_parall_clust_split_alltimes_allclust.append(var_parall)
                            var_perp_clust_split_alltimes_allclust.append(var_perp  )
                            var_tot_clust_split_alltimes_allclust.append(var_tot   )
                            
                            
                            var_parall_clust_split.append(var_parall)
                            var_perp_clust_split.append(var_perp  )
                            var_tot_clust_split.append(var_tot   )

                            data_son = np.array([time, state, ID_son, t_birth_son, centr_son[0], centr_son[1], num_vir_son, size_son, ID_fath_son, centr_fath_son[0], centr_fath_son[1], num_vir_fath_son, size_fath_son, vel_clust, clust_label_prec, var_parall, var_perp, var_tot, num_coords_son, num_coords_fath_son, clust_dir[0], clust_dir[1]])
                            
                            data_track_tot.append(data_son)
                            
                            
                            
                            data_track_tot_tmp=np.asarray(data_track_tot)
                            #~ print data_track_tot_tmp.shape
                            past_times = data_track_tot_tmp[data_track_tot_tmp[:,2] == ID_son, 0] - time
                            past_clust_dirs_x = data_track_tot_tmp[data_track_tot_tmp[:,2] == ID_son, -2]
                            past_clust_dirs_y = data_track_tot_tmp[data_track_tot_tmp[:,2] == ID_son, -1]
                            mask_past_dir=(past_clust_dirs_x!=0) | (past_clust_dirs_y!=0)
                            past_clust_dirs_x = past_clust_dirs_x[mask_past_dir]
                            past_clust_dirs_y = past_clust_dirs_y[mask_past_dir]
                            past_times = past_times[mask_past_dir]
                            
                            
                            
            
                            
                            if len(xs_all_centr_all_time[int(ID_son) - 1])==1:
                                ts_all_centr_all_time[int(ID_son) - 1]=[ts_all_centr_all_time[int(ID_son) - 1][0], time]
                                xs_all_centr_all_time[int(ID_son) - 1]=[xs_all_centr_all_time[int(ID_son) - 1][0], centr_son[0]]
                                ys_all_centr_all_time[int(ID_son) - 1]=[ys_all_centr_all_time[int(ID_son) - 1][0], centr_son[1]]
                            else:
                                ts_all_centr_all_time[int(ID_son) - 1].extend([time])
                                xs_all_centr_all_time[int(ID_son) - 1].extend([centr_son[0]])
                                ys_all_centr_all_time[int(ID_son) - 1].extend([centr_son[1]])
                
                
                            IDs_split[idx_son]=ID_son
                            t_birth_split[idx_son]=t_birth_son
                            
                            # father data, that stay the same
                            IDs_fath_split[idx_son]      =ID_fath_son
                            num_coords_fath_split[idx_son]  =num_coords_fath_son
                            num_vir_fath_split[idx_son]  =num_vir_fath_son
                            sizes_fath_split[idx_son]    =size_fath_son
                            centroids_fath_split[idx_son]=centr_fath_son
                            
                        
                        else: # splitting
                        
            
                            for idx_son in idx_sons: # all the new clusters that have this prec as closest
                 
                                last_centr_ID+=1
                 
                    
                                ID_son=last_centr_ID
                                
                                t_birth_son=time
                                
                                num_coords_son=nums_coords_split[idx_son]
                                num_vir_son=nums_split[idx_son]
                                size_son=sizes_split[idx_son]
                                centr_son=centroids_split[idx_son]
                                
                                ID_fath_son=IDs_split_prec[i_c_prec] # father data is from current  (prec) cluster
                                num_coords_fath_son=num_coords_split_prec[i_c_prec]
                                num_vir_fath_son=num_vir_split_prec[i_c_prec]
                                size_fath_son =sizes_split_prec[i_c_prec]
                                centr_fath_son=centr_prec
                                
                                vel_clust = np.sqrt((centr_prec[0] - centr_son[0])**2 + (centr_prec[1] - centr_son[1])**2)/(time - times[i-1])
                                
                                clust_dir=np.array([centr_son[0] - centr_prec[0] , centr_son[1] - centr_prec[1]])/np.sqrt((centr_prec[0] - centr_son[0])**2 + (centr_prec[1] - centr_son[1])**2)
                
                                clust_label_prec =clust_labels_split_prec[i_c_prec]
                                
                                #~ print "splitting check "
                                #~ print clust_label_prec, num_vir_fath_son, (labels_split_prec==clust_label_prec).sum()
                
                                
                                
                                vars_dist = clusters_dist_vars(data_viruses_list[i-1], labels_split_prec, clust_label_prec, clust_dir, num_vir_tot_list[i-1])
                                
                                var_parall= vars_dist[0]
                                var_perp  = vars_dist[1]
                                var_tot   = vars_dist[2]
                                            
                                
                                
                                vel_clust_split_alltimes_allclust.append(vel_clust)
                    
                                var_parall_clust_split_alltimes_allclust.append(var_parall)
                                var_perp_clust_split_alltimes_allclust.append(var_perp  )
                                var_tot_clust_split_alltimes_allclust.append(var_tot   )
                                    
                                
                                vel_clust_split.append(vel_clust)
                    
                                var_parall_clust_split.append(var_parall)
                                var_perp_clust_split.append(var_perp  )
                                var_tot_clust_split.append(var_tot   )
                                
                                
                                
                                
                                state=idx_sons.shape[0]
                    
                                data_son = np.array([time, state, ID_son, t_birth_son,  centr_son[0], centr_son[1], num_vir_son, size_son, ID_fath_son, centr_fath_son[0], centr_fath_son[1], num_vir_fath_son, size_fath_son, vel_clust, clust_label_prec, var_parall, var_perp, var_tot, num_coords_son, num_coords_fath_son, clust_dir[0], clust_dir[1]])
                                
                                data_track_tot.append(data_son)
                    
                                
                                #xs_all_centr_all_time.append([centr_son[0]])
                                #ys_all_centr_all_time.append([centr_son[1]])
                                ts_all_centr_all_time.append([time, time])
                                xs_all_centr_all_time.append([centr_fath_son[0], centr_son[0]])
                                ys_all_centr_all_time.append([centr_fath_son[1], centr_son[1]])
                        
                    
                                IDs_split[idx_son]=ID_son
                                t_birth_split[idx_son]=t_birth_son
                                
                                # father data, that stay the same
                                IDs_fath_split[idx_son]      =ID_fath_son
                                num_coords_fath_split[idx_son]  =num_coords_fath_son
                                num_vir_fath_split[idx_son]  =num_vir_fath_son
                                sizes_fath_split[idx_son]    =size_fath_son
                                centroids_fath_split[idx_son]=centr_fath_son
                    
                
                                num_split+=1
                                
                            num_split-=1 # because they were originated from 1 parent
                
                    
                    #~ print nums_split.size, num_vir_split_prec.size, nums_split.size - num_vir_split_prec.size, num_split - num_ext
                    
                    if  nums_split.size - num_vir_split_prec.size != num_split - num_ext:
                        print "ERROR IN CLUSTERS SPLIT AND EXTINCTIONS"
                        sys.exit()
                        
                
                
                num_ext_tot[i]=num_ext
                num_split_tot[i]=num_split
                
                if i >0:
                    num_ext_per_cl_tot[i]=num_ext/float(    centroids_split_prec.shape[0])
                    num_split_per_cl_tot[i]=num_split/float(centroids_split_prec.shape[0])
                else:
                    num_ext_per_cl_tot[i]=0
                    num_split_per_cl_tot[i]=0
                    
                
                inter_centroids_dists_split=np.array([0])
                if centroids_split.shape[0]>1:
                    inter_centroids_dists_split=pdist(centroids_split)
                    
                    
                
                
                
                
                
                #~ print labels_split.shape, centroids_split.shape, nums_split.shape, sizes_split.shape
                
                vel_clust_split=np.asarray(vel_clust_split)
                
                var_parall_clust_split=np.asarray(var_parall_clust_split)
                var_perp_clust_split  =np.asarray(var_perp_clust_split  )
                var_tot_clust_split   =np.asarray(var_tot_clust_split   )
                
                vel_clust_split_alltimes.append(vel_clust_split.mean())
                
                var_parall_clust_split_alltimes.append(var_parall_clust_split.mean())
                var_perp_clust_split_alltimes.append(var_perp_clust_split.mean())
                var_tot_clust_split_alltimes.append(var_tot_clust_split.mean())
                
                
                size_clusters_split.append(sizes_split.mean())
                inter_clusters_dist_split.append(inter_centroids_dists_split.mean())
                inter_clusters_dist_max_split.append(np.amax(inter_centroids_dists_split))
                
                num_coords_clusters_split.append(nums_coords_split.mean())
                num_vir_clusters_split.append(nums_split.mean())
                n_clusters_split.append(centroids_split.shape[0])
                
                #density_split=np.where(sizes_split>0, nums_split/sizes_split, 0.)
                densities_split=nums_split[(sizes_split>0.1) & (nums_split>10)]/sizes_split[(sizes_split>0.1) & (nums_split>10)]
                density_clusters_split.append(densities_split.mean())
                
                
                centroids_split_prec          =    centroids_split.copy()
                num_vir_split_prec            =    nums_split.copy()
                num_coords_split_prec            =    nums_coords_split.copy()
                sizes_split_prec              =    sizes_split.copy()
                IDs_split_prec                =    IDs_split.copy()
                IDs_fath_split_prec           =    IDs_fath_split.copy()
                num_coords_fath_split_prec       =    num_coords_fath_split.copy()
                num_vir_fath_split_prec       =    num_vir_fath_split.copy()
                sizes_fath_split_prec         =    sizes_fath_split.copy()
                centroids_fath_split_prec     =    centroids_fath_split.copy()
                t_birth_split_prec            =    t_birth_split.copy()
                
                labels_split_prec            =    labels_split.copy()
                clust_labels_split_prec            =    clust_labels_split.copy()
                
                
            
            
            
            
            
            
            
            data_track_tot	    =np.asarray(data_track_tot             )
            print data_track_tot.shape
            
            #data_fin =data_fin.T
            #
            #print data_fin
            #print data_fin.shape
            
            outfile_track='{inp}/clusters_track_real_{real}.dat'.format(inp=dir_in, real=real)
            
            with open(outfile_track, 'w') as f_handle:
                #np.savetxt(f_handle, data_track_tot, delimiter=" ",  newline=" ",  fmt="%s")
                #f_handle.write("\n")
                np.savetxt(f_handle, data_track_tot, fmt='%15.5f')
            
            
                
            
            
            
            
            
            times_diff_fity_split_alltimes_allclust        =np.asarray(times_diff_fity_split_alltimes_allclust                     )
            diff_fity_split_alltimes_allclust              =np.asarray(diff_fity_split_alltimes_allclust                           )
            diff_fitx_split_alltimes_allclust              =np.asarray(diff_fitx_split_alltimes_allclust                           )
            ang_past_dir_shortmem_split_alltimes_allclust  =np.unwrap(np.asarray(ang_past_dir_shortmem_split_alltimes_allclust               ))
            ang_past_dir_split_alltimes_allclust           =np.unwrap(np.asarray(ang_past_dir_split_alltimes_allclust                        ))
            
            
            times_sharp_turn_split_alltimes_allclust =np.asarray(times_sharp_turn_split_alltimes_allclust              )
            times_grad_split_alltimes_allclust       =np.asarray(times_grad_split_alltimes_allclust                    )
            grad_fittest_split_alltimes_allclust     =np.asarray(grad_fittest_split_alltimes_allclust              )
            shortmem_gradx_bulk_split_alltimes_allclust       =np.asarray(shortmem_gradx_bulk_split_alltimes_allclust                )
            shortmem_grady_fittest_split_alltimes_allclust    =np.asarray(shortmem_grady_fittest_split_alltimes_allclust             )
            shortmem_grady_bulk_split_alltimes_allclust       =np.asarray(shortmem_grady_bulk_split_alltimes_allclust                )
            
            gradx_bulk_split_alltimes_allclust       =np.asarray(gradx_bulk_split_alltimes_allclust                )
            grady_fittest_split_alltimes_allclust    =np.asarray(grady_fittest_split_alltimes_allclust             )
            grady_bulk_split_alltimes_allclust       =np.asarray(grady_bulk_split_alltimes_allclust                )
            
            vel_clust_split_alltimes_allclust             =np.asarray(vel_clust_split_alltimes_allclust             )
            
            var_parall_clust_split_alltimes_allclust             =np.asarray(var_parall_clust_split_alltimes_allclust             )
            var_perp_clust_split_alltimes_allclust               =np.asarray(var_perp_clust_split_alltimes_allclust             )
            var_tot_clust_split_alltimes_allclust                =np.asarray(var_tot_clust_split_alltimes_allclust             )
            
            
            vel_clust_split_alltimes             =np.asarray(vel_clust_split_alltimes             )
            
            var_parall_clust_split_alltimes             =np.asarray(var_parall_clust_split_alltimes             )
            var_perp_clust_split_alltimes               =np.asarray(var_perp_clust_split_alltimes             )
            var_tot_clust_split_alltimes                =np.asarray(var_tot_clust_split_alltimes             )
            
            size_clusters_split             =np.asarray(size_clusters_split             )
            inter_clusters_dist_split       =np.asarray(inter_clusters_dist_split       )
            inter_clusters_dist_max_split       =np.asarray(inter_clusters_dist_max_split       )
            density_clusters_split          =np.asarray(density_clusters_split          )
            n_clusters_split                =np.asarray(n_clusters_split                )
            num_vir_clusters_split          =np.asarray(num_vir_clusters_split          )
            num_coords_clusters_split          =np.asarray(num_coords_clusters_split          )
            
            num_vir_o_nclust_tot_list = np.asarray([np.sum(i) for i in num_vir_tot_list])/n_clusters_split.astype(float)
            
            frac_1cl=(n_clusters_split[:]==1).sum()/float(n_clusters_split[:].shape[0])
            
                  
               
            
            
            
            data = np.array([times, n_clusters_dbscan_list, n_clusters_dbscan_est_list, n_clusters_dbscan_opt_kvar_list, eps_dbscan_est_list, eps_dbscan_opt_kvar_list,  size_clusters_dbscan_list, size_clusters_dbscan_est_list, size_clusters_dbscan_opt_kvar_list, density_clusters_dbscan_list, density_clusters_dbscan_est_list, density_clusters_dbscan_opt_kvar_list, n_clusters_dbscan_opt_dispersion_list, n_clusters_dbscan_opt_dispersion_CV_list, n_clusters_dbscan_opt_kCV_list, mean_dispersion_dbscan_list, mean_dispersion_dbscan_est_list, mean_dispersion_dbscan_opt_list, mean_dispersion_CV_dbscan_list, mean_dispersion_CV_dbscan_est_list, mean_dispersion_CV_dbscan_opt_list, mean_kCV_dbscan_list, mean_kCV_dbscan_est_list, mean_kCV_dbscan_opt_list, mean_kvar_dbscan_list, mean_kvar_dbscan_est_list, mean_kvar_dbscan_opt_list, inter_clusters_dist_dbscan_opt_kvar_list, inter_clusters_dist_max_dbscan_opt_kvar_list, size_clusters_split, inter_clusters_dist_split, inter_clusters_dist_max_split, density_clusters_split, n_clusters_split, num_vir_clusters_split, num_coords_clusters_split, num_ext_tot, num_split_tot, num_ext_per_cl_tot, num_split_per_cl_tot, num_vir_o_nclust_tot_list, vel_clust_split_alltimes, var_parall_clust_split_alltimes, var_perp_clust_split_alltimes, var_tot_clust_split_alltimes])
            
            data_allclust = np.array([vel_clust_split_alltimes_allclust, var_parall_clust_split_alltimes_allclust, var_perp_clust_split_alltimes_allclust, var_tot_clust_split_alltimes_allclust])
            
            data_allclust_grads = np.array([times_grad_split_alltimes_allclust, grad_fittest_split_alltimes_allclust, gradx_bulk_split_alltimes_allclust, grady_fittest_split_alltimes_allclust, grady_bulk_split_alltimes_allclust, shortmem_gradx_bulk_split_alltimes_allclust, shortmem_grady_fittest_split_alltimes_allclust, shortmem_grady_bulk_split_alltimes_allclust])
            
            
            
            data_allclust_diffy = np.array([times_diff_fity_split_alltimes_allclust, diff_fity_split_alltimes_allclust, diff_fitx_split_alltimes_allclust, ang_past_dir_shortmem_split_alltimes_allclust, ang_past_dir_split_alltimes_allclust])
            
             
            
            
            
            
            
            
            
            data = data.T
                
                
                
            file_out='{inp}/space_clustering_fct_time_real_{real}.dat'.format(inp=dir_in, real=real)
            
            np.savez_compressed(file_out,data = data, data_allclust = data_allclust, data_allclust_grads = data_allclust_grads, frac_1cl=frac_1cl, times_sharp_turn_split_alltimes_allclust= times_sharp_turn_split_alltimes_allclust)
            
            
            
        
        
            # OUTDATED
            end2end_coarse_gr_x_fct_time_list           = []
            end2end_coarse_gr_y_fct_time_list           = []
            end2end_coarse_gr_ang_fct_time_list         = []
            end2end_coarse_gr_ang_short_fct_time_list   = []
            time_diff_traj_list                         = []
            time_diff_traj_list_ang                     = []
            ttot=0
                
                
            diff_maxfit_fit_rw_fct_time    = 0
            err_diff_maxfit_fit_rw_fct_time   = 0
            chisq_maxfit_fit_rw_fct_time   = 0
        
            diff_ang_fit_rw_fct_time    = 0
            err_diff_ang_fit_rw_fct_time   = 0
            chisq_ang_fit_rw_fct_time   = 0
        
            diff_ang_short_fit_rw_fct_time    = 0
            err_diff_ang_short_fit_rw_fct_time   = 0
            chisq_ang_short_fit_rw_fct_time   = 0
            
            
            
    
                
               
                          
            dir_out_data=dir_in
            
                
                
            avg_n_clusters_dbscan_opt_kvar_list           = n_clusters_dbscan_opt_kvar_list           .mean()  
            avg_size_clusters_split                       = size_clusters_split                       .mean()  
            avg_n_clusters_split                          = n_clusters_split                          .mean()  
            avg_num_vir_clusters_split                    = num_vir_clusters_split                    .mean()  
            avg_num_coords_clusters_split                 = num_coords_clusters_split                 .mean()  
            avg_num_ext_tot                               = num_ext_tot                               .mean()  
            avg_num_split_tot                             = num_split_tot                             .mean()  
            avg_num_ext_per_cl_tot                        = num_ext_per_cl_tot                        .mean()  
            avg_num_split_per_cl_tot                      = num_split_per_cl_tot                      .mean()  
            avg_num_vir_o_nclust_tot_list                 = num_vir_o_nclust_tot_list                 .mean()  
            avg_vel_clust_split_alltimes_allclust         = vel_clust_split_alltimes_allclust         .mean()           
            avg_var_parall_clust_split_alltimes_allclust  = var_parall_clust_split_alltimes_allclust  .mean()           
            avg_var_perp_clust_split_alltimes_allclust    = var_perp_clust_split_alltimes_allclust    .mean()           
            avg_var_tot_clust_split_alltimes_allclust     = var_tot_clust_split_alltimes_allclust     .mean()           
            avg_inter_clusters_dist_split                 = inter_clusters_dist_split                 .mean()
            max_inter_clusters_dist_split                 = inter_clusters_dist_split    .max()   
            max_inter_clusters_dist_max_split             = inter_clusters_dist_max_split.max()
            var_vel_clust_split_alltimes_allclust         = vel_clust_split_alltimes_allclust         .var()   
            var_var_parall_clust_split_alltimes_allclust  = var_parall_clust_split_alltimes_allclust  .var()   
            var_var_perp_clust_split_alltimes_allclust    = var_perp_clust_split_alltimes_allclust    .var()   
            var_var_tot_clust_split_alltimes_allclust     = var_tot_clust_split_alltimes_allclust     .var()   
            rate_num_ext_tot                              = np.cumsum(num_ext_tot)         [-1]/float( times[-1])    
            rate_num_split_tot                            = np.cumsum(num_split_tot)       [-1]/float( times[-1])      
            rate_num_ext_per_cl_tot                       = np.cumsum(num_ext_per_cl_tot)  [-1]/float( times[-1])        
            rate_num_split_per_cl_tot                     = np.cumsum(num_split_per_cl_tot)[-1]/float( times[-1]) 
            
            num_sharp_turns  = times_sharp_turn_split_alltimes_allclust.size/float( times[-1]) 
            rate_sharp_turns = times_sharp_turn_split_alltimes_allclust.size/float( times[-1]) 
            
            avg_grad_fittest_split_alltimes_allclust =  grad_fittest_split_alltimes_allclust.mean()
            var_grad_fittest_split_alltimes_allclust =  grad_fittest_split_alltimes_allclust.var()
            
            
            
            # OUTDATED STATS
            avg_shortmem_err_ys   = 0
            avg_shortmem_err_svsx = 0
            
            frac_shortmem_err_ys  =0
            frac_shortmem_err_svsx=0
            
            

            data = np.array([frac_1cl, avg_n_clusters_dbscan_opt_kvar_list, avg_size_clusters_split, avg_n_clusters_split, avg_num_vir_clusters_split, avg_num_coords_clusters_split, avg_num_ext_tot, avg_num_split_tot, avg_num_ext_per_cl_tot, avg_num_split_per_cl_tot, avg_num_vir_o_nclust_tot_list, avg_vel_clust_split_alltimes_allclust, avg_var_parall_clust_split_alltimes_allclust, avg_var_perp_clust_split_alltimes_allclust, avg_var_tot_clust_split_alltimes_allclust, avg_inter_clusters_dist_split, max_inter_clusters_dist_split, max_inter_clusters_dist_max_split, var_vel_clust_split_alltimes_allclust, var_var_parall_clust_split_alltimes_allclust, var_var_perp_clust_split_alltimes_allclust, var_var_tot_clust_split_alltimes_allclust, rate_num_ext_tot, rate_num_split_tot, rate_num_ext_per_cl_tot, rate_num_split_per_cl_tot, num_sharp_turns, rate_sharp_turns, avg_shortmem_err_ys, avg_shortmem_err_svsx, frac_shortmem_err_ys, frac_shortmem_err_svsx,    diff_maxfit_fit_rw_fct_time, err_diff_maxfit_fit_rw_fct_time, chisq_maxfit_fit_rw_fct_time, diff_ang_fit_rw_fct_time, err_diff_ang_fit_rw_fct_time, chisq_ang_fit_rw_fct_time, diff_ang_short_fit_rw_fct_time, err_diff_ang_short_fit_rw_fct_time, chisq_ang_short_fit_rw_fct_time ])
            
            
            file_out='{inp}/global_features_virclust_{real}.txt'.format(inp=dir_out_data, real=real)
        
            with open(file_out,'a') as f_handle:
                np.savetxt(f_handle, data, fmt='%15.15f', newline=" ")
                f_handle.write("\n")
        
            
            
            
            
            
            
            #~ PLOTS FIG 3B
            
            fig = plt.figure(figsize=(thisfigsize[0]/2.,thisfigsize[0]/2.))
            grid = gridspec.GridSpec(1, 1, left=0.12, right=0.97, top=0.93, bottom=0.13,
                 wspace=0.4, hspace=0.35)
            labeled_axes = []
            ax = plt.Subplot(fig, grid[0, 0])
            fig.add_subplot(ax)
            labeled_axes.append(ax)
            
            for i_centr, x_centr_serie in enumerate(xs_all_centr_all_time):
                x_centr_serie=np.asarray(x_centr_serie)
                y_centr_serie=np.asarray(ys_all_centr_all_time[i_centr])
                    
                
                
                #~ print i_centr, x_centr_serie.shape, y_centr_serie.shape
                #~ print x_centr_serie
                #~ print y_centr_serie
                
                ax.plot(x_centr_serie, y_centr_serie, linestyle='-', color='g')
                #~ ax.plot(x_centr_serie, y_centr_serie, linestyle='-', color='grey', linewidth=0.5)
            
            
            #    ax.plot(IS_x, IS_y, linestyle='-', color='r', label='average IS')
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            
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
            
            if maxdiff>2*recog_width:
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
              loc='upper right', bbox_to_anchor=(2.5, 1.2))
            mpsetup.despine(ax) 
            
            
            #### finish figure ####
            labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
                 xycoords=('axes fraction'), fontweight = 'bold')
            #    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
            #mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
            out_file='{out}/vir_alltimes_cartoon{real}.pdf'.format(out=dir_out_plots, real=real)
            #    print out_file, dpi=300
            fig.savefig(out_file)
            fig.clf()
            plt.close('all')
                
 
 
            ## plot  number clusters as fct of time
            
            fig = plt.figure(figsize=thisfigsize)
            grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.22,
                         wspace=0.4, hspace=0.35)
            labeled_axes = []
            ax = plt.Subplot(fig, grid[0, 0])
            fig.add_subplot(ax)
            labeled_axes.append(ax)
            ax.plot(times, n_clusters_split, linestyle='-', color='r', label='dbscan split')
            ax.set_xlabel('time (y)')
            ax.set_ylabel('number of clusters')
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
            out_file='{out}/space_clusters_1curve_{real}.png'.format(out=dir_out_plots, real=real)
            #    print out_file
            fig.savefig(out_file)
            
            
            
            
            
            
            ## plot  number clusters as fct of time
            
            fig = plt.figure(figsize=thisfigsize)
            grid = gridspec.GridSpec(1, 2, left=0.15, right=0.97, top=0.91, bottom=0.2,
                         wspace=0.4, hspace=0.35)
            labeled_axes = []
            ax = plt.Subplot(fig, grid[0, 0])
            fig.add_subplot(ax)
            labeled_axes.append(ax)
            ax.plot(times, n_clusters_dbscan_est_list, linestyle='-', color='b', label='dbscan eps est')
            ax.plot(times, n_clusters_dbscan_opt_dispersion_CV_list, linestyle='-', color='m', label='dbscan opt dispersion CV')
            ax.plot(times, n_clusters_dbscan_opt_dispersion_list, linestyle='-', color='g', label='dbscan opt dispersion')
            ax.plot(times, n_clusters_dbscan_opt_kCV_list, linestyle='-', color='k', label='dbscan opt kth \n dist CV')
            #ax.plot(times, n_clusters_dbscan_list, linestyle='-', color='r', label='dbscan')
            ax.plot(times, n_clusters_dbscan_opt_kvar_list, linestyle='-', color='y', label='dbscan opt kth \n dist var')
            ax.plot(times, n_clusters_split, linestyle='-', color='r', label='dbscan split')
            ax.set_xlabel('time (y)')
            ax.set_ylabel('number of clusters')
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
            out_file='{out}/space_clusters_{real}.png'.format(out=dir_out_plots, real=real)
            #    print out_file
            fig.savefig(out_file)
            
            
            
            
            fig.clf()
            plt.close('all')
            
            
            
            
            
            gc.collect()
            
            
            
            
            
