#!/bin/bash

# load current environment variables to context of the job
#SBATCH -p Short
#SBATCH --export=[ALL] 
# combine error and normal output into a single file 
# output in specified dir 
#SBATCH -o /home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/logs_plot_coarse_gr_mem_1_F0_1/slurm-%j.out
#SBATCH -n 1 # number of cores
#SBATCH --tasks-per-node=2 # cpu allocated for that job on the node(s) it's sent to. This means that if there are 28 cores, if I put 2 I'll be sending 14 jobs per node (if not parallelized and - n 1)
#SBATCH --mem-per-cpu=20000
#SBATCH -t 18-2:00 # time (D-HH:MM)

##### RUNS A BUNCH OS ANALYSIS AND PLOTTING PYTHON SCRIPTS ON THE OUTPUT OF THE MODEL SIMULATION


:'
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
'

 
py_script_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/python_scripts/coarse_grained_plots # scripts directory
script_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/scripts_coarse_grained # scripts directory
#cprogram_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/cprograms/coarse_grained_coevo # scripts directory
   
curr_dir=`pwd`



if (( $# != 2 ))
then

    echo "Usage $0, <input_dir> <param_dir> "
	
else

    echo  $SLURM_SUBMIT_HOST
    
    echo

    data_dir_fin=$1 # directory with input files
    param_dir=$2 # directory with input files
    #file_plot_log=${dir_io}log_plots.txt


    ulimit -S -v 4500000 # maximum RAM after which process is killed , in kb


    cd $py_script_dir
    

    
		    echo "PLOTTING AVERAGE DYNAMICS"
		    echo
		    
		    python plot_avg_dyn.py "$data_dir_fin" #>> "$file_plot_log"
		    echo 
            
            framedirs=${data_dir_fin}/realizations/realization_1/frames
            #~ framedirs=${data_dir_fin}/realizations/realization_${r}/frames


            echo
            cd ${data_dir_fin}/realizations/realization_1/
            pwd
            echo "unzipping"
            if [ ! -n "$(ls -A $framedirs 2>/dev/null)" ] && [ -f "${framedirs}.tar.gz" ]
            then
              echo "empty (or does not exist), untar" 
                cp "${framedirs}.tar.gz" "${framedirs}_bkp.tar.gz"
              rm -rf "$framedirs"
              tar -zxf "${framedirs}.tar.gz"
            else
              echo "contains files (or is a file)"
            fi
               
            cd $py_script_dir
            
     
     
		    echo "VIRAL CLUST" #>> "$file_plot_log"
		    echo #>> "$file_plot_log"
		    
		    python viral_clustering.py "$data_dir_fin" #>> "$file_plot_log"
		    echo 
            
            echo
            echo "zipping"
            cd "${framedirs}/.."
            framedirs_here=frames
            pwd
           if [ ! -n "$(ls -A $framedirs_here 2>/dev/null)" ] 
            then
              echo "empty (or does not exist)" 
            else
              echo "contains files (or is a file), zip"
                rm -rf "${framedirs_here}.tar.gz"
                tar -zcf "${framedirs_here}.tar.gz" "$framedirs_here"
                rm -rf "$framedirs_here"
              
              
            fi
            echo
            
            cd $py_script_dir
            
        
		    
		    echo "PERS LEN" #>> "$file_plot_log"
		    echo #>> "$file_plot_log"
		    
		    python persistence_length.py "$data_dir_fin" #>> "$file_plot_log"
		    echo 
		    
            
		    echo "RSPEC" #>> "$file_plot_log"
		    echo #>> "$file_plot_log"
		    
		    python rspeciations.py "$data_dir_fin" #>> "$file_plot_log"
		    echo 
		    
         
		    cd $curr_dir
		    echo " "

    ulimit -S -v unlimited # maximum ram after which process is killed , in kb

    cout "END"
fi
