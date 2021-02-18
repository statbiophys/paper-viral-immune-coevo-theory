#!/bin/bash

# runs the model wanted. Wants the "input_file" and the number of realizations I want. The input files are kept in the same dir as the c programs: "../simulazioni/c_programs_modelli/nome_modello/input.dat
#puts the output in the required directory, and prints warnings in a specific file. Moreover copies the input file as backup.


# cycle on ppl_num, f_m, mu, sigma
if (( $# != 2 ))
then

    echo "Usage $0, <input_file> <n_r> "
	
else

    input=$1 # cerr<<"Usage "<<argv[0]<<" <mu> <recognition width> <sigma> <eps> <I> <people number> <coverage save rate> <out_dir>  <fraction infected> <number infections> <full_frame save time> <n_reals> <real>" <<endl;
    n_r=$2 # numero realizzazioni
    #n_full=$3
    
    #estrae i parametri dal file, li salvo in un vettore per comodità, assumo che sia rispettato l'ordine specificato nel main.cpp, n_rwal e real non ci sono
    j=0
    for i in `awk 'NR==2 {print $0;}' $input`
    do
	    
	     p[j]=$i
	     
	     ((j++))
	    
    done
    

    mu=${p[0]}
    rec_width=${p[1]}
    jump_size=${p[2]}	
    ppl_num=${p[3]}	
    data_dir=${p[4]}	
    maxinf=${p[5]}	
    save_frames=${p[6]}	
    initial_condition=${p[7]}
    save_final_configuration=${p[8]}
    kernel=${p[9]}
    save_phylo=${p[10]}
    fake_in_cond=${p[11]}
    mem_points=${p[12]}
    F0=${p[13]}

    if [ ! -d $data_dir ]; then
        mkdir $data_dir
    fi    
    

    data_dir=${data_dir}/mem_${mem_points}/
    
		
		
	 
    py_script_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/python_scripts/coarse_grained_plots # scripts directory
    script_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/scripts_coarse_grained # scripts directory
    cprogram_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/cprograms/coarse_grained_coevo # scripts directory
    
    #mkdir /home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/logs_scan_manyvir
    mkdir /home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/logs_coarse_gr_mem_1_F0_1
    mkdir /home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/logs_plot_coarse_gr_mem_1_F0_1
     
    #data_dir=/home/jaco/Documents/immunitary_systems_viruses_coevo/viral_immune_coevo/contagion_simulation_results/prova_folder/
    curr_dir=`pwd`
    
    declare -A par_id=() # associative array, requires bash 4
    count_tot=0
    rm -rf /users/marchi/.cache/matplotlib/.matplotlib_lock-*
    
    
    
 ##### CYCLES ON PARAMETERS TO SWEEP
 
    for mem_points in  1 
    do



    for F0 in 1. #  
    do



    for ppl_num in  1000000000000 # 100000000 # 1000000000000 # 1000000 # 10000000 
    #~ for ppl_num in  10000000000 # 100000000 # 1000000000000 # 1000000 # 10000000 
    #~ for ppl_num in  100000000 # 100000000 # 1000000000000 # 1000000 # 10000000 
    do
     
	for mu in   0.000001  0.00001  0.0001  0.001   0.01  # ppl 1000000000000
	#~ for mu in   0.000001  0.00001  0.0001  0.001   0.01   0.1   # ppl 10000000000
	#~ for mu in   0.0001  0.001   0.01   0.1   # ppl 100000000
	do
	
	    unset rec_width_mu
	    rec_width_mu='()'
	    
	    	
	    if [[ $mu == 0.1 ]] ; then
	    
	    

		j=0
	    
		for rec_width in  500. 700. 1000. 1500. 2000. 3000. 5000.  # ppl 10000000000
		#~ for rec_width in  300. 500. 700. 1000. 1500. 2000. 3000. # ppl 100000000
		do
		     rec_width_mu[j]=$rec_width
		     
		     echo $j, ${rec_width_mu[j]}, $mu
		     
		     ((j++))		
		done
	    
	    	
	    elif [[ $mu == 0.01 ]] ; then
	    
	    

		j=0
	    
		for rec_width in  500. 700. 1000. 1500. 2000. 3000. 5000.  # ppl 1000000000000
		#~ for rec_width in 200. 300. 500. 700. 1000. 1500. 2000. 3000. 5000.  # ppl 10000000000
		#~ for rec_width in 200. 300. 500. 700. 1000.  # ppl 100000000
		do
		     rec_width_mu[j]=$rec_width
		     
		     echo $j, ${rec_width_mu[j]}, $mu
		     
		     ((j++))		
		done
	    
	    elif [[ $mu == 0.001 ]] ; then
	    
		#sigma_mu='()'
		j=0
	    
        for rec_width in 150. 200. 300. 500. 700. 1000. 1500. 2000. 3000. 5000.  # ppl 1000000000000
        #~ for rec_width in 200. 300. 500. 700. 1000. 1500. 2000. 3000. 5000.  # ppl 10000000000
		#~ for rec_width in 70. 80. 100. 150. 200. 300. 500. # ppl 100000000
		do
		     rec_width_mu[j]=$rec_width
		     echo $j, ${rec_width_mu[j]}, $mu
		     
		     ((j++))		
		done
		
	    elif [[ $mu == 0.0001 ]] ; then
	    
		#sigma_mu='()'
		j=0
	    
        for rec_width in 150. 200. 300. 500. 700. 1000. 1500. 2000. 3000. 5000.  # ppl 1000000000000
        #~ for rec_width in 100. 150. 200. 300. 500. 700. 1000.  # ppl 10000000000
		#~ for rec_width in  40. 50. 60. 70. 80. # ppl 100000000
		do
		     rec_width_mu[j]=$rec_width
		     echo $j, ${rec_width_mu[j]}, $mu
		     
		     ((j++))		
		done
	    
	    elif [[ $mu == 0.00001 ]]; then
	    
		#sigma_mu='()'
		j=0
	    
        for rec_width in 100. 150. 200. 300. 500. 700. 1000. 1500. 2000. 3000.  # ppl 1000000000000
		#~ for rec_width in 100. 150. 200. 300. 500. 700. 1000.  # ppl 10000000000
		do
		     rec_width_mu[j]=$rec_width
		     
		     echo $j, ${rec_width_mu[j]}, $mu
		     ((j++))		
		done
	    
	    elif [[ $mu == 0.000001 ]]; then
	    
		#sigma_mu='()'
		j=0
	    
        for rec_width in 100. 150. 200. 300. 500. 700. 1000. 1500. 2000. 3000.  # ppl 1000000000000
		#~ for rec_width in 30. 50. 70. 100. 150. 200. 300. 500. 700. # ppl 10000000000
		do
		     rec_width_mu[j]=$rec_width
		     
		     echo $j, ${rec_width_mu[j]}, $mu
		     ((j++))		
		done
	    fi
	    
	    echo "${rec_width_mu[@]}"
	    
	    
	    
	    
		#for sigma in  0.01 0.03 0.05 0.07 0.1 0.3 0.5 0.7 1. 10. 100.
		#for sigma in 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.01 0.03 0.05  0.1  0.5  100.
		for rec_width in "${rec_width_mu[@]}"
		#for sigma in 0.002  0.004  0.006 
		do
        
        ### HANDLE FOLDERS STRUCTURE

            param_dir_tmp="D_2_pplnum_${ppl_num}_maxt_${maxinf}_mu_${mu}_rec_width_${rec_width}_jump_size_${jump_size}_F0_${F0}" #
            
            param_dir_d=`echo $param_dir_tmp | sed 's/\./d/g' | sed 's/\s\s*/_/g'` # substitutes dots with d, and spaces and tabs with _
            
            param_dir=`echo $param_dir_d | sed 's/d0\+_/d_/g' | sed 's/d\(0*[123456789][123456789]*\)0\+_/d\1_/g'`
            
            data_dir_fin=${data_dir}/${param_dir}/
            if [ ! -d $data_dir ]; then
                mkdir $data_dir
             fi    
            
            cd $py_script_dir
             if [ ! -d $data_dir_fin ]; then
                mkdir $data_dir_fin
             fi
        
            if [ ! -d ${data_dir_fin}/realizations ]; then
                mkdir ${data_dir_fin}/realizations
             fi


		    		    
		    file_rerun_log=${data_dir_fin}log_reruns.txt
		    
		    #~ rerun_lines=`cat $file_rerun_log | wc -l`
		    #~ rerun_last5=`cat $file_rerun_log  | grep "." | tail -5 | head -1`


    
			    rm ${cprogram_dir}/submit_file_${param_dir}.dat
			    
			    for (( r=1; r<=$n_r; r++ )) # runs the program to create the required number of realizations
			    do
			    # cerr<<"Usage "<<argv[0]<<" <mu> <recognition width> <sigma> <eps> <I> <people number> <coverage save rate> <out_dir>  <fraction infected> <number infections> <full_frame save time> <n_reals> <real>" <<endl;
			    
			    
				echo "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0" >> ${cprogram_dir}/submit_file_${param_dir}.dat
				
				((count++))
				
				#echo "run "$count
				    
			    done
                
                
                    # RUN THE MODEL


                    par_id[$count_tot]=$(sbatch --exclude=zuzia --array=1-$n_r  -J par_1_1_${count_tot} ${script_dir}/arrayjob_coarse_grained_slurm.sh "$data_dir_fin" "$r" "$count" "$n_r" "$param_dir")
                    
                    echo "from bash ", ${par_id[$count_tot]}
                    
                    contagion_curr=${par_id[$count_tot]}
                    
                    echo "parscan id ", ${contagion_curr##* }
                    


	
		    
			    # PLOT AND ANALYZE
			    #
			    plot_curr=$(sbatch --exclude=zuzia -J plot_1_1_${count_tot}  --dependency=afterany:${contagion_curr##* } ${script_dir}/plot_all_slurm.sh "$data_dir_fin" "$param_dir")
			    #plot_curr=$(sbatch -J plot_${count_tot}  ${script_dir}/plot_all_slurm.sh "$data_dir_fin" "$param_dir")
			    
			    
			    # and now remove redundant files. Process wrapped in a script
			    
			    cd $script_dir
			    
			    #echo "REMOVE DATA"
			    #
			    #qsub -N rem_${count_tot} -hold_jid plot_${count_tot} ${script_dir}/remove_data.sh "$data_dir_fin"
			
				
			    count_tot_prec=$(($count_tot-5))
				
			    #contagion_prec=par_${count_tot}
			    
			    echo "saving ", ${par_id[$count_tot_prec]}
			    contagion_prec=${par_id[$count_tot_prec]}
			    
			    cd $curr_dir
			    echo " "
			    echo "ppl_num = "$ppl_num" mem_points = "$mem_points" f0 = "$f0" mu = "$mu" rec_width = "$rec_width", run "$count_tot
			    ((count_tot++))
		done   
	    done
	    
	done
		
    done
    done

fi
