#!/bin/bash
# load current environment variables to context of the job
#SBATCH -p Infinite
#SBATCH --export=[ALL] 
# combine error and normal output into a single file 
# output in specified dir 
#SBATCH -o /home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/logs_coarse_gr_mem_1_F0_1/slurm-%A_%a.out
#SBATCH -n 1 # number of cores
#SBATCH --tasks-per-node=1 # cpu allocated for that job on the node(s) it's sent to. This means that if there are 28 cores, if I put 2 I'll be sending 14 jobs per node (if not parallelized and - n 1)
#SBATCH --mem-per-cpu=20000
#SBATCH -t 140-2:00 # time (D-HH:MM)


##### RUNS THE MODEL SIMULATION AND  HANDLES REINITIALIZATIONS UPON EXTINCTIONS AND EXPLOSIONS


if (( $# != 5 ))
then

    echo "Usage $0, <input_file> <r> <count> <n_r> <param_dir> "
	
else

    data_dir_fin=$1 # directory with input files
    r=$2 # directory with input files
    count=$3 # directory with input files
    n_r=$4 # directory with input files
    param_dir=$5 # directory with input files
    
    ulimit -S -v 2500000


    py_script_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/python_scripts/coarse_grained_plots # scripts directory
    script_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/scripts_coarse_grained # scripts directory
    cprogram_dir=/home/marchi/immunitary_systems_viruses_coevo/viral_immune_coevo/cprograms/coarse_grained_coevo # scripts directory
    
    INFILE=${cprogram_dir}/submit_file_${param_dir}.dat
    #~ params=`awk "NR==$SLURM_ARRAY_TASK_ID" $INFILE`
    params=`awk "NR==1" $INFILE`
    
	echo "$INFILE" "$SLURM_ARRAY_TASK_ID"  $HOSTNAME
	echo "$params" 

    j=0
    for i in $params
    do
	    
	     p[j]=$i
	     
	     echo ${p[j]}
	     
	     ((j++))
	    
    done
    

    mu=${p[0]}
    rec_width=${p[1]}
    jump_size=${p[2]}	
    ppl_num=${p[3]}	
    data_dir=${p[4]}	
    maxinf=${p[5]}	
    save_frames=${p[6]}	
    n_r=${p[7]}
    r=${p[8]}
    initial_condition=${p[9]}
    save_final_configuration=${p[10]}
    kernel=${p[11]}
    save_phylo=${p[12]}
    fake_in_cond=${p[13]}
    mem_points=${p[14]}
    F0=${p[15]}

    
# "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0"




    IC_dir_fin=${data_dir_fin}/realizations/realization_${r}/"fullstatus_last"


    
    
    curr_dir=`pwd`
    
    count=0
    
    cd $py_script_dir
 
 
	
	###########   FIRST RUN        #######################
	
    
	
	#initial_condition="true" # this doesn't start from IC
	initial_condition="false" # this doesn't start from IC
      
  
	echo "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points"  "$F0" 
	
	
	${cprogram_dir}/coarse_grained_coevo  "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0"  # >> ${cprogram_dir}/warn_tmp.txt
	 #DEBUGGING
    #gdb  -ex=r --args  ${cprogram_dir}/coarse_grained_coevo  "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0"  #>> ${cprogram_dir}/warn_tmp.txt
    #gdb  --command=${cprogram_dir}/gdbscript.gdb --args  ${cprogram_dir}/coarse_grained_coevo  "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0"  #>> ${cprogram_dir}/warn_tmp.txt
	
	#~ mv ${cprogram_dir}/warn_tmp.txt ${data_dir_fin}/realizations/realization_${r}/warnings.txt # 1 0.1 0.1 0.1 100000 20000000 
	
	#~ rm -f ${data_dir_fin}/realizations/realization_${r}/initial_condition.dat
	
	#~ if [ "$save_phylo" = "false" ] ; then
	    #~ echo 'remove phylogeny file!'
	    #~ rm -f ${data_dir_fin}/realizations/realization_${r}/viral_phylogeny*
	#~ fi

	
      
	((count++))
	
	echo "run "$count
	    
    #done
    cd $py_script_dir





    file_rerun_log=${data_dir_fin}/log_reruns_${r}.txt

    echo "CHECK FIRST RUN" >> "$file_rerun_log"
    
    
   ####  CHECK FOR EXTINCTIONS AND EXPLOSIONS IN FIRST RUN
    
    #which python
    #python -V
    #
    #
    #
    #python steady_state_check_nophylo.py "$data_dir_fin" >> "$file_rerun_log"
    #
    #echo >> "$file_rerun_log"
    #
    #rerun=`grep "." "$file_rerun_log" | tail -1`    
    #extinct=`grep "." "$file_rerun_log" | tail -2 | head -1`    
    
    
    tfin=0
        
    if [ -s ${data_dir_fin}/realizations/realization_${r}/"evo_mean_stats_real_${r}.dat" ]; then
       
       evolines=`cat ${data_dir_fin}/realizations/realization_${r}/"evo_mean_stats_real_${r}.dat" | sed '/^\s*#/d;/^\s*$/d' | wc -l`
       
       if (( $evolines > 1 )); then
           tfin=`cat ${data_dir_fin}/realizations/realization_${r}/"evo_mean_stats_real_${r}.dat" |  tail -1 | awk '{printf "%i", $1;}' `
       fi
    fi
    
    rerun="False"
    extinct="False"
    explode="False"
    
    misstime=`echo "$maxinf"  "$tfin" | awk '{ printf "%i", $1 - $2; }'`
    
    if (( $misstime > 5000 )); then
    
        if [ -s ${data_dir_fin}/realizations/realization_${r}/"extinct_file.txt" ]; then
            rerun="True"
            extinct="True"
        fi

        if [ -s ${data_dir_fin}/realizations/realization_${r}/"expl_file.txt" ]; then
            rerun="True"
            explode="True"
        fi
    fi
    
    rm -f ${data_dir_fin}/realizations/realization_${r}/"extinct_file.txt"
    rm -f ${data_dir_fin}/realizations/realization_${r}/"expl_file.txt"
    
    
    echo "from bash: rerun "$rerun", extinct "$extinct", explode "$explode", maxinf "$maxinf", tfin "$tfin >> "$file_rerun_log"
    echo "from bash: rerun "$rerun", extinct "$extinct", explode "$explode", maxinf "$maxinf", tfin "$tfin
    
    #mv ${data_dir_fin}/realizations/realization_${r}/*.pdf ${IC_dir_fin}/.
    
    echo  >> "$file_rerun_log"
    echo
    
    file_rec_events=${data_dir_fin}/"cumulative_events_rerun_${r}.txt"

    #echo "# 1 time_check     2 runcount      3 time_tot_elaps     4 extinctcount    5 explodecount    6 probextinct    7 probexplode" > $file_rec_events
    
    tprec=0
    
    runcount=0
    tottime_erased=0
    extinctcount=0
    probextinct=0
    explodecount=0
    probexplode=0
    extinctcount_same=0
    explodecount_same=0
    
    if [ -s $file_rec_events ]; then
       
       lines=`cat $file_rec_events  | sed '/^\s*#/d;/^\s*$/d' | wc -l`
       
       if (( $lines > 0 )); then
            tprec=`cat $file_rec_events |  tail -1 | awk '{printf "%i", $1;}' `

            runcount=`cat $file_rec_events |  tail -1 | awk '{printf "%i", $2;}' `
            tottime_erased=`cat $file_rec_events |  tail -1 | awk '{printf "%i", $3 - $1;}' `
            extinctcount=`cat $file_rec_events |  tail -1 | awk '{printf "%i", $4;}' `
            explodecount=`cat $file_rec_events |  tail -1 | awk '{printf "%i", $5;}' `
       fi
    fi    
    
    
    #### RERUN AS MANY SIMULATION AS NECESSARY, ADVANCING THE INITIALIZATION FILE, TRACKING THE CUMULATIVE CONSECUTIVE EXPLOSIONS AND EXTINCTIONS, AND REMOVING ALL OUTPUT DATA AFTER THE REINITIALIZATION POINT 
    
    while [ "$rerun" = "True" ]; do
  	((runcount++))

	
        echo "RUN "$runcount >> "$file_rerun_log"
        echo "RUN "$runcount
		cd $py_script_dir
        
        #### CHOOSE THE INITIALIZATION FILE,
	    
		initial_condition="false" # this  start from IC
        
            
        tprec_prec=$tprec
        tprec=0
        
        if [ -d $IC_dir_fin ]; then
            
            if [ -s ${IC_dir_fin}/"miscellaneous.dat" ]; then
            
               
               lines=`cat ${IC_dir_fin}/"miscellaneous.dat"  | sed '/^\s*#/d;/^\s*$/d' | wc -l`
               
               if (( $lines > 0 )); then
                    echo "use last checkpoint  "  >> "$file_rerun_log"
               
                   tprec=`cat ${IC_dir_fin}/"miscellaneous.dat" |  tail -1 | awk '{printf "%i", $1;}' `
               fi
            fi
    
        fi 
        
        throwtime=`echo "$tfin"  "$tprec" | awk '{ printf "%i", $1 - $2; }'`
        
        if (( $throwtime < 1000 )); then
        
            echo "use previous checkpoint  "  >> "$file_rerun_log"
            echo "use previous checkpoint  "
        
            rm -rf $IC_dir_fin
            
            cp -r ${data_dir_fin}/realizations/realization_${r}/"fullstatus_sec_last" $IC_dir_fin
            
                
            tprec=0
            
            if [ -d $IC_dir_fin ]; then
                    
                if [ -s ${IC_dir_fin}/"miscellaneous.dat" ]; then
                   
                   lines=`cat ${IC_dir_fin}/"miscellaneous.dat"  | sed '/^\s*#/d;/^\s*$/d' | wc -l`
                   
                   if (( $lines > 0 )); then
                       tprec=`cat ${IC_dir_fin}/"miscellaneous.dat" |  tail -1 | awk '{printf "%i", $1;}' `
                   fi
                fi
        
            fi 
        
        fi
        
        if (( $tprec > 0 )); then
        
           initial_condition="true"
           
        fi
        
        ####  ABORT IF THE CUMULATIVE CONSECUTIVE EXPLOSIONS AND EXTINCTIONS IS 20
        

		if (( $tprec_prec == $tprec )); then
        
            if (( $extinctcount_same > 20 )); then
                echo "too many extinction, aborting script"  >> "$file_rerun_log"
                echo "too many extinction, aborting script" 
                probextinct=1
                break 2  # break [n] sintax
            fi
            if (( $explodecount_same > 20 )); then
                echo "too many expl, aborting script"  >> "$file_rerun_log"
                echo "too many expl, aborting script"
                probexplode=1
                break 2  # break [n] sintax
            fi
        else
            extinctcount_same=0
            explodecount_same=0
		fi	
        
        #### REMOVE ALL OUTPUT DATA AFTER THE REINITIALIZATION POINT 
        
        echo $IC_dir_fin  >> "$file_rerun_log"
        cat ${IC_dir_fin}/"miscellaneous.dat"  >> "$file_rerun_log"
        cat ${IC_dir_fin}/"miscellaneous.dat"  | sed '/^\s*#/d;/^\s*$/d' | wc -l   >> "$file_rerun_log"

        tottime_erased=`echo "$tfin"  "$tprec" "$tottime_erased" | awk '{ printf "%i", $3 + $1 - $2; }'`
        #tottime_erased=$(( $tottime_erased + $tfin - $tprec ))
        
        echo "checkpoint time "$tprec", previous "$tprec_prec", tfin "$tfin", initial_condition "$initial_condition", tottime_erased "$tottime_erased  >> "$file_rerun_log"
        echo "checkpoint time "$tprec", previous "$tprec_prec", tfin "$tfin", initial_condition "$initial_condition", tottime_erased "$tottime_erased
        
        throwndirs=${data_dir_fin}/realizations/realization_${r}/thrown_run_${runcount}
        
        mkdir $throwndirs
        
        for framefile in `ls ${data_dir_fin}/realizations/realization_${r}/frames/*antigenic_space_*`; do
            
            filename=$(basename -- "$framefile")
            extension="${filename##*.}"
            filename="${filename%.*}"
            
            timeframe="${filename##*_}"

            if (( $timeframe > $tprec )); then
            
               #~ rm -f $framefile
               mv $framefile ${throwndirs}/.
               
            fi            
        done
        
        
        # run clustering on thrown frames
        
        python viral_clustering_throwns.py "$throwndirs" "$r" # >> "$file_rerun_log"
            
        echo
        echo "zipping"
        cd "${throwndirs}/.."
        pwd
        rm -rf "${throwndirs}.tar.gz"
        tar -zcvf "${throwndirs}.tar.gz" "$throwndirs"
        rm -rf "$throwndirs"
        echo
        
        cd $py_script_dir
        
        
        
        for timefile in `ls ${data_dir_fin}/realizations/realization_${r}/*.dat`; do # carefull to separate properly dat and txt
            if [[ ! $timefile =~ "antigenic_space_" ]];then
             
                mv  $timefile ${timefile}_tmp
                
                if [[ $timefile =~ "evo_mean_" ]];then
                    awk -v thr=$tprec '(/#/ || $1<=thr) {print $0;}' ${timefile}_tmp > $timefile # this is the most important
                else
                    awk -v thr=$tprec 'NR<=thr {print $0;}' ${timefile}_tmp > $timefile
                fi
                
                rm -f ${timefile}_tmp
            fi
        done
        
        
        ####  TRACK THE CUMULATIVE CONSECUTIVE EXPLOSIONS AND EXTINCTIONS 
		
		if [ "$extinct" = "True" ]
		then
		    ((extinctcount++))
		    echo "precedent run went extinct. Number extinctions = "$extinctcount  >> "$file_rerun_log"
		    echo "precedent run went extinct. Number extinctions = "$extinctcount
            if (( $tprec_prec == $tprec )); then
                ((extinctcount_same++))
                 echo "extinctcount_same = "$extinctcount_same  >> "$file_rerun_log"
                 echo "extinctcount_same = "$extinctcount_same
            fi
		fi
		
		if [ "$explode" = "True" ]
		then
		    ((explodecount++))
		    echo "precedent run exploded. Number expl = "$explodecount  >> "$file_rerun_log"
		    echo "precedent run exploded. Number expl = "$explodecount
            if (( $tprec_prec == $tprec )); then
                ((explodecount_same++))
                 echo "explodecount_same = "$explodecount_same  >> "$file_rerun_log"
                 echo "explodecount_same = "$explodecount_same
            fi
		fi
        
  
        echo "$tprec"  "$runcount" "$tottime_erased" "$extinctcount" "$explodecount" "$probextinct" "$probexplode" | awk 'BEGIN{OFS="\t"} { print $1, $2, $3 + $1, $4, $5, $6, $7; }' >> $file_rec_events


        rm -f ${data_dir_fin}/realizations/realization_${r}/"extinct_file.txt"
        rm -f ${data_dir_fin}/realizations/realization_${r}/"expl_file.txt"
        
            
        #### RERUN
              
        echo "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points"  "$F0" 
        
        ${cprogram_dir}/coarse_grained_coevo  "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0" # >> ${cprogram_dir}/warn_tmp.txt
        
        #DEBUGGING
        #gdb  -ex=r --args  ${cprogram_dir}/coarse_grained_coevo  "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0"  #>> ${cprogram_dir}/warn_tmp.txt
        #gdb  --command=${cprogram_dir}/gdbscript.gdb --args  ${cprogram_dir}/coarse_grained_coevo  "$mu"  "$rec_width" "$jump_size" "$ppl_num" "$data_dir" "$maxinf" "$save_frames" "$n_r" "$r" "$initial_condition" "$save_final_configuration" "$kernel"  "$save_phylo" "$fake_in_cond" "$mem_points" "$F0"  #>> ${cprogram_dir}/warn_tmp.txt

        #~ mv ${cprogram_dir}/warn_tmp.txt ${data_dir_fin}/realizations/realization_${r}/warnings.txt # 1 0.1 0.1 0.1 100000 20000000 

		
		((count++))
		
		echo "run "$count
		    
	    #done
	    cd $py_script_dir
	
        echo "CHECK RUN "$runcount >> "$file_rerun_log"
        
         ####  CHECK FOR EXTINCTIONS AND EXPLOSIONS 
	
        
        tfin=0
            
        if [ -s ${data_dir_fin}/realizations/realization_${r}/"evo_mean_stats_real_${r}.dat" ]; then
           
           evolines=`cat ${data_dir_fin}/realizations/realization_${r}/"evo_mean_stats_real_${r}.dat" | sed '/^\s*#/d;/^\s*$/d' | wc -l`
           
           if (( $evolines > 1 )); then
               tfin=`cat ${data_dir_fin}/realizations/realization_${r}/"evo_mean_stats_real_${r}.dat" |  tail -1 | awk '{printf "%i", $1;}' `
           fi
        fi
        
        rerun="False"
        extinct="False"
        explode="False"
            
        misstime=`echo "$maxinf"  "$tfin" | awk '{ printf "%i", $1 - $2; }'`
        
        if (( $misstime > 5000 )); then
        
            if [ -s ${data_dir_fin}/realizations/realization_${r}/"extinct_file.txt" ]; then
                rerun="True"
                extinct="True"
            fi
    
            if [ -s ${data_dir_fin}/realizations/realization_${r}/"expl_file.txt" ]; then
                rerun="True"
                explode="True"
            fi
        fi
        
        echo "from bash: rerun "$rerun", extinct "$extinct", explode "$explode  >> "$file_rerun_log"
        
	    echo >> "$file_rerun_log"
	    
        echo "from bash: rerun "$rerun", extinct "$extinct", explode "$explode
	
	    #mv ${data_dir_fin}/realizations/realization_${r}/*.pdf ${IC_dir_fin}/.
	    
	    echo
	    
    done
 
    #echo
    #echo "zipping"
    #cd "${IC_dir_fin}/.."
    #pwd
    #rm -rf "${param_dir}.tar.gz"
    #tar -zcvf "${param_dir}.tar.gz" "$param_dir"
    #rm -rf "$param_dir"
    #
    
    echo
    
    ### SAVE RECORD OF EXTINCTION AND EXPLOSIONS

 
    echo "$tprec"  "$runcount" "$tottime_erased" "$extinctcount" "$explodecount" "$probextinct" "$probexplode" | awk 'BEGIN{OFS="\t"} { print $1, $2, $3 + $1, $4, $5, $6, $7; }' >> $file_rec_events

 
    echo "from bash: rerun "$rerun", extinct "$extinct
 
    
    
    echo "total number of runs "$runcount", total number of extinctions "$extinctcount", probextinct "$probextinct", total number of explosions "$explodecount", probextinct "$probexplode  >> "$file_rerun_log"
    echo "total number of runs "$runcount", total number of extinctions "$extinctcount", probextinct "$probextinct", total number of explosions "$explodecount", probextinct "$probexplode
 
    cd $curr_dir
    echo " "
    echo "ppl_num = "$ppl_num" mu = "$mu" rec_width = "$rec_width" sigma = "$sigma" eps = "$eps" save_cov = "$save_cov", run "$count
            

    echo "END"
    
    ulimit -S -v unlimited

fi




