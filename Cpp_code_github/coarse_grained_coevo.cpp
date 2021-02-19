#include "classi.h"
#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath>
#include <unistd.h> 
#include <algorithm>
#include <regex>
//#include <filesystem>
//namespace fs = std::filesystem;

//simulation of the viral-immune co-evolutionary dynamics.  Take parameters as input, and saves a number of .dat and .txt files to keep track of the systems state as well as of the fitness aproximation errors (see Methods and Supplementary Information)

using namespace std;

int dtorcalls=0;

random_device rd;
mt19937 mt(rd());

//double timescale_inf=3;


int main(int argc, char **argv) {
//cerr<<argc<<endl;
	
	if(argc!=17){
		cerr<<"Usage "<<argv[0]<<" <mu> <recognition width> <jump_size> <people number>  <out_dir>  <simulation time> <full_frame save time> <n_reals> <real> <initial_condition> <save_final_configuration> <kernel> <save_phylo> <fake_in_cond_str> <mem_points> <F0>   " <<endl;
		return -1;
	}
	
    int mem_points=atoi(argv[15]); // 
    
    //double F0 = 1.;
    double F0 =atof(argv[16]);
    
    
    
    bool phylo_subsample=true; 
	

    double mu=atof(argv[1]); //virus mutation rate 
    
    double t_I=1; // infection timescale
    
    double recog_width=atof(argv[2]);  // recognition characteristic distance in cross reactivity
    double latt_width=1.; // lattice width, 
    double jump_size=atof(argv[3]); //  characteristic jump distanc for viral mutations
    
    double L=100; 
    
    long int people_number=atol(argv[4]); // number of hosts 
    int DIM=2; // dimension of antigenic space

    int maxtime=atoi(argv[6]);// maximum experiment time, in cycles 
    int save_full_time=atoi(argv[7]); 	// Every how much time should I print viruses snapshots
    int save_time=10;	// Every how much time should I print avg stats
    
    int save_full_IS_time=save_full_time*100; // Every how much time should I print receptors snapshots
    int save_full_simulation=10000; // every how many cycles save full simulation status as checkpoint

    
    
    int n_real=atoi(argv[8]);
    int real=atoi(argv[9]); // cycle on realizations on an external loop
 
    cout << people_number << endl;
 
 
 
 
 
    string initial_condition (argv[10]); // load file with initial condition??
    bool in_cond;
    if (initial_condition=="true"){
	in_cond=true;
    }
    else if (initial_condition=="false"){
	in_cond=false;
    }
    
    
    else{
	    cout << "wrong string for in_cond, either false or true!" << endl;
	    exit(1);
	
    }
    cout << in_cond << endl;
    
 
 
    string save_final_configuration (argv[11]); // save final configuration??
    bool save_fin;
    if (save_final_configuration=="true"){
	save_fin=true;
    }
    else if (save_final_configuration=="false"){
	save_fin=false;
    }
    
    
    else{
	    cout << "wrong string for save_final_configuration, either false or true!" << endl;
	    exit(1);
	
    }
    cout << save_fin << endl;
    
    
    string save_vir_phylo (argv[13]); //  this is outdated, disregard
    bool vir_phylo;
    if (save_vir_phylo=="true"){
	vir_phylo=true;
    }
    else if (save_vir_phylo=="false"){
	vir_phylo=false;
    }
    
    
    else{
	    cout << "wrong string for save_final_configuration, either false or true!" << endl;
	    exit(1);
	
    }
    cout << vir_phylo << endl;
    
    string fake_in_cond_str (argv[14]); // outdated
    bool fake_in_cond;
    if (fake_in_cond_str=="true"){
	fake_in_cond=true;
    }
    else if (fake_in_cond_str=="false"){
	fake_in_cond=false;
    }
    
    
    else{
	    cout << "wrong string for fake_in_cond_str, either false or true!" << endl;
	    exit(1);
	
    }
    cout << fake_in_cond << endl;
    
    string kernell (argv[12]); // shape of recognition kernel. Implemented and debugged for exponential
    cout << kernell << endl;
    
    // files and dirs
    string out_dir (argv[5]); // output directory 
    
    // ---------------------------------------------------------------------------------  
    // ------------------ HANDLE OUTPUT DIRECTORY AND PARAMETERS BACKUP  ------------------   
    // ---------------------------------------------------------------------------------  
    
        
    //string mu_str=to_string(mu);
    string mu_str=to_string_with_precision(mu,10);
    mu_str.erase ( mu_str.find_last_not_of('0') + 1, string::npos );
    
    
    string F0_str=to_string(F0);
    F0_str.erase ( F0_str.find_last_not_of('0') + 1, string::npos );
    
    
    
    string recog_width_str=to_string(recog_width);
    recog_width_str.erase ( recog_width_str.find_last_not_of('0') + 1, string::npos );
    
    string latt_width_str=to_string(jump_size);
    latt_width_str.erase ( latt_width_str.find_last_not_of('0') + 1, string::npos );
    
    string people_number_str=to_string(people_number);
    string DIM_str=to_string(DIM);
    string maxinfections_str=to_string(maxtime);

    
    string param_dir ("D_" + DIM_str +"_pplnum_"+ people_number_str +"_maxt_"+ maxinfections_str +"_mu_"+ mu_str + "_rec_width_"+ recog_width_str +"_jump_size_"+ latt_width_str  +"_F0_"+ F0_str  +"/"); // output directory 

    replace( param_dir.begin(), param_dir.end(), '.', 'd'); 

// 
    
    out_dir = out_dir + param_dir;


//    in_cond_dir = in_cond_dir + param_dir_IC; // created by shell script

   
//--------------------------------------------------------------------------------------------DEALS WITH OUTPUT FILES   
    
    struct stat sb;

    string out_dir_parbackup=out_dir;
    

    string real_dir ("realizations/");

    out_dir = out_dir + real_dir;
    
 
    string out_dir_bas=out_dir;
    
    //---------------------------------------CYCLE ON REALIZATIONS-----------------
    
    //for(int real=1; real <= n_real; real++){
    
        string real_str=to_string(real);
	
	string real_n_dir ("realization_" + real_str +"/");
    
	out_dir = out_dir_bas + real_n_dir;
	
	cout<< out_dir << endl;
	
    if(! in_cond){
        if ((stat(out_dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
            cout << "out dir already exists, proceding with the deletion" << endl;
            system(("rm -rf "+out_dir).c_str());// linux command
        }
        
        if (!(stat(out_dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
        {
            int status;
            
            cout << "out dir does not exists" << endl;
        
            status = mkdir(out_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == status)
            {
            printf("Error creating directory! \n");
            exit(1);
            }
            else cout << "out dir succesfully created" << endl;
        }
        else{
            cout << "out dir already exists, deletion didn't work!" << endl;
            exit(1);
        }	
        
    }
    
    string fr ("frames/");
    
    string out_dir_frames = out_dir + fr;

    if (!(stat(out_dir_frames.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
    {
        int status;
        
        cout << "out dir frames does not exists" << endl;
    
        status = mkdir(out_dir_frames.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == status)
        {
        printf("Error creating directory frames! \n");
        exit(1);
        }
        else cout << "out dir frames succesfully created" << endl;
    }
    
	string saves_f ("saves_backup.txt");
	
	fstream saves_file;
	saves_file.open(out_dir + saves_f, fstream::out); 
	
	saves_file<<"# 1  kernel  2 save_phylo"<<endl;
	
	saves_file <<setw(19) <<  kernell  <<setw(19) <<  save_vir_phylo  << endl;


	 
        string evo_mean_f ("evo_mean_stats_real_" + real_str +".dat");
	
        //string traj_mean_f ("traj_mean_stats_real_" + real_str +".dat");
	
        string viral_phylogeny_f ("viral_phylogeny_real_" + real_str +".dat");
  
	fstream evo_mean_file;
	
	// fstream traj_mean_file;
	
	fstream viral_phylogeny_file;
    
    
    
    
    string FFT_e_f ("FFT_errors.dat");
    fstream FFT_e_file;
    
    string FFT_farfield_e_f ("FFT_comb_farfield_errors.dat");
    fstream FFT_farfield_e_file;
    
    string FFT_farfield_theta_e_f ("FFT_comb_farfield_theta_errors.dat");
    fstream FFT_farfield_theta_e_file;
    
    string farfield_e_f ("farfield_errors.dat");
    fstream farfield_e_file;
    

      
	    evo_mean_file.open(out_dir + evo_mean_f, fstream::out | fstream::app); 
            
            
        
	    FFT_e_file.open(out_dir + FFT_e_f, fstream::out | fstream::app); 
	    FFT_farfield_e_file.open(out_dir + FFT_farfield_e_f, fstream::out | fstream::app); 
	    FFT_farfield_theta_e_file.open(out_dir + FFT_farfield_theta_e_f, fstream::out | fstream::app); 
	    farfield_e_file.open(out_dir + farfield_e_f, fstream::out | fstream::app); 
	    //evo_mean_file.open(out_dir + evo_mean_f, fstream::app); 
	    //evo_mean_file.open(out_dir + evo_mean_f, fstream::out); 
	    
	    cout << out_dir + evo_mean_f << endl;
	    //evo_mean_file <<"# prova"<<endl;
	    
	    //traj_mean_file.open(out_dir + traj_mean_f, fstream::out | fstream::app); 
	    
	    if (vir_phylo){ 
		//viral_phylogeny_file.open(out_dir + viral_phylogeny_f, fstream::app); 
		viral_phylogeny_file.open(out_dir + viral_phylogeny_f, fstream::out | fstream::app); 
//		if(in_cond && phylo_subsample){
//		    viral_phylogeny_file<<"# 1 ID 1"<<setw(30)<<" 2 ID 2 " <<setw(30) <<" 3 branch_length "<<setw(30)<<" 4 time 1->2"<<setw(30)<<"5 1_x"<<setw(30)<<"6 1_y"<<setw(30)<<"7 2_x"<<setw(30)<<"8 2_y"<<endl;
//		}
	    }

	    
	    
	if(!in_cond){
	    
	    cout << "not incond" << endl;
            
	    
	    //evo_mean_file.open(out_dir + evo_mean_f, fstream::out); 
	    
	    //evo_mean_file<<"# 1 infection"<<setw(30)<<" 2 time "<<setw(30)<<" 3 num_infected"<<setw(30)<<"4 frac_infected"<<setw(30)<<"5 avg_return_time"<<setw(30)<<"6 avg_infection_time"<<setw(30)<<"7 non_naive"<<setw(30)<<"8 avg_infection_num"<<setw(30)<<"9 temp_avg_passed_infections"<<setw(30)<<"10 temp_avg_succesf_passed_infections"<<setw(30)<<"11 temp_avg_p_f"<<setw(30)<<"12 temp_avg_M"<<setw(30)<<"13 temp_avg_Rnot"<<setw(30)<<"14 num_jump_instant"<<setw(30)<<"15 avg_jump_hist"<<setw(30)<<"16 avg_viral_strain_pop_size"<<setw(30)<<"17 more_1_inf"<<setw(30)<<"18 memory"<<setw(30)<<"19 mem capacity"<<endl;
    
        evo_mean_file<<"# 0 time"<<setw(30)<<" 1 num_x"<<setw(30)<<" 2 num_y "<<setw(30)<<" 3 area"<<setw(30)<<"4 number IS coord"<<setw(30)<<"5 number viruses coord"<<setw(30)<<" 6 tot number IS"<<setw(30)<<"7 tot number viruses"<<setw(30)<<"8 avg viral fitness"<<setw(30)<<"9 var viral fitness" <<setw(30)<<"10 start_x" <<setw(30)<<"11 start_y" <<setw(30)<<"12 avg_vir_x" <<setw(30)<<"13 avg_vir_y"<<setw(30)<<"14 maxfit" <<setw(30)<<"15 number IS coord update" <<setw(30)<<"16 tot number IS update"  <<setw(30)<<"17 min vir x"   <<setw(30)<<"18 min vir y"   <<setw(30)<<"19 max vir x"   <<setw(30)<<"20 max vir y"   <<setw(30)<<"21 min IS x "   <<setw(30)<<"22 min IS y "   <<setw(30)<<"23 max IS x "   <<setw(30)<<"24 max IS y "   <<setw(30)<<"25 mean_displ_x_ "   <<setw(30)<<"26 mean_displ_y_ "   <<setw(30)<<"27 mean_displ_tot_ "   <<setw(30)<<"28 var_displ_x_ "   <<setw(30)<<"29 var_displ_y_ "   <<setw(30)<<"30 count_displ_ "   <<setw(30)<<"31 x1_max_fitn "   <<setw(30)<<"32 y1_max_fitn " <<endl;
        

     
        FFT_e_file<<"# 0 FFT_fitn_maxvir"<<setw(30)<<" 1 full_fitn_maxvir"<<setw(30)<<" 2 abs_diff_fitn_maxvir " <<"# 3 FFT_fitn_maxIS"<<setw(30)<<" 4 full_fitn_maxIS"<<setw(30)<<" 5 abs_diff_fitn_maxIS " <<"# 6 FFT_fitn_maxfit"<<setw(30)<<" 7 full_fitn_maxfit"<<setw(30)<<" 8 abs_diff_fitn_maxfit " <<endl;
        
        FFT_farfield_e_file<<"# 0 FFT_fitn_maxvir"<<setw(30)<<" 1 full_fitn_maxvir"<<setw(30)<<" 2 abs_diff_fitn_maxvir " <<"# 3 FFT_fitn_maxIS"<<setw(30)<<" 4 full_fitn_maxIS"<<setw(30)<<" 5 abs_diff_fitn_maxIS " <<"# 6 FFT_fitn_maxfit"<<setw(30)<<" 7 full_fitn_maxfit"<<setw(30)<<" 8 abs_diff_fitn_maxfit " <<endl;
        FFT_farfield_theta_e_file<<"# 0 FFT_fitn_maxvir"<<setw(30)<<" 1 full_fitn_maxvir"<<setw(30)<<" 2 abs_diff_fitn_maxvir " <<"# 3 FFT_fitn_maxIS"<<setw(30)<<" 4 full_fitn_maxIS"<<setw(30)<<" 5 abs_diff_fitn_maxIS " <<"# 6 FFT_fitn_maxfit"<<setw(30)<<" 7 full_fitn_maxfit"<<setw(30)<<" 8 abs_diff_fitn_maxfit " <<endl;
        
     
        farfield_e_file<<"# 0 farfield_fitn_maxvir"<<setw(30)<<" 1 full_fitn_maxvir"<<setw(30)<<" 2 abs_diff_fitn_maxvir " <<"# 3 farfield_fitn_maxIS"<<setw(30)<<" 4 full_fitn_maxIS"<<setw(30)<<" 5 abs_diff_fitn_maxIS " <<"# 6 farfield_fitn_maxfit"<<setw(30)<<" 7 full_fitn_maxfit"<<setw(30)<<" 8 abs_diff_fitn_maxfit " <<endl;
        
        
     
    
    
            
	    
	    //traj_mean_file.open(out_dir + traj_mean_f, fstream::out); 
	    
	    //traj_mean_file<<"# 1 infection"<<setw(30)<<" 2 virus_avg_x "<<setw(30)<<" 3 virus_avg_y"<<setw(30)<<"4 virus_avg_dist_first"<<setw(30)<<"5 IS_avg_x"<<setw(30)<<"6 IS_avg_y"<<setw(30)<<"7 IS_avg_dist_first"<<setw(30)<<"8 time"<<setw(30)<<"9 avg_viral_dmean"<<endl;
    
	    
	    
	    //viral_phylogeny_file.open(out_dir + viral_phylogeny_f, fstream::out); 
	    if (vir_phylo)

		viral_phylogeny_file<<"# 1 ID 1"<<setw(30)<<" 2 ID 2 " <<setw(30) <<" 3 branch_length "<<setw(30)<<" 4 time 1->2"<<setw(30)<<"5 1_x"<<setw(30)<<"6 1_y"<<setw(30)<<"7 2_x"<<setw(30)<<"8 2_y"<<endl;
		
		
	}
       
 
   


    int infection=0;
    double t=0;
    double t_init=0;
	//bool initiated=false; // if the actual experiment is started, if I want to wait termalization. I need it if I don't want to save final config and restart from that, for example with pop not fixed. Criterium should be to wAIT A LIL LONGER THAT FAKE ic directed jumps
	
    
	
	double prob_mut=(1-exp(-mu*t_I)); // per vir per day
	
	cout << prob_mut <<endl;
	    
        
 
 	int time_ratio_printed           =-1;
 	int time_ratio_IS_printed        =-1;
	int time_ratio_printed_avg       =-1;
	int time_ratio_printed_fullstatus=0;
	//int time_ratio_changed_params=-1;
	int last_echo_infection=0;
	//int last_vir_ID; 
    
    antigenic_space world;

    // ---------------------------------------------------------------------------------   
    //----------------------------------------LOAD INITIAL CONDITION ---------------------------
    // ---------------------------------------------------------------------------------  
    
    
    if(in_cond){

	    cout<<"LOADING INITIAL CONFIGURATION"<<endl;

	    fstream initial_condition_file;


        string initial_condition_dir ("fullstatus_last/"); // 
        
        initial_condition_dir= out_dir + initial_condition_dir;
        
        
        string initial_condition_f ("miscellaneous.dat"); // 

	    initial_condition_file.open(initial_condition_dir + initial_condition_f, fstream::in); 
        
         cout << initial_condition_dir + initial_condition_f << endl;


        double t_prec;
    
        long int numvirtot;
        
        
        int num_x;//grid dimensions, number of entries
        int num_y;
        
        int start_x; // coordinates corresponding to low left corner of current grid, in grid number units
        int start_y;
        
        int min_vir_x; // extremes of viral cloud, to be checked and changed whenever a site is emptied or colonized
        int min_vir_y;
        int max_vir_x;
        int max_vir_y;
        int min_IS_x;
        int min_IS_y;
        int max_IS_x;
        int max_IS_y;
        
        int pad_size; 
            
	    if(initial_condition_file.is_open())
	    {
	    
		
		string line;
		while (getline(initial_condition_file, line))
		{
		    if (line[0] != '#'){
			
			istringstream iss(line);
		    

			if (!(iss >> t_prec >> numvirtot >>  num_x >>  num_y >>  start_x >>  start_y >>  min_vir_x >>  min_vir_y >>  max_vir_x >>  max_vir_y >>  min_IS_x >>  min_IS_y >>  max_IS_x >>  max_IS_y >>  pad_size )) {
			    
			    cout << "Error loading initial_condition_file"<<endl;
			    cout<<t_prec<<" "<< numvirtot<<" "  << num_x<<" " << num_y<<" " << start_x <<" " << start_y <<endl;
			    
			    exit(1); 
			    
			} // error
	    
			    
		    }
		
		}
        }
	    else
	    {
		cout << "Error opening initial_condition_file";
		exit(1);
	    }
        
        initial_condition_file.close();
        
        
    
        
       
        vector<vector<antigenic_coordinate>> w;
        //deque<deque<antigenic_coordinate>> w;
        vector<vector<int>> vir_coords;
        vector<vector<int>> IS_coords;
        vector<vector<long int>> IS_coords_update;
     
        // grid initialization
        antigenic_coordinate empty;
        
        
        vector<antigenic_coordinate> empty_col(num_y, empty);
        
        
        for(int i=0; i < num_x; i++){
            w.push_back(empty_col);
        }    

        //while(true){}; 

        initial_condition_f = ("vir_coords.dat"); //

	    initial_condition_file.open(initial_condition_dir + initial_condition_f, fstream::in); 


        cout<<"before loading virs"<<endl;

	    if(initial_condition_file.is_open())
	    {
	    
		
		string line;
		while (getline(initial_condition_file, line))
		{
		    if (line[0] != '#'){
			
			istringstream iss(line);
		    
    
            int xi;
            int yi;
            long int num_virs;
            double fitness ;
            long int num_IS;
            
        
			if (!(iss >> xi >>  yi >>  num_IS >>  num_virs >>  fitness)) {
			    
			    cout << "Error loading initial_condition_file"<<endl;
			    cout<<xi<<" "<< yi <<endl;
			    
			    exit(1); 
			    
			} // error
            
            vector<int> coord {xi + start_x, yi + start_y};
            vir_coords.push_back(coord);
            
            w[xi][yi].set_num_virs(num_virs);
            w[xi][yi].set_fitness(fitness);
        
		    }
		
		}
		
		cout<<"after LOADING"<<endl;
        //while(true){}; 
	    
	    
	    }
	    else
	    {
		cout << "Error opening initial_condition_file";
		exit(1);
	    }
    
        initial_condition_file.close();


        initial_condition_f = ("IS_coords.dat"); // 

	    initial_condition_file.open(initial_condition_dir + initial_condition_f, fstream::in); 


        cout<<"before loading IS"<<endl;

	    if(initial_condition_file.is_open())
	    {
	    
		
		string line;
		while (getline(initial_condition_file, line))
		{
		    if (line[0] != '#'){
			
			istringstream iss(line);
		    
    
            int xi;
            int yi;
            long int num_IS;
            
        
			if (!(iss >> xi >>  yi >>  num_IS )) {
			    
			    cout << "Error loading initial_condition_file"<<endl;
			    cout<<xi<<" "<< yi <<endl;
			    
			    exit(1); 
			    
			} // error
            
            vector<int> coord {xi + start_x, yi + start_y};
            IS_coords.push_back(coord);
            
            w[xi][yi].set_num_IS(num_IS);
        
		    }
		
		}
		
		cout<<"after LOADING"<<endl;
	    
	    
	    }
	    else
	    {
		cout << "Error opening initial_condition_file";
		exit(1);
	    }
    
        initial_condition_file.close();

    

        initial_condition_f = ("IS_upd_coords.dat"); // 
	    initial_condition_file.open(initial_condition_dir + initial_condition_f, fstream::in); 


        cout<<"before loading IS upd"<<endl;

	    if(initial_condition_file.is_open())
	    {
	    
		
		string line;
		while (getline(initial_condition_file, line))
		{
		    if (line[0] != '#'){
			
			istringstream iss(line);
		    
    
            int xi;
            int yi;
            long int num_IS;
            long int num_IS_upd;
            
        
			if (!(iss >> xi >>  yi >>  num_IS  >>  num_IS_upd )) {
			    
			    cout << "Error loading initial_condition_file"<<endl;
			    cout<<xi<<" "<< yi <<endl;
			    
			    exit(1); 
			    
			} // error
            
            num_IS_upd=0;
            //vector<long int> coord {xi + start_x, yi + start_y, num_IS_upd};
            //IS_coords_update.push_back(coord);
            
		    }
		
		}
		
		cout<<"after LOADING"<<endl;
	    
	    
	    }
	    else
	    {
		cout << "Error opening initial_condition_file";
		exit(1);
	    }
    
        initial_condition_file.close();
        
        cout<<"vir_coords "<<vir_coords.size()<<" vir_coords "<<IS_coords.size()<<" IS_coords_update "<<IS_coords_update.size()<<" IS_coords "<<IS_coords.size()<<" num_x*num_y "<<num_x*num_y<<endl;
        cout<<"w "<<w.size()<<" w[0] "<<w[0].size()<<" start_x "<<start_x <<" start_y "<<start_y<<" num_x "<<num_x<<" num_y "<<num_y<<endl;
           
           
    
        cout<<"w "<<w.size()<<" w[0] "<<w[0].size()<<" start_x "<<start_x <<" start_y "<<start_y<<" num_x "<<num_x<<" num_y "<<num_y<<endl;
        
    
    
        cout<< "w "<<w.size() << "  "<<w.end() - w.begin()<<endl;;
        
        //while(true){}; 
        
    
        world= antigenic_space(vir_coords,  IS_coords, IS_coords_update, DIM, recog_width, latt_width,  jump_size,  people_number,  mem_points,  kernell,  mu,  out_dir,   num_x,   num_y,   start_x,   start_y, w, F0, pad_size, min_vir_x, min_vir_y, max_vir_x, max_vir_y, min_IS_x, min_IS_y, max_IS_x, max_IS_y);// 
        
        cout<<"after ctor"<<endl;
        
        //while(true){};
        
        t+=t_prec; // this has nothing to do with steady state, t_init stays 0
        

                
        time_ratio_printed           =int((t_prec + t_I)/save_full_time);
        time_ratio_IS_printed        =int((t_prec + t_I)/save_full_IS_time);
        time_ratio_printed_avg       =int((t_prec + t_I)/save_time);
        time_ratio_printed_fullstatus=int((t_prec + t_I)/save_full_simulation);

    
    
        
	    cout<<"initial condition initialized"<<endl;
	    
	}
    else{
         cout << "START creation"<< endl;
    
        world = create_world(DIM, recog_width,  latt_width, jump_size,  people_number,  mem_points,  kernell,  mu,  out_dir, F0);
         cout << "END creation"<< endl;

        double maxfiterr_o_std=0;
        world.update_fitness_IS_update_switch(true, maxfiterr_o_std);

        world.print_avg_to_file(evo_mean_file, t);
        world.print_to_file(t);
        
    }
    
    //while(true){};

    //t+= t_inf;      



 
    // ---------------------------------------------------------------------------------  
    // ------------------------ MONTE CARLO SIMULATION  ---------------------------------   
    // ---------------------------------------------------------------------------------  
      


    cout << "START SIMULATION with maximum time " << maxtime <<" REALIZATION = "<< real << endl;
// while(true){
//     
//     
// };        

	string param_f ("parameters_backup.dat");
	
	fstream param_file;
	param_file.open(out_dir_parbackup + param_f, fstream::out); 
    
	param_file<<"# 1 <mu> 	2 <recognition width>	3 <jump_size>	4 <people number>		5  <simulation time>	6 <full_frame save time>	7	 <n_reals>	8 <initial_condition>	9  <fake_initial_condition>	10 <phylo_subsample>	11 <mem_points>	12 <t_init>	13 <F0>"<<endl;
	
	param_file << fixed  << setprecision(5)<<setw(19) << mu <<setprecision(5)<<setw(19) << recog_width <<setprecision(5)<<setw(19) << jump_size <<setprecision(0)<<setw(19) <<  people_number <<setprecision(0)<<setw(19) <<  maxtime + int(t_init) <<setprecision(0)<<setw(19) <<  save_full_time <<setprecision(0)<<setw(19) <<  n_real <<setw(19) <<  initial_condition <<setw(19) <<  fake_in_cond_str <<setw(19) <<  BoolToString(phylo_subsample)<<setprecision(0)<<setw(19) <<  mem_points  <<setprecision(0)<<setw(19) <<  int(t_init)   <<setprecision(5)<<setw(19) <<  F0 << endl;

	param_file.close();
	    
        
        int count_full_upd=0;
        
        //while((t-t_init<maxtime || !initiated) && infected.size() >0)
        while(t-t_init<maxtime && world.get_num_vir_coords() >0)
        {
	   //cout <<"START loop " <<endl;

            infection++;
            bool print_progress=false;

            t+=t_I; 
            
            bool benchmark_update_fitness=(int(t)%10000==1);
            
            
            world.infection_cycle(count_full_upd, benchmark_update_fitness); // THE SHIT HAPPENS HERE

            //double infected_frac = infected.size()/double(people_number);
            
            
            if (world.get_num_vir_coords() > 0){
        
                //double t_next=infected[0].get_t_inf_end();
                double t_next= t+t_I;
                    
             
                
                if(t>t_next){
                    cout<< "Errore nel time ordering "<<endl;
                    exit(1);
                }
                
                int time_ratio_next=int(t_next/save_full_time);
                int time_ratio_IS_next=int(t_next/save_full_IS_time);
                //int time_ratio=int(t/save_full_time);
            
                int time_ratio_next_avg=int(t_next/save_time);
                
                int time_ratio_next_fullstatus=int(t_next/save_full_simulation);
        
            
                     //cout<<"AFTER INFECTION before printing"<< endl;
                        
//---------------------------------------------------PRINT DATA-------------------------------------------------   
                 cout << "printing time  " << t << " of maximum time " << maxtime << " inf " << infection<<endl;
                // Once every save_rate, I save data |
                if(time_ratio_next_avg > time_ratio_printed_avg){  //  
                    if(infection - last_echo_infection>1000 || infection==1) {
                        print_progress=true;
                        last_echo_infection=infection;
                    }
                    
                    time_ratio_printed_avg=time_ratio_next_avg;
                    
                    
                    world.print_avg_to_file(evo_mean_file, t);
                    
                    cout<<"Printed "<<endl;
                }
     
     
     //-----------------------------------------------------------------------------------------------------------------
                
                // Once every save_rate, I save full viral snapshots|
                if(time_ratio_next > time_ratio_printed ){  // 
                    time_ratio_printed=time_ratio_next;
                    cout << "printing GLOBAL infection  " << infection << " time " << t <<endl;
            
            
                    world.print_to_file(t);

                    cout <<"after printing global "<<endl;
            
                }
            
                // Once every save_rate_IS, I save IS subsample |
                if(time_ratio_IS_next > time_ratio_IS_printed ){  //
                    time_ratio_IS_printed=time_ratio_IS_next;
                    cout << "printing GLOBAL IS infection  " << infection << " time " << t <<endl;
            
            
                    world.print_to_file_IS_update(t);

                    cout <<"after printing global "<<endl;
            
                }
            
            
                // Once every save_full_simulation, I save checkpoint |
                if(time_ratio_next_fullstatus > time_ratio_printed_fullstatus  ){  // 
                    
                    bool first_save= (time_ratio_printed_fullstatus==0);
                     
                    time_ratio_printed_fullstatus=time_ratio_next_fullstatus;
                    cout << "printing FULL status, infection  " << infection << " time " << t <<endl;
            
                    if (world.get_numvirtot() > 1000 && world.get_numvirtot() < people_number/10)
                        world.print_to_file_fullstatus(first_save, t);

                    cout <<"after printing global "<<endl;
            
                }
            }
                 //cout<<"AFTER INFECTION after printing tot "<< infection << " "<< maxinfections<< " "<<infected.size()<< endl;
                 //cout<<"AFTER INFECTION after printing tot "<< (infection<maxinfections)  << " "<< (infected.size() >0) << " "<<(infection<maxinfections && infected.size() >0) << endl;
        }
        
        cout<<"AFTER all infections "<<  endl;
	
        
        FFT_e_file.close();
        FFT_farfield_e_file.close();
        FFT_farfield_theta_e_file.close();
        farfield_e_file.close();
        evo_mean_file.close();
        if (vir_phylo)
            viral_phylogeny_file.close();
    
	cout<<"END "<<  endl;

	return 0;
}

