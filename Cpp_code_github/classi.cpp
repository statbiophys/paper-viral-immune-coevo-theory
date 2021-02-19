#include "classi.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <numeric>
#include "fastsum_NFFT_simple.h"



using namespace std;

antigenic_coordinate::antigenic_coordinate() {


    
    num_virs_=0;    
    num_IS_=0;    
    fitness_=0.;    




};

antigenic_coordinate::antigenic_coordinate(long int num_virs, long int num_IS, double fitness) {


	
	num_virs_=num_virs;
	num_IS_=num_IS;
	fitness_=fitness;
   
};

antigenic_coordinate::antigenic_coordinate(const antigenic_coordinate & c) {
 
   
   
	//coord_=c.coord_;
	//first_nns_=c.first_nns_;
	num_virs_=c.num_virs_;
	num_IS_=c.num_IS_;
	fitness_=c.fitness_;


};

antigenic_coordinate & antigenic_coordinate::operator= (const antigenic_coordinate& c){ 


	num_virs_=c.num_virs_;
	num_IS_=c.num_IS_;
	fitness_=c.fitness_;

	return *this;

};




antigenic_space::antigenic_space () {
    //cout << "entered constructor antigenic_space "  << endl;
    
    DIM_           = 0;
    vir_coords_    = vector<vector<int>>();
    IS_coords_     = vector<vector<int>>();
    IS_update_coords_     = vector<vector<long int>>();
    recog_width_   = 0.;
    latt_width_    = 0.;
    jump_size_    = 0.;
    people_number_ = 0;
    mem_points_    = 0;
    F0_            = 0;
    kernel_        = "";
    mu_            = 0;
    out_dir_       = "";
    num_x_         = 0;
    num_y_         = 0;
    start_x_         = 0;
    start_y_         = 0;
    //world_         = deque<deque<antigenic_coordinate>>();
    world_         = vector<vector<antigenic_coordinate>>();
    
    min_vir_x_         = 0;
    min_vir_y_         = 0;
    max_vir_x_         = 0;
    max_vir_y_         = 0;
    min_IS_x_          = 0;
    min_IS_y_          = 0;
    max_IS_x_          = 0;
    max_IS_y_          = 0;
    
    pad_size_          = 0 ;
    
    numvirtot_         = 0;
    
    displ_tot_    =0;
    displ_x_    =0;
    displ_y_    =0;
    displ_x_sq_ =0;
    displ_y_sq_ =0;
    count_displ_=0;
    
    
    
    
        //std::cout << "Host ctor @" << this << std::endl; 

};


antigenic_space::antigenic_space (vector<vector<int>> & vir_coords,  vector<vector<int>> & IS_coords,  vector<vector<long int>> & IS_update_coords, int DIM, double recog_width, double latt_width, double jump_size, long int people_number, int mem_points, string kernell, double mu, string out_dir,  int num_x,  int num_y,  int start_x,  int start_y, vector<vector<antigenic_coordinate>> & world, double F0, int pad_size, int min_vir_x, int min_vir_y, int max_vir_x, int max_vir_y, int min_IS_x, int min_IS_y, int max_IS_x, int max_IS_y) {
    //cout << "entered constructor antigenic_space "  << endl;
    
    cout<<"constructor of copies"<<endl;
    
    DIM_           = DIM;
    vir_coords_    = vir_coords;
    IS_coords_     = IS_coords;
    IS_update_coords_     = IS_update_coords;
    recog_width_   = recog_width;
    latt_width_    = latt_width;
    if(latt_width_!=1){
        cout << "latt_width_ not 1! "<< latt_width_ << endl;
        exit(1); 
    }        
    jump_size_    = jump_size;
    people_number_ = people_number;
    mem_points_    = mem_points;
    F0_            = F0;
    kernel_        = kernell;
    mu_            = mu;
    out_dir_       = out_dir;
    num_x_         = num_x;
    num_y_         = num_y;
    start_x_         = start_x;
    start_y_         = start_y;
    
    cout<<world.size()<<endl;
    cout<<world[0].size()<<endl;
    world_         = world;
    //deque<deque<antigenic_coordinate>>().swap( world)       ;
    vector<vector<antigenic_coordinate>>().swap( world)       ;
    cout<<world_.size()<<endl;
    cout<<world_[0].size()<<endl;

    
    min_vir_x_         = min_vir_x;
    min_vir_y_         = min_vir_y;
    max_vir_x_         = max_vir_x;
    max_vir_y_         = max_vir_y;
    min_IS_x_          = min_IS_x ;
    min_IS_y_          = min_IS_y ;
    max_IS_x_          = max_IS_x ;
    max_IS_y_          = max_IS_y ;
    
    pad_size_          = pad_size ;
    
    numvirtot_         = get_tot_vir();
    
    displ_tot_    =0;
    displ_x_    =0;
    displ_y_    =0;
    displ_x_sq_ =0;
    displ_y_sq_ =0;
    count_displ_=0;
    
    
    
        cout<<"out of ctor"<<endl;

        //std::cout << "Host ctor @" << this << std::endl; 

};

antigenic_space::antigenic_space(const antigenic_space & h){
  
    DIM_          = h.DIM_          ;
    vir_coords_   = h.vir_coords_   ;
    IS_coords_    = h.IS_coords_    ;
    IS_update_coords_    = h.IS_update_coords_    ;
    recog_width_  = h.recog_width_  ;
    latt_width_   = h.latt_width_   ;
    jump_size_   = h.jump_size_   ;
    people_number_= h.people_number_;
    mem_points_   = h.mem_points_   ;
    F0_            = h.F0_   ;
    kernel_       = h.kernel_       ;
    mu_           = h.mu_           ;
    out_dir_      = h.out_dir_      ;
    num_x_        = h.num_x_         ;
    num_y_        = h.num_y_         ;
    start_x_      = h.start_x_      ;
    start_y_      = h.start_y_      ;
    world_        = h.world_;
    numvirtot_        = h.numvirtot_;
    
    min_vir_x_        = h.min_vir_x_;
    min_vir_y_        = h.min_vir_y_;
    max_vir_x_        = h.max_vir_x_;
    max_vir_y_        = h.max_vir_y_;
    min_IS_x_         = h.min_IS_x_ ;
    min_IS_y_         = h.min_IS_y_ ;
    max_IS_x_         = h.max_IS_x_ ;
    max_IS_y_         = h.max_IS_y_ ;
    
    pad_size_         = h.pad_size_ ;
    
    displ_tot_             = h.displ_tot_     ;
    displ_x_             = h.displ_x_     ;
    displ_y_             = h.displ_y_     ;
    displ_x_sq_          = h.displ_x_sq_  ;
    displ_y_sq_          = h.displ_y_sq_  ;
    count_displ_         = h.count_displ_ ;
  
};


antigenic_space & antigenic_space::operator= (const antigenic_space& h){ 
  
    DIM_          = h.DIM_          ;
    vir_coords_   = h.vir_coords_   ;
    IS_coords_    = h.IS_coords_    ;
    IS_update_coords_    = h.IS_update_coords_    ;
    recog_width_  = h.recog_width_  ;
    latt_width_   = h.latt_width_   ;
    jump_size_   = h.jump_size_   ;
    people_number_= h.people_number_;
    mem_points_   = h.mem_points_   ;
    F0_           = h.F0_   ;
    kernel_       = h.kernel_       ;
    mu_           = h.mu_           ;
    out_dir_      = h.out_dir_      ;
    num_x_        = h.num_x_;
    num_y_        = h.num_y_;
    start_x_      = h.start_x_;
    start_y_      = h.start_y_;
    world_        = h.world_;
    numvirtot_        = h.numvirtot_;

 
    min_vir_x_        = h.min_vir_x_;
    min_vir_y_        = h.min_vir_y_;
    max_vir_x_        = h.max_vir_x_;
    max_vir_y_        = h.max_vir_y_;
    min_IS_x_         = h.min_IS_x_ ;
    min_IS_y_         = h.min_IS_y_ ;
    max_IS_x_         = h.max_IS_x_ ;
    max_IS_y_         = h.max_IS_y_ ;
    
    pad_size_         = h.pad_size_ ;
    
    displ_tot_             = h.displ_tot_     ;
    displ_x_             = h.displ_x_     ;
    displ_y_             = h.displ_y_     ;
    displ_x_sq_          = h.displ_x_sq_  ;
    displ_y_sq_          = h.displ_y_sq_  ;
    count_displ_         = h.count_displ_ ;


         //std::cout << "Host ctor @" << this << std::endl; 
   
	return *this;

};



    

long int antigenic_space::get_tot_IS() const
{
    long int tot_IS=0;
       
    for(int i=0; i < IS_coords_.size(); i++){
        int x =  IS_coords_[i][0] - start_x_;
        int y =  IS_coords_[i][1] - start_y_;
        
        tot_IS+=world_[x][y].get_num_IS() ;
    };
    
    return tot_IS;
    
};


    

long int antigenic_space::get_tot_IS_update() const //
{
    long int tot_IS=0;
    long int tot_IS_abs=0;
       
    for(int i=0; i < IS_update_coords_.size(); i++){
        
        tot_IS+=IS_update_coords_[i][2];
        tot_IS_abs+= abs(IS_update_coords_[i][2]); // 
    };
    
    
    if (tot_IS!=0){ // NUM IS CONSEVED
        cout << "nonzero update IS! "<< tot_IS << endl;
        exit(1); 
    }
    
    return tot_IS_abs;
    
};




long int antigenic_space::get_tot_vir() const
{
    long int tot_vir=0;
       
    for(int i=0; i < vir_coords_.size(); i++){ 
        int x =  vir_coords_[i][0] - start_x_;
        int y =  vir_coords_[i][1] - start_y_;
        
        tot_vir+=world_[x][y].get_num_vir() ;
    };
    
    return tot_vir;
    
};




vector<double> antigenic_space::get_avg_fitness() const
{
    double avg_fitness=0.;
    double avg_fitness_sq=0.;
    double maxfit=-DBL_MAX;
    
    long int tot_vir=0;
       
    for(int i=0; i < vir_coords_.size(); i++){ 
        int x =  vir_coords_[i][0] - start_x_;
        int y =  vir_coords_[i][1] - start_y_;
        
        tot_vir+=world_[x][y].get_num_vir() ;
                
        avg_fitness   +=(world_[x][y].get_fitness() )*(world_[x][y].get_num_vir());
        avg_fitness_sq+=pow((world_[x][y].get_fitness() ),2.)*(world_[x][y].get_num_vir());
        
        if (world_[x][y].get_fitness() > maxfit) maxfit=world_[x][y].get_fitness();
    };
    
    vector<double> fitn_stat { avg_fitness/tot_vir, avg_fitness_sq/tot_vir - pow(avg_fitness/tot_vir, 2.), maxfit};
    return fitn_stat;
    
};


int antigenic_space::get_maxobs(int obs) const
{
    double max=-DBL_MAX;
    int idx_max= numeric_limits<int>::min();
    
    
       
    for(int i=0; i < vir_coords_.size(); i++){ 
        int x =  vir_coords_[i][0] - start_x_;
        int y =  vir_coords_[i][1] - start_y_;
        
        if (obs == 2){
            if (world_[x][y].get_fitness() > max){ 
                max=world_[x][y].get_fitness();
                idx_max=i;
            }
        }
        else if (obs == 1){
            if (world_[x][y].get_num_IS() > max){ 
                max=world_[x][y].get_num_IS();
                idx_max=i;
            }
        }
        else if (obs == 0){
            if (world_[x][y].get_num_vir() > max){ 
                max=world_[x][y].get_num_vir();
                idx_max=i;
            }
        }
        else{
            cout << "obs not implemented! "<< obs << endl;
            exit(1); 
        }
    };
    
    return idx_max;
    
};


vector<int> antigenic_space::get_maxobs() const
{
    long int max_nv=numeric_limits<long int>::min();
    long int max_nis=numeric_limits<long int>::min();
    double max_f=-DBL_MAX;
    int idx_max_nv= numeric_limits<int>::min();
    int idx_max_nis= numeric_limits<int>::min();
    int idx_max_f= numeric_limits<int>::min();
    
    
       
    for(int i=0; i < vir_coords_.size(); i++){ 
        int x =  vir_coords_[i][0] - start_x_;
        int y =  vir_coords_[i][1] - start_y_;
        
            if (world_[x][y].get_fitness() > max_f){ 
                max_f=world_[x][y].get_fitness();
                idx_max_f=i;
            }
            if (world_[x][y].get_num_IS() > max_nis){ 
                max_nis=world_[x][y].get_num_IS();
                idx_max_nis=i;
            }
            if (world_[x][y].get_num_vir() > max_nv){ 
                max_nv=world_[x][y].get_num_vir();
                idx_max_nv=i;
            }
    };
    
    vector<int> idxs_max { idx_max_nv, idx_max_nis, idx_max_f};
    return idxs_max;
    
};


vector<double> antigenic_space::get_avg_vir_coord() const
{
    double avg_x=0.;
    double avg_y=0.;
    
    long int tot_vir=0;
       
    for(int i=0; i < vir_coords_.size(); i++){ 
        int x =  vir_coords_[i][0] - start_x_;
        int y =  vir_coords_[i][1] - start_y_;
        
        tot_vir+=world_[x][y].get_num_vir() ;
                
        avg_x+=idx_to_real_x(x)*(world_[x][y].get_num_vir());
        avg_y+=idx_to_real_y(y)*(world_[x][y].get_num_vir());
    }; 
    

    vector<double> avg_coords { avg_x/tot_vir, avg_y/tot_vir };
    return avg_coords;
        
};



void antigenic_space::print_to_file(double t) const
{


    string fr ("frames/");
    
    
    string t_str=to_string(int(t));
    //t_str.erase ( t_str.find_last_not_of('0') + 1, string::npos );
    replace( t_str.begin(), t_str.end(), '.', 'd'); 
    
    string antigenic_space_f ("antigenic_space_time_"+ t_str +".dat");

    string outfile= out_dir_+ fr + antigenic_space_f;

    fstream outstream;
    
    outstream.open(outfile, fstream::out); 
    outstream<<"# 1 x"<<setw(30)<<" 2 y "<<setw(30)<<" 3 number IS"<<setw(30)<<"4 number viruses"<<setw(30)<<"5 viral fitness" <<endl;
            
    for(int i=0; i < vir_coords_.size(); i++){ 
        

        int xi= vir_coords_[i][0] - start_x_;
        int yi= vir_coords_[i][1] - start_y_;
            
			
        long int num_virs   = world_[xi][yi].get_num_vir();
        long int num_IS     = world_[xi][yi].get_num_IS();
        double fitness = world_[xi][yi].get_fitness();
        
        double x= idx_to_real_x(xi);
        double y= idx_to_real_y(yi);

        outstream << fixed  << setprecision(5)<<setw(19) << x <<setprecision(5)<<setw(19) << y <<setprecision(0)<<setw(19) << num_IS <<setprecision(0)<<setw(19) <<  num_virs <<setprecision(5)<<setw(19) <<  fitness <<endl; 
    }



    outstream.close();

};


void antigenic_space::print_to_file_IS_update(double t) const // IS subsample
{
    
    string fr ("frames/");


    string t_str=to_string(int(t));
    //t_str.erase ( t_str.find_last_not_of('0') + 1, string::npos );
    replace( t_str.begin(), t_str.end(), '.', 'd'); 
    
    string antigenic_space_f ("IS_update_antigenic_space_time_"+ t_str +".dat");

    string outfile= out_dir_+ fr + antigenic_space_f;

    fstream outstream;
    
    outstream.open(outfile, fstream::out); 
    outstream<<"# 1 x"<<setw(30)<<" 2 y "<<setw(30)<<" 3 number IS"<<setw(30)<<"4 number viruses"<<setw(30)<<"5 viral fitness" <<setw(30)<<"6 number IS update" <<endl;
            
    for(int i=0; i < IS_update_coords_.size(); i++){ 
        

        long int xi= IS_update_coords_[i][0] - start_x_;
        long int yi= IS_update_coords_[i][1] - start_y_;
        long int num_IS_upd= IS_update_coords_[i][2];
            
			
        long int num_virs   = world_[xi][yi].get_num_vir();
        long int num_IS     = world_[xi][yi].get_num_IS();
        double fitness = world_[xi][yi].get_fitness();
        
        double x= idx_to_real_x(xi);
        double y= idx_to_real_y(yi);

        outstream << fixed  << setprecision(5)<<setw(19) << x <<setprecision(5)<<setw(19) << y <<setprecision(0)<<setw(19) << num_IS <<setprecision(0)<<setw(19) <<  num_virs <<setprecision(5)<<setw(19) <<  fitness  <<setprecision(0)<<setw(19) <<  num_IS_upd <<endl; 
    }


//    
    outstream.close();

};

void antigenic_space::print_to_file_fullstatus(bool first_save, double t) const
{
    
    string fullstatus_last_dir     ("fullstatus_last/");
    string fullstatus_sec_last_dir ("fullstatus_sec_last/");
    
    fullstatus_last_dir    =out_dir_+ fullstatus_last_dir;
    fullstatus_sec_last_dir=out_dir_+ fullstatus_sec_last_dir;
    
    
    if(!first_save){
            
        
        struct stat sb_sec;

        if ((stat(fullstatus_sec_last_dir.c_str(), &sb_sec) == 0 && S_ISDIR(sb_sec.st_mode))){
            cout << "fullstatus_sec_last_dir already exists, proceding with the deletion" << endl;
            system(("rm -rf "+fullstatus_sec_last_dir).c_str());//require myrm.sh to be installed in bin, move to trash
        }
        
        system(("cp -r "+ fullstatus_last_dir +" "+ fullstatus_sec_last_dir ).c_str());//not portable to windows
        
        
    }
 
    
    struct stat sb;

    if ((stat(fullstatus_last_dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
        cout << "fullstatus_last_dir already exists, proceding with the deletion" << endl;
        system(("rm -rf "+fullstatus_last_dir).c_str());//require myrm.sh to be installed in bin, move to trash
    }
    
    if (!(stat(fullstatus_last_dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
    {
        int status;
        
        cout << "out dir does not exists" << endl;
    
        status = mkdir(fullstatus_last_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
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
    
 
    
    string vir_coords_f ("vir_coords.dat");

    string outfile= fullstatus_last_dir+ vir_coords_f;

    fstream outstream_vir;
    
    outstream_vir.open(outfile, fstream::out); 
    //outstream<<"# 1 x"<<setw(30)<<" 2 y "<<setw(30)<<" 3 number IS"<<setw(30)<<"4 number viruses"<<setw(30)<<"5 viral fitness" <<endl;
             
    for(int i=0; i < vir_coords_.size(); i++){ 
        

        int xi= vir_coords_[i][0] - start_x_;
        int yi= vir_coords_[i][1] - start_y_;
            
			
        long int num_virs   = world_[xi][yi].get_num_vir();
        long int num_IS     = world_[xi][yi].get_num_IS();
        double fitness = world_[xi][yi].get_fitness();
        
        //double x= idx_to_real_x(xi);
        //double y= idx_to_real_y(yi);
        
        int x= xi;
        int y= yi;

        outstream_vir << fixed  << setprecision(0)<<setw(19) << x <<setprecision(0)<<setw(19) << y <<setprecision(0)<<setw(19) << num_IS <<setprecision(0)<<setw(19) <<  num_virs <<setprecision(14)<<setw(19) <<  fitness <<endl; 
    }
    
    outstream_vir.close();


    string IS_coords_f ("IS_coords.dat");

    outfile= fullstatus_last_dir+ IS_coords_f;

    fstream outstream_IS;
    
    outstream_IS.open(outfile, fstream::out); 
             
    for(int i=0; i < IS_coords_.size(); i++){
        int x =  IS_coords_[i][0] - start_x_;
        int y =  IS_coords_[i][1] - start_y_;
        long int num_IS =world_[x][y].get_num_IS() ;
    
        outstream_IS << fixed  << setprecision(0)<<setw(19) << x <<setprecision(0)<<setw(19) << y <<setprecision(0)<<setw(19) << num_IS <<endl; 
    }
    outstream_IS.close();



    string IS_upd_coords_f ("IS_upd_coords.dat");

    outfile= fullstatus_last_dir+ IS_upd_coords_f;

    fstream outstream_IS_upd;
    
    outstream_IS_upd.open(outfile, fstream::out); 
             
    for(int i=0; i < IS_update_coords_.size(); i++){ 
        

        long int xi= IS_update_coords_[i][0] - start_x_;
        long int yi= IS_update_coords_[i][1] - start_y_;
        long int num_IS_upd= IS_update_coords_[i][2];
        long int num_IS     = world_[xi][yi].get_num_IS();
        
        int x= xi;
        int y= yi;

        outstream_IS_upd << fixed  << setprecision(0)<<setw(19) << x <<setprecision(0)<<setw(19) << y <<setprecision(0)<<setw(19) << num_IS  <<setprecision(0)<<setw(19) <<  num_IS_upd <<endl; 
    }
    outstream_IS_upd.close();






    string misc_f ("miscellaneous.dat");

    outfile= fullstatus_last_dir+ misc_f;

    fstream outstream_misc;
    
    outstream_misc.open(outfile, fstream::out); 
    
    
    outstream_misc << fixed  << setprecision(5)<<setw(19) << t <<setprecision(0)<<setw(19) << numvirtot_ <<setprecision(0)<<setw(19) << num_x_  <<setprecision(0)<<setw(19) <<  num_y_  <<setprecision(0)<<setw(19) << start_x_ <<setprecision(0)<<setw(19) << start_y_  <<setprecision(0)<<setw(19) <<  min_vir_x_ <<setprecision(0)<<setw(19) << min_vir_y_ <<setprecision(0)<<setw(19) << max_vir_x_  <<setprecision(0)<<setw(19) <<  max_vir_y_  <<setprecision(0)<<setw(19) <<  min_IS_x_ <<setprecision(0)<<setw(19) << min_IS_y_ <<setprecision(0)<<setw(19) << max_IS_x_  <<setprecision(0)<<setw(19) <<  max_IS_y_ <<setprecision(0)<<setw(19) <<  pad_size_ <<endl; 

    outstream_misc.close();

    
    if(first_save){
        system(("cp -r "+ fullstatus_last_dir +" "+ fullstatus_sec_last_dir ).c_str());//not portable to windows
    }
 

};


void antigenic_space::print_avg_to_file(fstream& outstream, double t)
{


            
    int area=0;
    for(int i=0; i < world_.size(); i++){ 
        area+=world_[i].size();
    }
    
    if (area!=num_x_*num_y_){
        cout << "num_y*num_x different than grid area! "<< num_x_ << " * " << num_y_ << " area " << area << endl;
        exit(1); 
    }
    
	long int tot_num_IS_upd = get_tot_IS_update();  
    
    
	long int tot_num_IS = get_tot_IS();  
    if (tot_num_IS!=people_number_*mem_points_){
        cout << "people_number_*mem_points_ different than tot_num_IS! "<< people_number_ << " * " << mem_points_ << " tot_num_IS " << tot_num_IS << endl;
        exit(1); 
    }
    
    
      
	long int tot_num_vir = get_tot_vir();    
    if (tot_num_vir!=numvirtot_){
        cout << "inconsistent numvirtot_! "<< numvirtot_ << " vs " << tot_num_vir << endl;
        exit(1); 
    }
     
	vector<double> fitn_stat = get_avg_fitness();    
	double avg_fitness = fitn_stat[0];    
	double var_fitness = fitn_stat[1];    
	double maxfit = fitn_stat[2];    
     
	vector<double> avg_coords = get_avg_vir_coord();    
	double avg_x = avg_coords[0];    
	double avg_y = avg_coords[1];    
    
    
    double mean_displ_tot_ =displ_tot_/float(count_displ_);
    double mean_displ_x_ =displ_x_/float(count_displ_);
    double mean_displ_y_ =displ_y_/float(count_displ_);
    double var_displ_x_  =displ_x_sq_/float(count_displ_) - pow(mean_displ_x_,2.);
    double var_displ_y_  =displ_y_sq_/float(count_displ_) - pow(mean_displ_y_,2.);
    
    
    vector<int> idxs_max= get_maxobs(); //num virs, IS, fitn
    int idx_max_vir = idxs_max[2];
    int x1_max_fitn =  vir_coords_[idx_max_vir][0];
    int y1_max_fitn =  vir_coords_[idx_max_vir][1];
    
    

    
    outstream << fixed  << setprecision(5)<<setw(19) << t << setprecision(0)<<setw(19) << num_x_ <<setprecision(0)<<setw(19) << num_y_ <<setprecision(0)<<setw(19) << area <<setprecision(0)<<setw(19) <<  IS_coords_.size() <<setprecision(0)<<setw(19) <<  vir_coords_.size()  <<setprecision(0)<<setw(19) <<  tot_num_IS <<setprecision(0)<<setw(19) << tot_num_vir <<setprecision(7)<<setw(19) <<  avg_fitness <<setprecision(7)<<setw(19) <<  var_fitness  <<setprecision(0)<<setw(19) <<  start_x_ <<setprecision(0)<<setw(19) <<  start_y_ <<setprecision(5)<<setw(19) <<  avg_x <<setprecision(5)<<setw(19) <<  avg_y <<setprecision(7)<<setw(19) <<  maxfit <<setprecision(0)<<setw(19) <<  IS_update_coords_.size() <<setprecision(0)<<setw(19) <<  tot_num_IS_upd <<setprecision(0)<<setw(19) <<  min_vir_x_ <<setprecision(0)<<setw(19) <<  min_vir_y_ <<setprecision(0)<<setw(19) <<  max_vir_x_ <<setprecision(0)<<setw(19) <<  max_vir_y_ <<setprecision(0)<<setw(19) <<  min_IS_x_ <<setprecision(0)<<setw(19) <<  min_IS_y_ <<setprecision(0)<<setw(19) <<  max_IS_x_ <<setprecision(0)<<setw(19) <<  max_IS_y_  <<setprecision(5)<<setw(19) <<  mean_displ_x_ <<setprecision(5)<<setw(19) <<  mean_displ_y_ <<setprecision(5)<<setw(19) <<  mean_displ_tot_<<setprecision(5)<<setw(19) <<  var_displ_x_ <<setprecision(5)<<setw(19) <<  var_displ_y_ <<setprecision(0)<<setw(19) <<  count_displ_ <<setprecision(0)<<setw(19) <<  x1_max_fitn <<setprecision(0)<<setw(19) <<  y1_max_fitn <<endl;
    
    
    
    displ_tot_    =0;
    displ_x_    =0;
    displ_y_    =0;
    displ_x_sq_ =0;
    displ_y_sq_ =0;
    count_displ_=0;
    

};


   
    
double antigenic_space::idx_to_real_x(int i1) const
{
    
    double k=(i1 + start_x_)*latt_width_;

    return k;   
    
};

    
    
    
double antigenic_space::idx_to_real_y(int j1) const
{
    
    double k=(j1 + start_y_)*latt_width_;

    return k;   
    
};


    
    
double antigenic_space::distance(int i1, int j1, int i2, int j2) const
{
    double x1= idx_to_real_x(i1);
    double y1= idx_to_real_y(j1);
  
    double x2= idx_to_real_x(i2);
    double y2= idx_to_real_y(j2);
    
    double dist =sqrt( pow( x1 - x2, 2.) + pow( y1 - y2, 2.));
    

    return dist;   
    
};


    
    

    
    

double kern( double distance, double rec, string ker) 
{
    
    double k=0.;
 
	if (ker=="heavyside"){
	    if (distance<=rec) k=1.; //implementation of the theta. Change if needed
	}
	else if (ker=="gaussian"){
	    
	    k=exp(- pow((distance/rec),2.)/2.);
	}
	else if (ker=="exponential"){
	    
	    k=exp(- (distance/rec));
	}
	
	
	else{
		cout << "wrong string for kernel, either heavyside, gaussian or exponential!" << endl;
		exit(1);
	    
	}
    

    return k;   
    
};



vector<double> antigenic_space::update_fitness_fullconv(bool upd_fitn){
    cout <<  " update_fitness_fullconv "   << endl;



    vector<double> new_fitn_full (vir_coords_.size(), 0.); 
  
    //~ cout <<  " full convolution "   << endl; 
     
    chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();
    
    for(int j=0; j < vir_coords_.size(); j++){ 
        
        double conv=0.;
    
        int x1 =  vir_coords_[j][0] - start_x_;
        int y1 =  vir_coords_[j][1] - start_y_;
           
        for(int i=0; i < IS_coords_.size(); i++){ 
         
            int x2 =  IS_coords_[i][0] - start_x_;
            int y2 =  IS_coords_[i][1] - start_y_;
            
            conv += (world_[x2][y2].get_num_IS())*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);
            
        }
        
        conv= conv/float(people_number_);
        
        double fitn = fitness_fct(conv, mem_points_, F0_);
        
       
        
        if (fitn<-1) fitn = -1;
        
        new_fitn_full[j] += fitn;
        
        //vir_coords_[j].set_fitness(fitn)
       if(upd_fitn)  world_[x1][y1].set_fitness(fitn);
        
    }
    
    chrono::steady_clock::time_point  end = std::chrono::steady_clock::now();  

  
    //~ cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
    
    //~ cout <<  " updated "   << endl;

    return new_fitn_full;
};


// convolve only the updated coordinates
vector<double> antigenic_space::update_fitness_fullconv_update(bool upd_fitn, bool combined_prep, bool combined_prep_newvirs, bool approx_fullIS, vector<vector<long int>> IS_update_coords_subset, vector<vector<int>> IS_coords_subset, double conv_precomp_update, double conv_precomp_tot){
    //cout << endl;
    
    cout <<  " update_fitness_fullconv UPDATE "   << endl;
    
    if (!combined_prep){
        IS_update_coords_subset=IS_update_coords_; 
        IS_coords_subset=IS_coords_;
        conv_precomp_update=0.;
        conv_precomp_tot=0.;
    }
    else if(!combined_prep_newvirs){
        IS_coords_subset=IS_coords_;
        conv_precomp_tot=0.;
    }
    
    //~ cout << IS_update_coords_subset.size() << " " <<IS_coords_subset.size()<<endl;
    //~ cout << conv_precomp_update << " " <<conv_precomp_tot<<endl;
    
    if (conv_precomp_update==0 || IS_update_coords_subset.size()==0){
        IS_update_coords_subset=IS_update_coords_;
        conv_precomp_update=0.;
        //if (approx_fullIS) IS_update_coords_subset= IS_coords_; // in this case _update label is wrong
        if (approx_fullIS){ // ugly but safer than reinterpret_cast
            vector<vector<long int>> veclong (IS_coords_.size(), vector<long int>(2));
            for (int i=0; i< IS_coords_.size(); i++) veclong[i]=vector<long int>((IS_coords_[i]).begin(), (IS_coords_[i]).end());
            IS_update_coords_subset= veclong; // in this case _update label is wrong
        }
    }
    if (conv_precomp_tot==0 || IS_coords_subset.size()==0){
        IS_coords_subset=IS_coords_;
        conv_precomp_tot=0.;
        
    }
    if (isnan(conv_precomp_tot)|| isnan(conv_precomp_update)){
 		cout << "nan precomp!" <<     conv_precomp_update   <<   "  "   <<  conv_precomp_tot  << endl;
		exit(1);
        
    }
    
    
    int count_old_vir=0;
    int count_new_vir=0;

    vector<double> new_fitn_full (vir_coords_.size(), 0.); 
  
    //~ cout <<  " full convolution "   << endl; 
     
    chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();

    for(int j=0; j < vir_coords_.size(); j++){ 
        
        double conv=0.;
        double fitn=0.;
    
        int x1 =  vir_coords_[j][0] - start_x_;
        int y1 =  vir_coords_[j][1] - start_y_;
        
        double fitness_prec=  world_[x1][y1].get_fitness();
        
        if (fitness_prec !=0 || approx_fullIS){ //  
            count_old_vir+=1;
               
            for(int i=0; i < IS_update_coords_subset.size(); i++){ 
             
                long int x2 =  IS_update_coords_subset[i][0] - start_x_;
                long int y2 =  IS_update_coords_subset[i][1] - start_y_;
                        
                long int num_upd;
                if(approx_fullIS)  num_upd = world_[x2][y2].get_num_IS();
                else num_upd =  IS_update_coords_subset[i][2];
            
                
                conv += (num_upd)*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);
                
            }
            
            
            double conv_tot= conv + conv_precomp_update;
            
            conv_tot= conv_tot/float(people_number_);
            
            if(!approx_fullIS) conv_tot+= fitness_inv_fct(fitness_prec, mem_points_, F0_);
            
            fitn = fitness_fct(conv_tot, mem_points_, F0_);
            
            //cout <<  " fitn "  << fitn -  fitness_prec << endl; 
            //cout  << IS_update_coords_subset.size() << endl; 
            
        }
        
        else if(! approx_fullIS){
            count_new_vir+=1;
            
            for(int i=0; i < IS_coords_subset.size(); i++){ 
             
                int x2 =  IS_coords_subset[i][0] - start_x_;
                int y2 =  IS_coords_subset[i][1] - start_y_;
                
                conv += (world_[x2][y2].get_num_IS())*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);
                
            }
            
            conv= conv + conv_precomp_tot;
            conv= conv/float(people_number_);

            fitn = fitness_fct(conv, mem_points_, F0_);

        }
        
        
        if (fitn<-1) fitn = -1;
        
  
        new_fitn_full[j] += fitn;
        
        
        if(upd_fitn) world_[x1][y1].set_fitness(fitn);
    
    } 
    
       
    chrono::steady_clock::time_point  end = std::chrono::steady_clock::now();  

  
    //~ cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
    

    
    //~ // I could reimplement on world with cutoff at like 7 rec_width_
    //~ cout <<   " old "   <<  count_old_vir   << endl;
    //~ cout <<  " new fullconv "   <<    count_new_vir   << endl;
    //~ cout <<  " updated "   << endl;

    return new_fitn_full;


};


// COMPUTE convolution on updated coordinates after farfield pre-computation
vector<double> antigenic_space::update_fitness_fullconv_update_precomp_thetabins(bool upd_fitn, bool combined_prep, bool combined_prep_newvirs, bool approx_fullIS, vector<vector<long int>> IS_update_coords_subset, vector<vector<int>> IS_coords_subset, vector<vector<double>> conv_precomp_update, vector<vector<double>> conv_precomp_tot){
    //cout << endl;
    
    cout <<  " update_fitness_fullconv "   << endl;
    
    if (!combined_prep){
        cout << "error, thetabins must be precomputed! " << endl;
        exit(1); 
    }
    else if(!combined_prep_newvirs){
        IS_coords_subset=IS_coords_;
        conv_precomp_tot=vector<vector<double>>();
    }
    
    //~ cout << IS_update_coords_subset.size() << " " <<IS_coords_subset.size()<<endl;
    //~ cout << conv_precomp_update.size() << " " <<conv_precomp_tot.size()<<endl;
    
    if (conv_precomp_update.size()==0 || IS_update_coords_subset.size()==0){
        IS_update_coords_subset=IS_update_coords_;
        conv_precomp_update=vector<vector<double>>();
        if (approx_fullIS){ // ugly but safer than reinterpret_cast
            vector<vector<long int>> veclong (IS_coords_.size(), vector<long int>(2));
            for (int i=0; i< IS_coords_.size(); i++) veclong[i]=vector<long int>((IS_coords_[i]).begin(), (IS_coords_[i]).end());
            IS_update_coords_subset= veclong; // in this case _update label is wrong
        }
        //if (approx_fullIS) IS_update_coords_subset= IS_coords_; // in this case _update label is wrong
    }
    if (conv_precomp_tot.size()==0 || IS_coords_subset.size()==0){
        IS_coords_subset=IS_coords_;
        conv_precomp_tot=vector<vector<double>>();
        
    }
    
        
    const double center_rect_x_vir= (min_vir_x_ + (max_vir_x_-min_vir_x_)/2.)*latt_width_;
    const double center_rect_y_vir= (min_vir_y_ + (max_vir_y_-min_vir_y_)/2.)*latt_width_;
  
    const double center_rect_x_IS= (min_IS_x_ + (max_IS_x_-min_IS_x_)/2.)*latt_width_;
    const double center_rect_y_IS= (min_IS_y_ + (max_IS_y_-min_IS_y_)/2.)*latt_width_;
     
    const double x2_ref_center_IS_notnorm = center_rect_x_IS - center_rect_x_vir; //center rect IS direction, unit cector
    const double y2_ref_center_IS_notnorm = center_rect_y_IS - center_rect_y_vir;
    const double x2_ref_center_IS = x2_ref_center_IS_notnorm/(sqrt(pow(x2_ref_center_IS_notnorm,2.) + pow(y2_ref_center_IS_notnorm,2.)));
    const double y2_ref_center_IS = y2_ref_center_IS_notnorm/(sqrt(pow(x2_ref_center_IS_notnorm,2.) + pow(y2_ref_center_IS_notnorm,2.)));
    
    const double theta_IS_ref=atan2(y2_ref_center_IS,x2_ref_center_IS);   

    int count_old_vir=0;
    int count_new_vir=0;

    vector<double> new_fitn_full (vir_coords_.size(), 0.); 
  
    //~ cout <<  " full convolution "   << endl; 
     
    chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();

    for(int j=0; j < vir_coords_.size(); j++){ 
        
        double conv=0.;
        double fitn=0.;
    
        int x1 =  vir_coords_[j][0] - start_x_;
        int y1 =  vir_coords_[j][1] - start_y_;
        
        double fitness_prec=  world_[x1][y1].get_fitness();

        double x1_ref_center_rv = idx_to_real_x(x1) - center_rect_x_vir; //specific IS direction
        double y1_ref_center_rv = idx_to_real_y(y1) - center_rect_y_vir;
        
        double dist_vir_center_rv= distance_real(x1_ref_center_rv, y1_ref_center_rv, 0., 0.);
       
        //theta w.r.t. IS ref
        // project change basis to ref and find theta in that basis
        double x2_ref_center_basis_ISc = x1_ref_center_rv* x2_ref_center_IS + y1_ref_center_rv *y2_ref_center_IS ; 
        double y2_ref_center_basis_ISc = -x1_ref_center_rv* y2_ref_center_IS + y1_ref_center_rv *x2_ref_center_IS ;
        
        double theta;
        if(y2_ref_center_basis_ISc==0 && x2_ref_center_basis_ISc==0) theta=0; 
        
        else theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
        
        if (fitness_prec !=0 || approx_fullIS){ // virus was already there
            count_old_vir+=1;
        
            for(int i=0; i < conv_precomp_update.size(); i++){ 
                double theta_bin=conv_precomp_update[i][0];
                
                if (theta_bin!=100){
            
                    //double theta=atan2(y1_ref_center_rv,x1_ref_center_rv);
                    
                    double x2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*cos(theta_bin);//ptojection on bin average direction, in the rotated basis, but only needed for distance, basis choice does not matter (orthonormal)
                    double y2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*sin(theta_bin);
                    
                    int sign_proj= 1;
                    if (cos(theta - theta_bin)<0) sign_proj= -1;
                    
                    double dist_vir_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.); // parallel bin
                    
                    if (isnan(conv_precomp_update[i][2])){
                        cout << "nan precomp update!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                        exit(1);
                        
                    }
                    
                    //cout << "nan precomp update! " <<     conv_precomp_update[i][2] <<endl;
                  
                    conv += conv_precomp_update[i][2]*exp(sign_proj*dist_vir_center/recog_width_);
                                
                    if (isnan(conv)){
                        cout << "nan conv precomp update!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                        exit(1);
                        
                    }
                }
    
            }
  
            for(int i=0; i < IS_update_coords_subset.size(); i++){ 
             
                long int x2 =  IS_update_coords_subset[i][0] - start_x_;
                long int y2 =  IS_update_coords_subset[i][1] - start_y_;
                        
                long int num_upd;
                if(approx_fullIS)  num_upd = world_[x2][y2].get_num_IS();
                else num_upd =  IS_update_coords_subset[i][2];
            
                
                conv += (num_upd)*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);// /float(people_number_)
              
                if (isnan(conv)){
                    cout << "nan conv update!" <<     i   <<   "  "   <<  num_upd << endl;
                    exit(1);
                    
                }
                
            }
            
            double conv_tot= conv/float(people_number_); //conv_precomp_update +
            if (! approx_fullIS) conv_tot+= fitness_inv_fct(fitness_prec, mem_points_, F0_); //conv_precomp_update +
            
            fitn = fitness_fct(conv_tot, mem_points_, F0_);
  
        }
        else if (! approx_fullIS){
            count_new_vir+=1;
           
            if(combined_prep_newvirs){
                
                for(int i=0; i < conv_precomp_tot.size(); i++){ 
                    double theta_bin=conv_precomp_tot[i][0];
                    
                    if (theta_bin!=100){
    
                        double x2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*cos(theta_bin);//ptojection on bin average direction
                        double y2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*sin(theta_bin);
                        
                        int sign_proj= 1;
                        if (cos(theta - theta_bin)<0) sign_proj= -1;
                        
                        double dist_vir_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.); // parallel bin
              
                        if (isnan(conv_precomp_tot[i][2])){
                            cout << "nan precomp tot!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                            exit(1);
                            
                        }
                        
                        //cout << "nan precomp tot! " <<     conv_precomp_tot[i][2] <<endl;
                        
                        conv += conv_precomp_tot[i][2]*exp(sign_proj*dist_vir_center/recog_width_);
                        
                        if (isnan(conv)){
                            cout << "nan conv precomp tot!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                            exit(1);
                            
                        }
                    }
                }        
            }
            
            for(int i=0; i < IS_coords_subset.size(); i++){ 
             
                int x2 =  IS_coords_subset[i][0] - start_x_;
                int y2 =  IS_coords_subset[i][1] - start_y_;
                
                conv += (world_[x2][y2].get_num_IS())*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);
            }
            conv=conv/float(people_number_);
            fitn = fitness_fct(conv, mem_points_, F0_);
        }
        
        
        if (fitn<-1) fitn = -1;
        
        new_fitn_full[j] += fitn;

                
        if (isnan(fitn)){
            cout << "nan fitn! " <<     j   <<   "  "   <<  fitness_prec  << endl;
            exit(1);
            
        }

        if(upd_fitn) world_[x1][y1].set_fitness(fitn);
    } 
    
       
    chrono::steady_clock::time_point  end = std::chrono::steady_clock::now();  

  
    //~ cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
    

    
    //~ // I could reimplement on world with cutoff at like 7 rec_width_
    //~ cout <<   " old "   <<  count_old_vir   << endl;
    //~ cout <<  " new fullconv "   <<    count_new_vir   << endl;
    //~ cout <<  " updated "   << endl;
    
    
    return new_fitn_full;


};


// OUTDATED, DISREGARD
void antigenic_space::update_fitness_IS_update_prep_farfield(vector<vector<long int>> & IS_update_coords_notapprox, double & conv_precomp_update, vector<int> & limits_IS_notapprox_update, int & count_old_vir, int & count_new_vir, vector<int> & idx_old_vir, vector<int> & idx_new_vir, bool approx_fullIS){//, bool combine_FFT
    cout <<  " update_fitness_farfield "   << endl;
    
    
    //if(combine_FFT){

    count_old_vir=0;
    count_new_vir=0;
    
    //vector<int> idx_old_vir;
    //vector<int> idx_new_vir;
    
    vector<vector<long int>> main_IS_toapprox=IS_update_coords_;
    //if(approx_fullIS) main_IS_toapprox=IS_coords_;
    if (approx_fullIS){ // ugly but safer than reinterpret_cast
        vector<vector<long int>> veclong (IS_coords_.size(), vector<long int>(2));
        for (int i=0; i< IS_coords_.size(); i++) veclong[i]=vector<long int>((IS_coords_[i]).begin(), (IS_coords_[i]).end());
        main_IS_toapprox= veclong; // in this case _update label is wrong
    }

    for(int j=0; j < vir_coords_.size(); j++){ 
        
    
        int x1 =  vir_coords_[j][0] - start_x_;
        int y1 =  vir_coords_[j][1] - start_y_;
        
        double fitness_prec=  world_[x1][y1].get_fitness();
                     
        if (fitness_prec !=0 ){ // virus was already there
            count_old_vir+=1;
                   
            idx_old_vir.push_back(j);
        }
        
        else{
            count_new_vir+=1;
            idx_new_vir.push_back(j);      
        }
        
    } 
        
    //}
       
    double center_rect_x_vir= (min_vir_x_ + (max_vir_x_-min_vir_x_)/2.)*latt_width_;
    double center_rect_y_vir= (min_vir_y_ + (max_vir_y_-min_vir_y_)/2.)*latt_width_;

    double width_rect_vir= (max_vir_x_-min_vir_x_)*latt_width_;
    double heigth_rect_vir= (max_vir_y_-min_vir_y_)*latt_width_;
    
    double maxdiag= sqrt(pow(heigth_rect_vir,2.)+pow(width_rect_vir,2.))/2.;
    
    limits_IS_notapprox_update = {numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::min(), numeric_limits<int>::min()}; // min_x, min_y max_x, max_y
    
    //~ cout <<  " update farfield precomputation "   << endl;
    
    
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    
    //vector<vector<long int>> IS_update_coords_notapprox;
    
    int count_far_field=0;
    int count_close_field=0;
    conv_precomp_update=0.;

    for(int i=0; i < main_IS_toapprox.size(); i++){ 
        
        bool approx = false;
     
        long int x2 =  main_IS_toapprox[i][0] - start_x_;
        long int y2 =  main_IS_toapprox[i][1] - start_y_;
        long int num_upd;
        if(approx_fullIS)  num_upd = world_[x2][y2].get_num_IS();
        else num_upd =  main_IS_toapprox[i][2];
        
        double x2_ref_center = idx_to_real_x(x2) - center_rect_x_vir;
        double y2_ref_center = idx_to_real_y(y2) - center_rect_y_vir;
        
        double dist_IS_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.);
        
    
        
        vector<double> intersections_rectvir = intersections_rect(x2_ref_center, y2_ref_center, heigth_rect_vir, width_rect_vir);
        
       
        double x2_perp_ref_center = y2_ref_center;
        double y2_perp_ref_center = - x2_ref_center;
        
        vector<double> intersections_perp_rectvir= intersections_rect(x2_perp_ref_center, y2_perp_ref_center, heigth_rect_vir, width_rect_vir);
         
        
        double dist_inter_center= distance_real(intersections_rectvir[0], intersections_rectvir[1], 0., 0.);
        double dist_inter_perp_center= distance_real(intersections_perp_rectvir[0], intersections_perp_rectvir[1], 0., 0.);
        
        // must be outside rectangle! And then 3 conditions
        
        bool out_rect=(fabs(x2_ref_center) > width_rect_vir/2.) || (fabs(y2_ref_center) > heigth_rect_vir/2.);
        
        bool separable_sqrt=(10. * pow(dist_inter_perp_center, 2.) < pow(dist_IS_center - dist_inter_center, 2.));
        
        bool separable_exp=(100. *pow(dist_inter_perp_center, 2.) < recog_width_ * fabs(dist_IS_center - dist_inter_center)) ;
        
        bool vir_indep=(20. * dist_inter_center < recog_width_ );
        
        bool vir_indep_close=(20. * (dist_IS_center +  maxdiag) < recog_width_ );
        
        
 
        
        if ( out_rect && separable_sqrt  && separable_exp && vir_indep ){
            
            conv_precomp_update += (num_upd)*kern(dist_IS_center, recog_width_, kernel_);// /float(people_number_)
            approx = true;
            count_far_field++;
            
        }
        //else if(vir_indep_close){
        //    conv_precomp_update += float(num_upd)/float(people_number_);
        //    approx = true;
        //    count_close_field++;
        //
        //}
        
        if(! approx){
            IS_update_coords_notapprox.push_back(main_IS_toapprox[i]);
            
            if (main_IS_toapprox[i][0] < limits_IS_notapprox_update[0]) limits_IS_notapprox_update[0] = main_IS_toapprox[i][0];
            if (main_IS_toapprox[i][1] < limits_IS_notapprox_update[1]) limits_IS_notapprox_update[1] = main_IS_toapprox[i][1];
            if (main_IS_toapprox[i][0] > limits_IS_notapprox_update[2]) limits_IS_notapprox_update[2] = main_IS_toapprox[i][0];
            if (main_IS_toapprox[i][1] > limits_IS_notapprox_update[3]) limits_IS_notapprox_update[3] = main_IS_toapprox[i][1];
            
        }
    }    
    
    
};



// precompute farfield approximation for convolution. 
void antigenic_space::update_fitness_IS_update_prep_farfield_thetabins(bool combined_prep_newvirs, vector<vector<long int>> & IS_update_coords_notapprox, vector<vector<int>> & IS_coords_notapprox, vector<vector<double>> & conv_precomp_update, vector<vector<double>> & conv_precomp_tot, vector<int> & limits_IS_notapprox_update, vector<int> & limits_IS_notapprox_tot, int & count_old_vir, int & count_new_vir, vector<int> & idx_old_vir, vector<int> & idx_new_vir, bool approx_fullIS){//, bool combine_FFT


    cout <<  " update_fitness_farfield theta "   << endl;
    
    
    //~ string FFT_e_f ("farfield_params.dat");

    //~ string outfile= out_dir_+ FFT_e_f;

    //~ fstream outstream;
    
    //~ outstream.open(outfile, fstream::out | fstream::app); 
    
    
    //~ outstream << fixed  <<setprecision(0)<<setw(19) << combined_prep_newvirs  <<setprecision(0)<<setw(19) << approx_fullIS  ; //  <<setprecision(0)<<setw(19) << thetabins 

        
    
    
    //if(combine_FFT){

    count_old_vir=0;
    count_new_vir=0;
    
    //vector<int> idx_old_vir;
    //vector<int> idx_new_vir;
    
    vector<vector<long int>> main_IS_toapprox=IS_update_coords_;
    if (approx_fullIS){ // 
        vector<vector<long int>> veclong (IS_coords_.size(), vector<long int>(2));
        for (int i=0; i< IS_coords_.size(); i++) veclong[i]=vector<long int>((IS_coords_[i]).begin(), (IS_coords_[i]).end());
        main_IS_toapprox= veclong; // 
    }
    
           
     double fitn_avg        =0.;
     double fitn_sq        =0.;
           
     double conv_avg        =0.;
     double conv_sq        =0.;
           
     long int tot_vir=0;

    for(int j=0; j < vir_coords_.size(); j++){ 
        
    
        int x1 =  vir_coords_[j][0] - start_x_;
        int y1 =  vir_coords_[j][1] - start_y_;
        
        double fitness_prec=  world_[x1][y1].get_fitness();
        
        if (fitness_prec !=0 ){ // virus was already there
            count_old_vir+=1;
                   
            idx_old_vir.push_back(j);
    
            long int num_vir=world_[x1][y1].get_num_vir();
        
            tot_vir+= num_vir;

            double conv_prec     =  fitness_inv_fct(fitness_prec, mem_points_, F0_);    

            conv_avg       += num_vir*conv_prec;
            conv_sq       += num_vir*pow(conv_prec,2.);

            fitn_avg       += num_vir*fitness_prec;
            fitn_sq       += num_vir*pow(fitness_prec,2.);
         
        }
        
        else{
            count_new_vir+=1;
            idx_new_vir.push_back(j);      
        }
        
    } 
      
     fitn_avg       /=tot_vir;
     fitn_sq        /=tot_vir;
     
     double fitn_stdev=sqrt(fitn_sq - pow(fitn_avg,2.) );
     
     const double abserr_thr= 0.001*fitn_stdev; //
       
     conv_avg       /=tot_vir;
     conv_sq        /=tot_vir;
     
     double conv_stdev=sqrt(conv_sq - pow(conv_avg,2.) );
     
     const double abserr_thr_conv= 0.001*conv_stdev/(1. + exp(2*F0_/mem_points_)*pow(conv_stdev/mem_points_,2.)); // max absolute convolution error w.r.t. stddev, in order to keep fitness error lower than 10^-2 stdev(fitness)
     
    //}
    
    // compute the minimal rectangle containing all viruses  and the one with IS
       
    const double center_rect_x_vir= (min_vir_x_ + (max_vir_x_-min_vir_x_)/2.)*latt_width_;
    const double center_rect_y_vir= (min_vir_y_ + (max_vir_y_-min_vir_y_)/2.)*latt_width_;

    const double width_rect_vir= (max_vir_x_-min_vir_x_)*latt_width_;
    const double heigth_rect_vir= (max_vir_y_-min_vir_y_)*latt_width_;
    
    const double maxdiag= sqrt(pow(heigth_rect_vir,2.)+pow(width_rect_vir,2.))/2.;

    const double center_rect_x_IS= (min_IS_x_ + (max_IS_x_-min_IS_x_)/2.)*latt_width_;
    const double center_rect_y_IS= (min_IS_y_ + (max_IS_y_-min_IS_y_)/2.)*latt_width_;
     
    const double x2_ref_center_IS_notnorm = center_rect_x_IS - center_rect_x_vir; //center rect IS direction, unit cector
    const double y2_ref_center_IS_notnorm = center_rect_y_IS - center_rect_y_vir;
    const double x2_ref_center_IS = x2_ref_center_IS_notnorm/(sqrt(pow(x2_ref_center_IS_notnorm,2.) + pow(y2_ref_center_IS_notnorm,2.)));
    const double y2_ref_center_IS = y2_ref_center_IS_notnorm/(sqrt(pow(x2_ref_center_IS_notnorm,2.) + pow(y2_ref_center_IS_notnorm,2.)));
    
    const double theta_IS_ref=atan2(y2_ref_center_IS,x2_ref_center_IS);   
    
    
 
 
    const double x2_perp_ref_center_IS_notnorm =   y2_ref_center_IS_notnorm;
    const double y2_perp_ref_center_IS_notnorm = - x2_ref_center_IS_notnorm;
    vector<double> intersections_perp_rectvir_IS_notnorm= intersections_rect(x2_perp_ref_center_IS_notnorm, y2_perp_ref_center_IS_notnorm, heigth_rect_vir, width_rect_vir);
    
    int sign_x;
    if (intersections_perp_rectvir_IS_notnorm[0]==0) sign_x = -1;
    else sign_x=intersections_perp_rectvir_IS_notnorm[0]/fabs(intersections_perp_rectvir_IS_notnorm[0]);
    
    int sign_y;
    if (intersections_perp_rectvir_IS_notnorm[1]==0) sign_y=-1;
    else sign_y=intersections_perp_rectvir_IS_notnorm[1]/fabs(intersections_perp_rectvir_IS_notnorm[1]);
    
    double corn_closest_x_1= sign_x*width_rect_vir/2.;
    double corn_closest_y_1= sign_y*heigth_rect_vir/2.;
    double corn_closest_x_2= -corn_closest_x_1;
    double corn_closest_y_2= -corn_closest_y_1;
    
    //distance between a "typical IS coordinate and the virus closest to it
    
    double dist_IS_corn1= distance_real(x2_ref_center_IS_notnorm, y2_ref_center_IS_notnorm, corn_closest_x_1, corn_closest_y_1);
    double dist_IS_corn2= distance_real(x2_ref_center_IS_notnorm, y2_ref_center_IS_notnorm, corn_closest_x_2, corn_closest_y_2);
    
    int x1_max_vir_tmp =  int((corn_closest_x_1 + center_rect_x_vir)/latt_width_) - start_x_; //  
    int y1_max_vir_tmp =  int((corn_closest_y_1 + center_rect_y_vir)/latt_width_) - start_y_;
    
    if (dist_IS_corn1 > dist_IS_corn2){
        x1_max_vir_tmp =  int((corn_closest_x_2 + center_rect_x_vir)/latt_width_) - start_x_; //  
        y1_max_vir_tmp =  int((corn_closest_y_2 + center_rect_y_vir)/latt_width_) - start_y_;
    }
    
    const int x1_max_vir = x1_max_vir_tmp; //  
    const int y1_max_vir = y1_max_vir_tmp;

    double conv_max_vir_tmp=0.;
    double fitn_max_vir_tmp=0.;

    
    for(int i=0; i < IS_coords_.size(); i++){ 
     
        int x2 =  IS_coords_[i][0] - start_x_;
        int y2 =  IS_coords_[i][1] - start_y_;
        
        conv_max_vir_tmp += (world_[x2][y2].get_num_IS())*kern(distance(x1_max_vir, y1_max_vir, x2, y2), recog_width_, kernel_);
        
    }
    
    conv_max_vir_tmp= conv_max_vir_tmp/float(people_number_);

    fitn_max_vir_tmp = fitness_fct(conv_max_vir_tmp, mem_points_, F0_); // fitness for the closest virus
    
    if (fitn_max_vir_tmp<-1) fitn_max_vir_tmp = -1;

    const double conv_max_vir=conv_max_vir_tmp;
    const double fitn_max_vir=fitn_max_vir_tmp;

    const double x1_ref_center_rv_max_vir = idx_to_real_x(x1_max_vir) - center_rect_x_vir; //specific IS direction
    const double y1_ref_center_rv_max_vir = idx_to_real_y(y1_max_vir) - center_rect_y_vir;
    
    const double dist_vir_center_rv_max_vir= distance_real(x1_ref_center_rv_max_vir, y1_ref_center_rv_max_vir, 0., 0.);
   
    //theta w.r.t. IS ref
    // project change basis to ref and find theta in that basis
    const double x2_ref_center_basis_ISc_max_vir = x1_ref_center_rv_max_vir* x2_ref_center_IS + y1_ref_center_rv_max_vir *y2_ref_center_IS ; 
    const double y2_ref_center_basis_ISc_max_vir = -x1_ref_center_rv_max_vir* y2_ref_center_IS + y1_ref_center_rv_max_vir *x2_ref_center_IS ;
    
    double theta_max_vir_tmp;
    if(y2_ref_center_basis_ISc_max_vir==0 && x2_ref_center_basis_ISc_max_vir==0) theta_max_vir_tmp=0; 
    else theta_max_vir_tmp=atan2(y2_ref_center_basis_ISc_max_vir,x2_ref_center_basis_ISc_max_vir); // it's out of vir rect, therefore can't be 0,0
    const double theta_max_vir=theta_max_vir_tmp;
    
    // rectangle containing the IS coordinates that can't be pecomputed
    limits_IS_notapprox_update = {numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::min(), numeric_limits<int>::min()}; // min_x, min_y max_x, max_y
    limits_IS_notapprox_tot = {numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::min(), numeric_limits<int>::min()}; // min_x, min_y max_x, max_y

    //~ cout <<  " update farfield theta precomputation "   << endl;
    
    double sepexp_molt=50.;
    
    double sepsqrt_molt=20.;
    
    
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    double theta_min = numeric_limits<double>::max();        
    double theta_max = - numeric_limits<double>::max();  
    
    double log_conv_approx_min = numeric_limits<double>::max();
    double log_conv_approx_max = - numeric_limits<double>::max();
    
    double abserr_tot = 0;  
    int not_abserr_red = 0;  
    
    vector<int> idx_approx_update;
    vector<double> conv_approx_vec;
    
    int Nbins=100;      
    
    // check if IS coord can be approximated and compute the min and max angles between any IS and the viruses center. Store approximated contribution to convolution 


    for(int i=0; i < main_IS_toapprox.size(); i++){ 
        
     
        long int x2 =  main_IS_toapprox[i][0] - start_x_;
        long int y2 =  main_IS_toapprox[i][1] - start_y_;
        long int num_upd;
        if(approx_fullIS)  num_upd = world_[x2][y2].get_num_IS();
        else num_upd =  main_IS_toapprox[i][2];
        
        double x2_ref_center = idx_to_real_x(x2) - center_rect_x_vir;
        double y2_ref_center = idx_to_real_y(y2) - center_rect_y_vir;
        
        double dist_IS_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.);
        
        vector<double> intersections_rectvir = intersections_rect(x2_ref_center, y2_ref_center, heigth_rect_vir, width_rect_vir);
        
       
        double x2_perp_ref_center = y2_ref_center;
        double y2_perp_ref_center = - x2_ref_center;
        
        vector<double> intersections_perp_rectvir= intersections_rect(x2_perp_ref_center, y2_perp_ref_center, heigth_rect_vir, width_rect_vir);
        
        
        double dist_inter_center= distance_real(intersections_rectvir[0], intersections_rectvir[1], 0., 0.);
        double dist_inter_perp_center= distance_real(intersections_perp_rectvir[0], intersections_perp_rectvir[1], 0., 0.);
        
        // must be outside rectangle! And then 3 conditions
        
        bool out_rect=(fabs(x2_ref_center) > width_rect_vir/2.) || (fabs(y2_ref_center) > heigth_rect_vir/2.);
        
        bool separable_sqrt=(sepsqrt_molt * pow(dist_inter_perp_center, 2.) < pow(dist_IS_center - dist_inter_center, 2.));
        
        bool separable_exp=(sepexp_molt *pow(dist_inter_perp_center, 2.) < recog_width_ * fabs(dist_IS_center - dist_inter_center));
                
                
        //theta w.r.t. IS ref
            
        // project change basis to ref and find theta in that basis
        double x2_ref_center_basis_ISc = x2_ref_center* x2_ref_center_IS + y2_ref_center *y2_ref_center_IS ; 
        double y2_ref_center_basis_ISc = -x2_ref_center* y2_ref_center_IS + y2_ref_center *x2_ref_center_IS ;
        
        double theta; // can be 0,0
        if(y2_ref_center_basis_ISc==0 && x2_ref_center_basis_ISc==0) theta=0; 
        else theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); 
    

        double x2_ref_center_max_vir = dist_vir_center_rv_max_vir*cos(theta_max_vir - theta)*cos(theta);//ptojection on bin average direction, in the rotated basis, but only needed for distance, basis choice does not matter (orthonormal)
        double y2_ref_center_max_vir = dist_vir_center_rv_max_vir*cos(theta_max_vir - theta)*sin(theta);
        
        int sign_proj= 1;
        if (cos(theta_max_vir - theta)<0) sign_proj= -1;
        
        double dist_vir_center= distance_real(x2_ref_center_max_vir, y2_ref_center_max_vir, 0., 0.); // parallel bin
    
        double conv_approx= ((num_upd)*kern(distance(x1_max_vir, y1_max_vir, x2, y2), recog_width_, kernel_) - (num_upd)*kern(dist_IS_center, recog_width_, kernel_)*exp(sign_proj*dist_vir_center/recog_width_))/float(people_number_);
    
        // convolution in topvir coord, removing approximation due to this IS coord alone, assuming small angle bins
        double conv_max_vir_approx= conv_max_vir - conv_approx;
         
        double abserr_upbound=fabs(fitness_fct(conv_max_vir_approx, mem_points_, F0_) - fitn_max_vir);
        
        bool abserr_red=(abserr_upbound < abserr_thr/1000.);
        
        
        if ( out_rect && separable_sqrt  && separable_exp && abserr_red){
            
            abserr_tot+=abserr_upbound; // fitness non linear, I can add errors like this only if very small
          
            double log_conv_approx;
            if (conv_approx==0 && conv_approx_vec.size()>0) log_conv_approx=log_conv_approx_min;
            else if (conv_approx==0 && conv_approx_vec.size()==0) log_conv_approx=log(fabs((num_upd/float(people_number_))*kern(dist_IS_center, recog_width_, kernel_)*exp(dist_inter_center/recog_width_)*pow(dist_inter_perp_center, 2.)/(2.*recog_width_ * fabs(dist_IS_center - dist_inter_center))));
            else log_conv_approx= log(fabs(conv_approx));
            
            idx_approx_update.push_back(i);
            conv_approx_vec.push_back(conv_approx);
            
            if (theta < theta_min) theta_min = theta;
            if (theta > theta_max) theta_max = theta;
             
            if (log_conv_approx < log_conv_approx_min) log_conv_approx_min = log_conv_approx;
            if (log_conv_approx > log_conv_approx_max) log_conv_approx_max = log_conv_approx;

        }
        else{ 
            IS_update_coords_notapprox.push_back(main_IS_toapprox[i]);
            
            
            if (main_IS_toapprox[i][0] < limits_IS_notapprox_update[0]) limits_IS_notapprox_update[0] = main_IS_toapprox[i][0];
            if (main_IS_toapprox[i][1] < limits_IS_notapprox_update[1]) limits_IS_notapprox_update[1] = main_IS_toapprox[i][1];
            if (main_IS_toapprox[i][0] > limits_IS_notapprox_update[2]) limits_IS_notapprox_update[2] = main_IS_toapprox[i][0];
            if (main_IS_toapprox[i][1] > limits_IS_notapprox_update[3]) limits_IS_notapprox_update[3] = main_IS_toapprox[i][1];
            
            if( out_rect && separable_sqrt  && separable_exp && !abserr_red) not_abserr_red++;
        }
                
    } 
    
    
    
    
    if(theta_min!=numeric_limits<double>::max() && theta_max != - numeric_limits<double>::max() && theta_min!=theta_max){
        
        theta_max =  theta_max + max((theta_max - theta_min)/(Nbins*10), pow(10., -9)); // theta max in last bin
        
        double theta_size=(theta_max - theta_min)/Nbins;
        
        if (theta_size < pow(10., -6)){
            theta_size = pow(10., -6);
            Nbins=ceil((theta_max - theta_min)/theta_size);
        } 
        
        
        // prepare a vector of angle bins where I put precomputed convolutions for each angle  
        
        int Nbin_conv=idx_approx_update.size()/1000;
        log_conv_approx_max =  log_conv_approx_max + (log_conv_approx_max - log_conv_approx_min)/(Nbin_conv*1000); // theta max in last bin

        double bin_logconv_size=(log_conv_approx_max - log_conv_approx_min)/Nbin_conv;
        vector<int> log_conv_count= vector<int>(Nbin_conv);
        
        vector<int> empty;
        vector<vector<int>> log_conv_bin_idxs (Nbin_conv, empty);
        
        // histogram convolution approximations 
        if (Nbin_conv>1 && abserr_tot > abserr_thr){
            
            
            for(int k=0; k < idx_approx_update.size(); k++){ 
                
                int i = idx_approx_update[k];
    
                double conv_approx= conv_approx_vec[k];
    
                double log_conv_approx;
                if (conv_approx==0) log_conv_approx=log_conv_approx_min;
                else log_conv_approx= log(fabs(conv_approx));
                                
                int key = int((log_conv_approx - log_conv_approx_min)/bin_logconv_size);
                log_conv_count[key] ++;
                
                log_conv_bin_idxs[key].push_back(k);
                
            }
            
            vector<int> idx_approx_update_prec=idx_approx_update;
            vector<int>().swap( idx_approx_update); 
            
            double approx_conv_cumul=conv_max_vir;
            
            for(int key=0; key < log_conv_count.size(); key++){ 
                if(log_conv_count[key]>0){
                    for(int j=0; j < log_conv_bin_idxs[key].size(); j++){ 
                        int k = log_conv_bin_idxs[key][j];
                        int i = idx_approx_update_prec[k];
                        
                        //now check cumulative fitness error 
                        double conv_approx= conv_approx_vec[k];
                        
                        double abserr_fitn_cumul=fabs(fitness_fct(approx_conv_cumul - conv_approx, mem_points_, F0_) - fitn_max_vir);
                
                        //if (abserr_fitn_cumul > abserr_thr) break;
                        if (abserr_fitn_cumul < abserr_thr) {
                            idx_approx_update.push_back(idx_approx_update_prec[k]);
                            approx_conv_cumul-= conv_approx;
                            
                        }
                        else{
                            IS_update_coords_notapprox.push_back(main_IS_toapprox[i]);
                            
                            if (main_IS_toapprox[i][0] < limits_IS_notapprox_update[0]) limits_IS_notapprox_update[0] = main_IS_toapprox[i][0];
                            if (main_IS_toapprox[i][1] < limits_IS_notapprox_update[1]) limits_IS_notapprox_update[1] = main_IS_toapprox[i][1];
                            if (main_IS_toapprox[i][0] > limits_IS_notapprox_update[2]) limits_IS_notapprox_update[2] = main_IS_toapprox[i][0];
                            if (main_IS_toapprox[i][1] > limits_IS_notapprox_update[3]) limits_IS_notapprox_update[3] = main_IS_toapprox[i][1];
                        
                        }
                    }
                }
            }
            
            
        
        int count_far_field=0;
        int count_close_field=0;
        conv_precomp_update= vector<vector<double>>(Nbins,vector<double>(3)); // theta_sum, count_IS_bin, conv_precomp.
        
        //~ cout <<  " pre second loop "   << endl;
        
        // compute prefactors
    
        for(int k=0; k < idx_approx_update.size(); k++){ 
            
            int i = idx_approx_update[k];
            
            long int x2 =  main_IS_toapprox[i][0] - start_x_;
            long int y2 =  main_IS_toapprox[i][1] - start_y_;
            long int num_upd;
            if(approx_fullIS)  num_upd = world_[x2][y2].get_num_IS();
            else num_upd =  main_IS_toapprox[i][2];
            
            double x2_ref_center = idx_to_real_x(x2) - center_rect_x_vir;
            double y2_ref_center = idx_to_real_y(y2) - center_rect_y_vir;
            
            double dist_IS_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.);
            
            
            double x2_ref_center_basis_ISc = x2_ref_center* x2_ref_center_IS + y2_ref_center *y2_ref_center_IS ; 
            double y2_ref_center_basis_ISc = -x2_ref_center* y2_ref_center_IS + y2_ref_center *x2_ref_center_IS ;
            
            double theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
                        
            //double theta_size=(theta_max - theta_min)/Nbins;
            int theta_bin=int((theta - theta_min)/theta_size);
            //cout <<  "SEC LOOP: theta_min "  <<  theta_min  <<" theta_max "  <<  theta_max  <<" theta_size "  <<  theta_size   <<" theta_IS_ref "  <<  theta_IS_ref   <<" theta_bin "  <<  theta_bin <<" theta "  <<  theta   << endl;
            
            conv_precomp_update[theta_bin][0]+=theta*abs(num_upd);
            
            //if (conv_precomp_update[theta_bin][1]!=0) cout << "non zero init " <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  conv_precomp_update[theta_bin][0] <<   "  "   <<  conv_precomp_update[theta_bin][1]  << endl;
            
            conv_precomp_update[theta_bin][1]+=abs(num_upd);
            
            if (isnan(conv_precomp_update[theta_bin][0]) || isnan(conv_precomp_update[theta_bin][1]) ){
                cout << "nan second loop precomp update! " <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  conv_precomp_update[theta_bin][0] <<   "  "   <<  conv_precomp_update[theta_bin][1]  << endl;
                exit(1);
            }
        }    
         //~ cout <<  " post second loop "   << endl;
       
    
        for(int i=0; i < conv_precomp_update.size(); i++){ 
            if( (long)(conv_precomp_update[i][1])<=0) conv_precomp_update[i][0]=100;
            else conv_precomp_update[i][0]/=conv_precomp_update[i][1];
            
            if (isnan(conv_precomp_update[i][0]) || isnan(conv_precomp_update[i][1]) ){
                cout << "nan avg theta loop precomp update! " <<     i   <<   "  "   <<  conv_precomp_update[i][0] <<   "  "   <<  conv_precomp_update[i][1]  << endl;
                exit(1);
                
            }

            //conv_precomp_update[i][1]=0;
        }
            
        //~ cout <<  " pre third loop "   << endl;
     
        for(int k=0; k < idx_approx_update.size(); k++){  
            //cout <<  " in third loop "   << endl;
            
            int i = idx_approx_update[k];
            
            bool approx = false;
            
         
            long int x2 =  main_IS_toapprox[i][0] - start_x_;
            long int y2 =  main_IS_toapprox[i][1] - start_y_;
            long int num_upd;
            if(approx_fullIS)  num_upd = world_[x2][y2].get_num_IS();
            else num_upd =  main_IS_toapprox[i][2];
                
            double x2_ref_center_rv = idx_to_real_x(x2) - center_rect_x_vir; //specific IS direction
            double y2_ref_center_rv = idx_to_real_y(y2) - center_rect_y_vir;
            
            double dist_IS_center_rv= distance_real(x2_ref_center_rv, y2_ref_center_rv, 0., 0.);
        
            double x2_ref_center_basis_ISc = x2_ref_center_rv* x2_ref_center_IS + y2_ref_center_rv *y2_ref_center_IS ; 
            double y2_ref_center_basis_ISc = -x2_ref_center_rv* y2_ref_center_IS + y2_ref_center_rv *x2_ref_center_IS ;
            
            double theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
            //double theta=atan2(y2_ref_center_rv,x2_ref_center_rv);
            
            if(theta>=theta_min && theta< theta_max){ //otherwise wasn't eligible
                
                //double theta_size=(theta_max - theta_min)/Nbins;
                int theta_bin=int((theta - theta_min)/theta_size);
                //cout <<  " i "  <<  i  <<" IS_update_coords_.size() "  <<  IS_update_coords_.size() << endl;
                //cout <<  " theta_min "  <<  theta_min  <<" theta_max "  <<  theta_max  <<" theta_size "  <<  theta_size   <<" theta_IS_ref "  <<  theta_IS_ref   <<" theta_bin "  <<  theta_bin <<" theta "  <<  theta     << endl;
                //cout  <<" conv_precomp_update[theta_bin][0] "  <<  conv_precomp_update[theta_bin][0]     << endl;
                
                if (cos(theta - conv_precomp_update[theta_bin][0])<0){
                    cout << "error, huge theta bin, projection wrong side IS upd! " << endl;
                    cout <<  " theta_min "  <<  theta_min  <<" theta_max "  <<  theta_max  <<" theta_size "  <<  theta_size   <<" theta_bin "  <<  theta_bin <<" theta "  <<  theta  <<" conv_precomp_update[theta_bin][0] "  <<  conv_precomp_update[theta_bin][0]    << endl;
                    exit(1); 
                }
                
                if (conv_precomp_update[theta_bin][0]==100){
                    cout << "empty theta bin! " << endl;
                    cout <<  " theta_min "  <<  theta_min  <<" theta_max "  <<  theta_max  <<" theta_size "  <<  theta_size   <<" theta_bin "  <<  theta_bin <<" theta "  <<  theta  <<" conv_precomp_update[theta_bin][0] "  <<  conv_precomp_update[theta_bin][0]   <<" conv_precomp_update[theta_bin][0] "  <<  conv_precomp_update.size()    << endl;
                    
                    cout <<  " i "  <<  i  <<" k "  <<  k  <<" num_upd "  <<  num_upd    <<" abs(num_upd) "  <<  abs(num_upd)   <<" x2_ref_center_basis_ISc "  <<  x2_ref_center_basis_ISc <<" y2_ref_center_basis_ISc "  <<  y2_ref_center_basis_ISc  <<"  dist_IS_center_rv "  <<  dist_IS_center_rv   <<"  conv_precomp_update[theta_bin][1] "  <<  conv_precomp_update[theta_bin][1]    <<" int conv_precomp_update[theta_bin][1] "  << int( conv_precomp_update[theta_bin][1])   <<" long conv_precomp_update[theta_bin][1] "  << (long)( conv_precomp_update[theta_bin][1])   << endl;
                    exit(1); 
                }
                
                double x2_ref_center = dist_IS_center_rv*cos(theta - conv_precomp_update[theta_bin][0])*cos(conv_precomp_update[theta_bin][0] + theta_IS_ref);//ptojection on bin average direction, in cardinal basis
                double y2_ref_center = dist_IS_center_rv*cos(theta - conv_precomp_update[theta_bin][0])*sin(conv_precomp_update[theta_bin][0] + theta_IS_ref);
                
                double dist_IS_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.);
                
                double dist_IS_center_perp= fabs(dist_IS_center_rv*sin(theta - conv_precomp_update[theta_bin][0]));
                
                vector<double> intersections_rectvir = intersections_rect(x2_ref_center, y2_ref_center, heigth_rect_vir, width_rect_vir);
                
               
                double x2_perp_ref_center = y2_ref_center;
                double y2_perp_ref_center = -x2_ref_center;
                
                vector<double> intersections_perp_rectvir= intersections_rect(x2_perp_ref_center, y2_perp_ref_center, heigth_rect_vir, width_rect_vir);
                 
                
                double dist_inter_center= distance_real(intersections_rectvir[0], intersections_rectvir[1], 0., 0.);
                double dist_inter_perp_center= distance_real(intersections_perp_rectvir[0], intersections_perp_rectvir[1], 0., 0.);
                
                // must be outside rectangle! And then 2 conditions
                
                bool out_rect=(fabs(x2_ref_center_rv) > width_rect_vir/2.) || (fabs(y2_ref_center_rv) > heigth_rect_vir/2.);
                
                bool separable_sqrt=(sepsqrt_molt * pow(dist_IS_center_perp + dist_inter_perp_center, 2.) < pow(dist_IS_center - dist_inter_center, 2.));
                
                bool separable_exp=(sepexp_molt *pow(dist_IS_center_perp + dist_inter_perp_center, 2.) < recog_width_ * fabs(dist_IS_center - dist_inter_center)) ;
                
                
                if ( out_rect && separable_sqrt  && separable_exp ){
                    
                    
                    //conv_precomp_update[theta_bin][1]+=abs(num_upd);
                    conv_precomp_update[theta_bin][2]+=(num_upd)*kern(dist_IS_center, recog_width_, kernel_);// /float(people_number_)
                           
                    if (isnan(conv_precomp_update[theta_bin][0]) || isnan(conv_precomp_update[theta_bin][1]) || isnan(conv_precomp_update[theta_bin][2]) ){
                        cout << "nan third loop fin precomp update! " <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  conv_precomp_update[theta_bin][0] <<   "  "   <<  conv_precomp_update[theta_bin][1]  <<   "  "   <<  conv_precomp_update[theta_bin][2]  << endl;
                        exit(1);
                        
                    }
                
                    approx = true;
                    count_far_field++;
                    
                }
            }
            if(! approx){
                IS_update_coords_notapprox.push_back(main_IS_toapprox[i]);
                
                
                
                if (main_IS_toapprox[i][0] < limits_IS_notapprox_update[0]) limits_IS_notapprox_update[0] = main_IS_toapprox[i][0];
                if (main_IS_toapprox[i][1] < limits_IS_notapprox_update[1]) limits_IS_notapprox_update[1] = main_IS_toapprox[i][1];
                if (main_IS_toapprox[i][0] > limits_IS_notapprox_update[2]) limits_IS_notapprox_update[2] = main_IS_toapprox[i][0];
                if (main_IS_toapprox[i][1] > limits_IS_notapprox_update[3]) limits_IS_notapprox_update[3] = main_IS_toapprox[i][1];
            
            }
            //cout <<  " out third loop "   << endl;
        }
       
        //~ cout <<   " FAR FIELD APPROX update:  main_IS_toapprox.size() "   <<  main_IS_toapprox.size()  << " IS_update_coords_notapprox.size() "   <<  IS_update_coords_notapprox.size()   << " count_far_field "   <<  count_far_field   << " count_close_field "   <<  count_close_field   << endl;
        //~ cout <<  " theta_min "  <<  theta_min  << " theta_max "  <<  theta_max  <<" theta_size "  <<  theta_size   <<" theta_IS_ref "  <<  theta_IS_ref   << endl;
        
        //~ outstream << fixed  << setprecision(0)<<setw(19) <<  main_IS_toapprox.size() << setprecision(0)<<setw(19) << IS_update_coords_notapprox.size() ;
        //~ outstream << fixed  << setprecision(10)<<setw(19) <<  theta_min << setprecision(10)<<setw(19) << theta_max  << setprecision(10)<<setw(19) << theta_size  << setprecision(10)<<setw(19) << theta_IS_ref << setprecision(0)<<setw(19) << Nbins ;


    }
    else {
        IS_update_coords_notapprox=vector<vector<long int>>();    
        limits_IS_notapprox_update = {numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::min(), numeric_limits<int>::min()}; // min_x, min_y max_x, max_y
        
        //~ outstream << fixed  << setprecision(0)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 ;
        //~ outstream << fixed  << setprecision(0)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 ;

        //~ outstream << fixed  << setprecision(0)<<setw(19) <<  -1 << setprecision(0)<<setw(19) << -1 ;
        //~ outstream << fixed  << setprecision(10)<<setw(19) <<  -1 << setprecision(10)<<setw(19) << -1  << setprecision(10)<<setw(19) << -1  << setprecision(10)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 ;

    }
    
 
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    
    //~ cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;   
    
    
    // REPEAT SAME PASSAGES FOR ALL POPULATED COORDINATES IF WE WANT TO DO THE FULL CONVOLUTION AFTER
    
    if (count_new_vir>0 && combined_prep_newvirs && !approx_fullIS){
        
        //~ cout <<  " IS tot farfield precomputation "   << endl;
        
        sepexp_molt =50.;
        sepsqrt_molt=20.;
        
        
        begin = std::chrono::steady_clock::now();
            
    
        double theta_min_tot = numeric_limits<double>::max();        
        double theta_max_tot = - numeric_limits<double>::max();  
        
        vector<int> idx_approx_tot;
    
        log_conv_approx_min = numeric_limits<double>::max();
        log_conv_approx_max = - numeric_limits<double>::max();
        
        abserr_tot = 0;  
        not_abserr_red = 0;  
        
        vector<double>().swap( conv_approx_vec); 
        
    
        for(int i=0; i < IS_coords_.size(); i++){ 
            
         
            long int x2 =  IS_coords_[i][0] - start_x_;
            long int y2 =  IS_coords_[i][1] - start_y_;
            long int num_upd =  world_[x2][y2].get_num_IS();
            
            double x2_ref_center = idx_to_real_x(x2) - center_rect_x_vir;
            double y2_ref_center = idx_to_real_y(y2) - center_rect_y_vir;
            
            double dist_IS_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.);
            
        
            
            vector<double> intersections_rectvir = intersections_rect(x2_ref_center, y2_ref_center, heigth_rect_vir, width_rect_vir);
            
           
            double x2_perp_ref_center = y2_ref_center;
            double y2_perp_ref_center = - x2_ref_center;
            
            vector<double> intersections_perp_rectvir= intersections_rect(x2_perp_ref_center, y2_perp_ref_center, heigth_rect_vir, width_rect_vir);
            
            
            double dist_inter_center= distance_real(intersections_rectvir[0], intersections_rectvir[1], 0., 0.);
            double dist_inter_perp_center= distance_real(intersections_perp_rectvir[0], intersections_perp_rectvir[1], 0., 0.);
            
            // must be outside rectangle! And then 3 conditions
            
            bool out_rect=(fabs(x2_ref_center) > width_rect_vir/2.) || (fabs(y2_ref_center) > heigth_rect_vir/2.);
            
            bool separable_sqrt=(sepsqrt_molt * pow(dist_inter_perp_center, 2.) < pow(dist_IS_center - dist_inter_center, 2.));
            
            bool separable_exp=(sepexp_molt *pow(dist_inter_perp_center, 2.) < recog_width_ * fabs(dist_IS_center - dist_inter_center));
                 
                
            
            //theta w.r.t. IS ref
                
            // project change basis to ref and find theta in that basis
            double x2_ref_center_basis_ISc = x2_ref_center* x2_ref_center_IS + y2_ref_center *y2_ref_center_IS ; 
            double y2_ref_center_basis_ISc = -x2_ref_center* y2_ref_center_IS + y2_ref_center *x2_ref_center_IS ;
            
            double theta; // can be 0,0
            if(y2_ref_center_basis_ISc==0 && x2_ref_center_basis_ISc==0) theta=0; 
            else theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); 
        
    
            double x2_ref_center_max_vir = dist_vir_center_rv_max_vir*cos(theta_max_vir - theta)*cos(theta);//ptojection on bin average direction, in the rotated basis, but only needed for distance, basis choice does not matter (orthonormal)
            double y2_ref_center_max_vir = dist_vir_center_rv_max_vir*cos(theta_max_vir - theta)*sin(theta);
            
            int sign_proj= 1;
            if (cos(theta_max_vir - theta)<0) sign_proj= -1;
            
            double dist_vir_center= distance_real(x2_ref_center_max_vir, y2_ref_center_max_vir, 0., 0.); // parallel bin
        
            double conv_approx= ((num_upd)*kern(distance(x1_max_vir, y1_max_vir, x2, y2), recog_width_, kernel_) - (num_upd)*kern(dist_IS_center, recog_width_, kernel_)*exp(sign_proj*dist_vir_center/recog_width_))/float(people_number_);
        
            // convolution in topvir coord, removing approximation due to this IS coord alone, assuming small angle bins
            double conv_max_vir_approx= conv_max_vir - conv_approx;
             
            double abserr_upbound=fabs(fitness_fct(conv_max_vir_approx, mem_points_, F0_) - fitn_max_vir);
            
            bool abserr_red=(abserr_upbound < abserr_thr/1000.);
            

            if ( out_rect && separable_sqrt  && separable_exp  && abserr_red){
                    
                idx_approx_tot.push_back(i);
                
                if (theta < theta_min_tot) theta_min_tot = theta;
                if (theta > theta_max_tot) theta_max_tot = theta;
                               
                abserr_tot+=abserr_upbound; // fitness non linear, I can add errors like this only if very small
              
                double log_conv_approx;
                if (conv_approx==0 && conv_approx_vec.size()>0) log_conv_approx=log_conv_approx_min;
                else if (conv_approx==0 && conv_approx_vec.size()==0) log_conv_approx=log(fabs((num_upd/float(people_number_))*kern(dist_IS_center, recog_width_, kernel_)*exp(dist_inter_center/recog_width_)*pow(dist_inter_perp_center, 2.)/(2.*recog_width_ * fabs(dist_IS_center - dist_inter_center))));
                else log_conv_approx= log(fabs(conv_approx));
                
                conv_approx_vec.push_back(conv_approx);
                
                if (log_conv_approx < log_conv_approx_min) log_conv_approx_min = log_conv_approx;
                if (log_conv_approx > log_conv_approx_max) log_conv_approx_max = log_conv_approx;
    
            }        
            else{ 
                IS_coords_notapprox.push_back(IS_coords_[i]);
                
                if( out_rect && separable_sqrt  && separable_exp && !abserr_red) not_abserr_red++;
                         
                if (IS_coords_[i][0] < limits_IS_notapprox_tot[0]) limits_IS_notapprox_tot[0] = IS_coords_[i][0];
                if (IS_coords_[i][1] < limits_IS_notapprox_tot[1]) limits_IS_notapprox_tot[1] = IS_coords_[i][1];
                if (IS_coords_[i][0] > limits_IS_notapprox_tot[2]) limits_IS_notapprox_tot[2] = IS_coords_[i][0];
                if (IS_coords_[i][1] > limits_IS_notapprox_tot[3]) limits_IS_notapprox_tot[3] = IS_coords_[i][1];
            }

        } 
        //theta_max_tot =  nextafter(theta_max_tot, numeric_limits<double>::max()); 
         
        //~ cout <<  " tot: abserr_thr "   <<  abserr_thr  << " not_abserr_red "   <<  not_abserr_red << " abserr_tot "   <<  abserr_tot << endl;
       //~ outstream << fixed  << setprecision(10)<<setw(19) << abserr_thr << setprecision(0)<<setw(19) << not_abserr_red <<setprecision(10)<<setw(19) << abserr_tot ;

        if(theta_min_tot!=numeric_limits<double>::max() && theta_max_tot != - numeric_limits<double>::max() && theta_min_tot!=theta_max_tot){
        
            
            double theta_size_prec=(theta_max - theta_min)/Nbins;
            
            int Nbins=100;
            if(theta_min!=numeric_limits<double>::max() && theta_max != - numeric_limits<double>::max() && theta_min!=theta_max)
                Nbins=ceil((theta_max_tot - theta_min_tot)/theta_size_prec);
            
            if(Nbins>1000)
                Nbins=1000;
            
            theta_max_tot =  theta_max_tot + max((theta_max_tot - theta_min_tot)/(Nbins*10), pow(10., -9)); 
            
            double theta_size=(theta_max_tot - theta_min_tot)/Nbins;
            
            
            if (theta_size < pow(10., -6)){
                theta_size = pow(10., -6);
                Nbins=ceil((theta_max_tot - theta_min_tot)/theta_size);
            } 
        
            //int Nbins=100;      
        

           
            //~ cout <<  "tot post first loop "   << endl;
            //~ cout <<  " IS_coords_notapprox.size() "   <<  IS_coords_notapprox.size()  << " idx_approx_tot.size() "   <<  idx_approx_tot.size()  << endl;
            //~ outstream << fixed  << setprecision(0)<<setw(19) << IS_coords_notapprox.size() << setprecision(0)<<setw(19) << idx_approx_tot.size() ;

            //cout <<  " theta_min_tot "  <<  theta_min_tot  << " theta_max_tot "  <<  theta_max_tot  <<" theta_size "  <<  theta_size   <<" theta_IS_ref "  <<  theta_IS_ref   << endl;
                
                    
            //~ cout <<  " tot pre bin convolution sort "   << endl; // countbin sort on convolution
            
            int Nbin_conv=idx_approx_tot.size()/1000;
            log_conv_approx_max =  log_conv_approx_max + (log_conv_approx_max - log_conv_approx_min)/(Nbin_conv*1000); // theta max in last bin
    
            double bin_logconv_size=(log_conv_approx_max - log_conv_approx_min)/Nbin_conv;
            vector<int> log_conv_count= vector<int>(Nbin_conv);
            
            
            vector<int> empty;
            vector<vector<int>> log_conv_bin_idxs (Nbin_conv, empty);
                
            if (Nbin_conv>1 && abserr_tot > abserr_thr){
                
                
                for(int k=0; k < idx_approx_tot.size(); k++){ 
                    
                    int i = idx_approx_tot[k];
        
                    double conv_approx= conv_approx_vec[k];
        
                    double log_conv_approx;
                    if (conv_approx==0) log_conv_approx=log_conv_approx_min;
                    else log_conv_approx= log(fabs(conv_approx));
                                    
                    int key = int((log_conv_approx - log_conv_approx_min)/bin_logconv_size);
                    log_conv_count[key] ++;
                    
                    log_conv_bin_idxs[key].push_back(k);
                    
                }
                
                vector<int> idx_approx_tot_prec=idx_approx_tot;
                vector<int>().swap( idx_approx_tot); 
                
                double approx_conv_cumul=conv_max_vir;
                
                for(int key=0; key < log_conv_count.size(); key++){ 
                    if(log_conv_count[key]>0){
                        for(int j=0; j < log_conv_bin_idxs[key].size(); j++){ 
                            int k = log_conv_bin_idxs[key][j];
                            int i = idx_approx_tot_prec[k];
                            
                            //now check cumulative fitness error 
                            double conv_approx= conv_approx_vec[k];
                            
                            double abserr_fitn_cumul=fabs(fitness_fct(approx_conv_cumul - conv_approx, mem_points_, F0_) - fitn_max_vir);
                    
                            //if (abserr_fitn_cumul > abserr_thr) break;
                            if (abserr_fitn_cumul < abserr_thr) {
                                idx_approx_tot.push_back(idx_approx_tot_prec[k]);
                                approx_conv_cumul-= conv_approx;
                                
                            }
                            else{
                                IS_coords_notapprox.push_back(IS_coords_[i]);
                                
                                if (IS_coords_[i][0] < limits_IS_notapprox_tot[0]) limits_IS_notapprox_tot[0] = IS_coords_[i][0];
                                if (IS_coords_[i][1] < limits_IS_notapprox_tot[1]) limits_IS_notapprox_tot[1] = IS_coords_[i][1];
                                if (IS_coords_[i][0] > limits_IS_notapprox_tot[2]) limits_IS_notapprox_tot[2] = IS_coords_[i][0];
                                if (IS_coords_[i][1] > limits_IS_notapprox_tot[3]) limits_IS_notapprox_tot[3] = IS_coords_[i][1];
                            }
                        }
                    }
                }
                
                
            }// else does not matter, already threw away big contributions.
    
            //~ cout <<  " IS_coords_notapprox.size() "   <<  IS_coords_notapprox.size()  << " idx_approx_tot.size() "   <<  idx_approx_tot.size()  << endl;
    
            //~ outstream << fixed  << setprecision(0)<<setw(19) << IS_coords_notapprox.size() << setprecision(0)<<setw(19) << idx_approx_tot.size() ;

            
            //vector<vector<long int>> IS_coords_notapprox;
            
            int count_far_field=0;
            int count_close_field=0;
            conv_precomp_tot= vector<vector<double>>(Nbins,vector<double>(3)); // theta_sum, count_IS_bin, conv_precomp.
            
            //~ cout <<  " tot pre second loop "   << endl;
      
            for(int k=0; k < idx_approx_tot.size(); k++){  
                //cout <<  " in third loop "   << endl;
                
                int i = idx_approx_tot[k];                
             
                long int x2 =  IS_coords_[i][0] - start_x_;
                long int y2 =  IS_coords_[i][1] - start_y_;
                long int num_upd =  world_[x2][y2].get_num_IS();
                
                double x2_ref_center = idx_to_real_x(x2) - center_rect_x_vir;
                double y2_ref_center = idx_to_real_y(y2) - center_rect_y_vir;
                
                double dist_IS_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.);
            
                double x2_ref_center_basis_ISc = x2_ref_center* x2_ref_center_IS + y2_ref_center *y2_ref_center_IS ; 
                double y2_ref_center_basis_ISc = -x2_ref_center* y2_ref_center_IS + y2_ref_center *x2_ref_center_IS ;
                
                double theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
                         
                //double theta_size=(theta_max_tot - theta_min_tot)/Nbins;
                int theta_bin=int((theta - theta_min_tot)/theta_size);
                                   
                if (isnan(conv_precomp_tot[theta_bin][0]) || isnan(conv_precomp_tot[theta_bin][1]) ){
                    cout << "nan second loop precomp tot! " <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  conv_precomp_tot[theta_bin][0] <<   "  "   <<  conv_precomp_tot[theta_bin][1]  << endl;
                    exit(1);
                    
                }
            
                conv_precomp_tot[theta_bin][0]+=theta*abs(num_upd);
                conv_precomp_tot[theta_bin][1]+=abs(num_upd);
            
            }    
             //~ cout <<  "tot post second loop "   << endl;

        
            for(int i=0; i < conv_precomp_tot.size(); i++){ 
                if((long)(conv_precomp_tot[i][1])<=0) conv_precomp_tot[i][0]=100;
                else conv_precomp_tot[i][0]/=conv_precomp_tot[i][1];
 
                if (isnan(conv_precomp_tot[i][0]) || isnan(conv_precomp_tot[i][1]) ){
                    cout << "nan avg theta loop precomp update! " <<     i   <<   "  "   <<  conv_precomp_tot[i][0] <<   "  "   <<  conv_precomp_tot[i][1]  << endl;
                    exit(1);
                    
                }
    
                conv_precomp_tot[i][1]=0;
            }
                
             //~ cout <<  "tot pre third loop "   << endl;

            for(int k=0; k < idx_approx_tot.size(); k++){  
                //cout <<  " in third loop "   << endl;
                
                int i = idx_approx_tot[k];                
                
                bool approx = false;
                
             
                long int x2 =  IS_coords_[i][0] - start_x_;
                long int y2 =  IS_coords_[i][1] - start_y_;
                long int num_upd =  world_[x2][y2].get_num_IS();
                
                double x2_ref_center_rv = idx_to_real_x(x2) - center_rect_x_vir; //specific IS direction
                double y2_ref_center_rv = idx_to_real_y(y2) - center_rect_y_vir;
                
                double dist_IS_center_rv= distance_real(x2_ref_center_rv, y2_ref_center_rv, 0., 0.);
                    
                double x2_ref_center_basis_ISc = x2_ref_center_rv* x2_ref_center_IS + y2_ref_center_rv *y2_ref_center_IS ; 
                double y2_ref_center_basis_ISc = -x2_ref_center_rv* y2_ref_center_IS + y2_ref_center_rv *x2_ref_center_IS ;
                
                double theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
                //double theta=atan2(y2_ref_center_rv,x2_ref_center_rv);
                
                if(theta>=theta_min_tot && theta< theta_max_tot){ //otherwise wasn't eligible
    
                        
                    //double theta_size=(theta_max_tot - theta_min_tot)/Nbins;
                    int theta_bin=int((theta - theta_min_tot)/theta_size);
                               
                    if (cos(theta - conv_precomp_tot[theta_bin][0])<0){
                        cout << "error, huge theta bin, projection wrong side IS tot! " << endl;
                        exit(1); 
                    }
                                    
                    if (conv_precomp_tot[theta_bin][0]==100){
                        cout << "empty theta bin IS tot! " << endl;
                        cout <<  " theta_min "  <<  theta_min  <<" theta_max "  <<  theta_max  <<" theta_size "  <<  theta_size   <<" theta_bin "  <<  theta_bin <<" theta "  <<  theta  <<" conv_precomp_tot[theta_bin][0] "  <<  conv_precomp_tot[theta_bin][0]    << endl;
                        exit(1); 
                    }
                
         
                    double x2_ref_center = dist_IS_center_rv*cos(theta - conv_precomp_tot[theta_bin][0])*cos(conv_precomp_tot[theta_bin][0] + theta_IS_ref);//ptojection on bin average direction
                    double y2_ref_center = dist_IS_center_rv*cos(theta - conv_precomp_tot[theta_bin][0])*sin(conv_precomp_tot[theta_bin][0] + theta_IS_ref);
                    
                    double dist_IS_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.);
                    
                    double dist_IS_center_perp= fabs(dist_IS_center_rv*sin(theta - conv_precomp_tot[theta_bin][0]));
                    
                    vector<double> intersections_rectvir = intersections_rect(x2_ref_center, y2_ref_center, heigth_rect_vir, width_rect_vir);
                    
                   
                    double x2_perp_ref_center = y2_ref_center;
                    double y2_perp_ref_center = - x2_ref_center;
                    
                    vector<double> intersections_perp_rectvir= intersections_rect(x2_perp_ref_center, y2_perp_ref_center, heigth_rect_vir, width_rect_vir);
                     
                    
                    double dist_inter_center= distance_real(intersections_rectvir[0], intersections_rectvir[1], 0., 0.);
                    double dist_inter_perp_center= distance_real(intersections_perp_rectvir[0], intersections_perp_rectvir[1], 0., 0.);
                    
                    // must be outside rectangle! And then 3 conditions
                    
                    bool out_rect=(fabs(x2_ref_center_rv) > width_rect_vir/2.) || (fabs(y2_ref_center_rv) > heigth_rect_vir/2.);
                    
                    bool separable_sqrt=(sepsqrt_molt * pow(dist_IS_center_perp + dist_inter_perp_center, 2.) < pow(dist_IS_center - dist_inter_center, 2.));
                    
                    bool separable_exp=(sepexp_molt *pow(dist_IS_center_perp + dist_inter_perp_center, 2.) < recog_width_ * fabs(dist_IS_center - dist_inter_center)) ;
                    
                    
                    if ( out_rect && separable_sqrt  && separable_exp ){
                        
                        
                        conv_precomp_tot[theta_bin][1]+=num_upd;
                        conv_precomp_tot[theta_bin][2]+=(num_upd)*kern(dist_IS_center, recog_width_, kernel_);// /float(people_number_)
                                          
                        if (isnan(conv_precomp_tot[theta_bin][0]) || isnan(conv_precomp_tot[theta_bin][1]) || isnan(conv_precomp_tot[theta_bin][2]) ){
                            cout << "nan third loop fin precomp tot! " <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  conv_precomp_tot[theta_bin][0] <<   "  "   <<  conv_precomp_tot[theta_bin][1]  <<   "  "   <<  conv_precomp_tot[theta_bin][2]  << endl;
                            exit(1);
                            
                        }

                        approx = true;
                        count_far_field++;
                        
                    }
                }
                if(! approx){
                    IS_coords_notapprox.push_back(IS_coords_[i]);
                    
                    if (IS_coords_[i][0] < limits_IS_notapprox_tot[0]) limits_IS_notapprox_tot[0] = IS_coords_[i][0];
                    if (IS_coords_[i][1] < limits_IS_notapprox_tot[1]) limits_IS_notapprox_tot[1] = IS_coords_[i][1];
                    if (IS_coords_[i][0] > limits_IS_notapprox_tot[2]) limits_IS_notapprox_tot[2] = IS_coords_[i][0];
                    if (IS_coords_[i][1] > limits_IS_notapprox_tot[3]) limits_IS_notapprox_tot[3] = IS_coords_[i][1];
                }
            } 
                
                
        }
        else {
            IS_coords_notapprox=vector<vector<int>>();    
            
            limits_IS_notapprox_tot = {numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::min(), numeric_limits<int>::min()}; // min_x, min_y max_x, max_y
           
           
        }     
         
        //~ end = std::chrono::steady_clock::now();
        
        //~ cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;   
        
    }
    else{
        //~ outstream << fixed  << setprecision(10)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 <<setprecision(10)<<setw(19) << -1 ;
    
        //~ outstream << fixed  << setprecision(0)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 ;
        //~ outstream << fixed  << setprecision(0)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 ;
    
        //~ outstream << fixed  << setprecision(0)<<setw(19) <<  -1 << setprecision(0)<<setw(19) << -1 ;
        //~ outstream << fixed  << setprecision(10)<<setw(19) <<  -1 << setprecision(10)<<setw(19) << -1  << setprecision(10)<<setw(19) << -1  << setprecision(10)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 ;

    }
    
    
    
    //~ outstream  <<endl;
    
};


// combines the farfield pre-computation and the exact convolution on updated coordinates
vector<double> antigenic_space::update_fitness_IS_update_combine_farfield_fullconv(bool upd_fitn, bool thetabins, bool combined_prep_newvirs, bool approx_fullIS){
    
    vector<vector<long int>> IS_update_coords_notapprox; 
    vector<vector<int>> IS_coords_notapprox;
    double conv_precomp_update=0.;
    double conv_precomp_tot=0.;
    vector<vector<double>> conv_precomp_update_theta;
    vector<vector<double>> conv_precomp_tot_theta;
    vector<int> limits_IS_notapprox_update {0,0,0,0};
    vector<int> limits_IS_notapprox_tot {0,0,0,0};
    int count_old_vir=-1;
    int count_new_vir=-1;
    vector<int> idx_old_vir;
    vector<int> idx_new_vir;
    
    vector<double> new_fitn;
    
    if(thetabins){
            
        update_fitness_IS_update_prep_farfield_thetabins(combined_prep_newvirs, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update_theta, conv_precomp_tot_theta, limits_IS_notapprox_update, limits_IS_notapprox_tot, count_old_vir, count_new_vir, idx_old_vir, idx_new_vir, approx_fullIS);
    
        
        new_fitn = update_fitness_fullconv_update_precomp_thetabins(upd_fitn, true, combined_prep_newvirs, approx_fullIS, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update_theta, conv_precomp_tot_theta);
    }
    else{
       
        update_fitness_IS_update_prep_farfield(IS_update_coords_notapprox, conv_precomp_update, limits_IS_notapprox_update, count_old_vir, count_new_vir, idx_old_vir, idx_new_vir, approx_fullIS);
        combined_prep_newvirs=false;
        
        new_fitn = update_fitness_fullconv_update(upd_fitn, true, combined_prep_newvirs, approx_fullIS, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update, conv_precomp_tot);
    }
        

    
    //~ string FFT_e_f ("farfield_errors.dat");

    //~ string outfile= out_dir_+ FFT_e_f;

    //~ fstream outstream;
    
    //~ outstream.open(outfile, fstream::out | fstream::app); 
     
    // check against a full conv done on the fitness nose and the coordinate with max number viruses, and save results to file
    
    vector<int> idxs_max= get_maxobs(); //num virs, IS, fitn
    
    vector<double> FFT_errors; //error fullconv to do vs actual. contains for each entry: fitn1, fitn2, abs(fitn1 - fitn2). Total 9 elements 
    
    for(int j=0; j < idxs_max.size(); j++){ 
        
        int idx = idxs_max[j];
        
        int x1 =  vir_coords_[idx][0] - start_x_;
        int y1 =  vir_coords_[idx][1] - start_y_;
        
        double fitness_FFT=  new_fitn[idx];
        
        FFT_errors.push_back(fitness_FFT);
        
    
        double conv=0.;
        double fitn=0.;
    
        
        for(int i=0; i < IS_coords_.size(); i++){ 
         
            int x2 =  IS_coords_[i][0] - start_x_;
            int y2 =  IS_coords_[i][1] - start_y_;
            
            conv += (world_[x2][y2].get_num_IS())*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);
            
        }
        
        conv= conv/float(people_number_);

        fitn = fitness_fct(conv, mem_points_, F0_);

        
        
        if (fitn<-1) fitn = -1;
        
         FFT_errors.push_back(fitn);
         
         FFT_errors.push_back(abs(fitness_FFT- fitn));
        
       
    
        //~ outstream << fixed  << setprecision(10)<<setw(19) << fitness_FFT << setprecision(10)<<setw(19) << fitn <<setprecision(10)<<setw(19) << abs(fitness_FFT- fitn) ;

        
    }
    
    
    //~ outstream  <<endl;
    
    return new_fitn;
};


// NFFT convolution

vector<double> antigenic_space::update_fitness_IS_update_NFFT(bool upd_fitn, bool combined_prep, bool combined_prep_newvirs, bool thetabins , bool FFT_newvirs_gen , bool approx_fullIS , vector<vector<long int>> IS_update_coords_subset , vector<vector<int>> IS_coords_subset, vector<vector<double>> conv_precomp_update, vector<vector<double>>  conv_precomp_tot, vector<int> limits_IS_notapprox_update, vector<int> limits_IS_notapprox_tot, int count_old_vir, int count_new_vir, vector<int> idx_old_vir, vector<int> idx_new_vir, bool swap_tot_upd_FFT){
    //, int min_IS_notapprox_x, int min_IS_notapprox_y, int max_IS_notapprox_x, int max_IS_notapprox_y

    cout <<  " update_fitness_fullconv update FFT  "   << endl;
    //~ cout <<  " initialize FFT plan "   << endl;
    
    
    if(swap_tot_upd_FFT){
        
        if(! thetabins || ! combined_prep || ! combined_prep_newvirs || approx_fullIS)
        {
            cout << "error, NFFT swap_tot_upd_FFT incompatible options! " << endl;
            exit(1); 
        }
        
        cout <<  " swap_tot_upd_FFT "   << endl;
        
        FFT_newvirs_gen=false;
    
    }

    
    
    

    
    if (approx_fullIS) FFT_newvirs_gen=false;
   
      
    if (!combined_prep){
        IS_coords_subset=IS_coords_;
        conv_precomp_tot=vector<vector<double>>();
        IS_update_coords_subset=IS_update_coords_;
        conv_precomp_update=vector<vector<double>>();
        
        if(combined_prep_newvirs) {    
            cout << "error, NFFT combined_prep_newvirs must be precomputed! " << endl;
            exit(1); 
        }
        if(thetabins) {    
            cout << "error, NFFT thetabins must be precomputed! " << endl;
            exit(1); 
        }
    }
    else if(idx_old_vir.size() + idx_new_vir.size()==0 || count_old_vir + count_new_vir ==0){
        cout << "error, NFFT precomputed wrong virus classification! " << endl;
        exit(1); 
    }
    else if(!combined_prep_newvirs){
        IS_coords_subset=IS_coords_;
        conv_precomp_tot=vector<vector<double>>();
    }
    
    //~ cout << IS_update_coords_subset.size() << " " <<IS_coords_subset.size()<<endl;
    //~ cout << conv_precomp_update.size() << " " <<conv_precomp_tot.size()<<endl;
    
    bool no_update_subset=false;
    
    if (conv_precomp_update.size()==0 || IS_update_coords_subset.size()==0 || limits_IS_notapprox_update.size()==0){ // if not prep, or did not find approximations
        IS_update_coords_subset=IS_update_coords_;
        conv_precomp_update=vector<vector<double>>();
        no_update_subset=true;
        if (approx_fullIS){ // ugly but safer than reinterpret_cast
            vector<vector<long int>> veclong (IS_coords_.size(), vector<long int>(2));
            for (int i=0; i< IS_coords_.size(); i++) veclong[i]=vector<long int>((IS_coords_[i]).begin(), (IS_coords_[i]).end());
            IS_update_coords_subset= veclong; // in this case _update label is wrong
            
            for (int i=0; i< IS_coords_.size(); i++) {
                for (int j=0; j< IS_coords_[0].size(); j++) {
                    if(IS_coords_[i][j]!=IS_update_coords_subset[i][j]){
                        cout << "error, fullIS approx no prepared, wrong IS_update_coords_subset i "  <<  i  << " j "   <<  j   << " IS_coords_[i][j] "   <<  IS_coords_[i][j]    << " IS_update_coords_subset[i][j] "   <<  IS_update_coords_subset[i][j]   << endl;
                        exit(1); 
                    }
                }
            }
            
        }
        
    }
    if (conv_precomp_tot.size()==0 || IS_coords_subset.size()==0 || limits_IS_notapprox_tot.size()==0){
        IS_coords_subset=IS_coords_;
        conv_precomp_tot=vector<vector<double>>();
        combined_prep_newvirs=false;
        
    }
    if (conv_precomp_tot.size()==0 && conv_precomp_update.size()==0){
        combined_prep=false;
        combined_prep_newvirs=false;
        thetabins=false;
    }
    
    if (combined_prep && !thetabins){
        if(conv_precomp_tot.size()!=1 || conv_precomp_update.size()!=1 || conv_precomp_tot[0].size()!=1 || conv_precomp_update[0].size()!=1 ){
            cout << "error, precomputed not theta conv preconmuted must be scalar! " << endl;
            exit(1); 
        }
    }
   
   

    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  
  
    if(!combined_prep){
    
        //int count_old_vir=0;
        //int count_new_vir=0;
        
        //vector<int> idx_old_vir;
        //vector<int> idx_new_vir;
        
        idx_old_vir=vector<int>();
        idx_new_vir=vector<int>();
        count_old_vir=0;
        count_new_vir=0;
    
        for(int j=0; j < vir_coords_.size(); j++){ 
            
        
            int x1 =  vir_coords_[j][0] - start_x_;
            int y1 =  vir_coords_[j][1] - start_y_;
            
            double fitness_prec=  world_[x1][y1].get_fitness();
                         
            if (fitness_prec !=0 ){ // virus was already there
                count_old_vir+=1;
                       
                idx_old_vir.push_back(j);
            }
            
            else{
                count_new_vir+=1;
                idx_new_vir.push_back(j);      
            }
            
        } 
    }
    if(!combined_prep || no_update_subset){   
         
        //min_IS_notapprox_x = numeric_limits<int>::max();        
        //min_IS_notapprox_y = numeric_limits<int>::max();        
        //max_IS_notapprox_x = numeric_limits<int>::min();        
        //max_IS_notapprox_y = numeric_limits<int>::min();  
        
        limits_IS_notapprox_update = {numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::min(), numeric_limits<int>::min()}; // min_x, min_y max_x, max_y
    
        for(int i=0; i < IS_update_coords_subset.size(); i++){ 
        
            if (IS_update_coords_subset[i][0] < limits_IS_notapprox_update[0]) limits_IS_notapprox_update[0] = IS_update_coords_subset[i][0];
            if (IS_update_coords_subset[i][1] < limits_IS_notapprox_update[1]) limits_IS_notapprox_update[1] = IS_update_coords_subset[i][1];
            if (IS_update_coords_subset[i][0] > limits_IS_notapprox_update[2]) limits_IS_notapprox_update[2] = IS_update_coords_subset[i][0];
            if (IS_update_coords_subset[i][1] > limits_IS_notapprox_update[3]) limits_IS_notapprox_update[3] = IS_update_coords_subset[i][1];
            
        }
        
    
    }
    
    if(!combined_prep || !combined_prep_newvirs){// && FFT_newvirs_gen
        //if (min_IS_x_ < min_IS_notapprox_x) min_IS_notapprox_x = min_IS_x_;
        //if (min_IS_y_ < min_IS_notapprox_y) min_IS_notapprox_y = min_IS_y_;
        //if (max_IS_x_ > max_IS_notapprox_x) max_IS_notapprox_x = max_IS_x_;
        //if (max_IS_y_ > max_IS_notapprox_y) max_IS_notapprox_y = max_IS_y_;
        limits_IS_notapprox_tot = {numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::min(), numeric_limits<int>::min()}; // min_x, min_y max_x, max_y

        limits_IS_notapprox_tot[0] = min_IS_x_;
        limits_IS_notapprox_tot[1] = min_IS_y_;
        limits_IS_notapprox_tot[2] = max_IS_x_;
        limits_IS_notapprox_tot[3] = max_IS_y_;
    }
    
    
    if(swap_tot_upd_FFT){
    
        vector<vector<long int>> veclong (IS_coords_subset.size(), vector<long int>(2));
        for (int i=0; i< IS_coords_subset.size(); i++) veclong[i]=vector<long int>((IS_coords_subset[i]).begin(), (IS_coords_subset[i]).end());
        
        vector<vector<int>> veci (IS_update_coords_subset.size(), vector<int>(3));
        for (int i=0; i< IS_update_coords_subset.size(); i++) veci[i]=vector<int>((IS_update_coords_subset[i]).begin(), (IS_update_coords_subset[i]).end());
        
        IS_update_coords_subset= veclong; // in this case _update label is wrong
        
        IS_coords_subset=veci;

        //swap(IS_update_coords_subset, IS_coords_subset);
        swap(conv_precomp_update, conv_precomp_tot);
        swap(limits_IS_notapprox_update, limits_IS_notapprox_tot);
        swap(idx_old_vir, idx_new_vir);
        swap(count_old_vir, count_new_vir);
    }

    

    // rescale space in a ball with radius 0.25
    const double min_x_real=min(min_vir_x_, limits_IS_notapprox_update[0])*latt_width_;
    const double min_y_real=min(min_vir_y_, limits_IS_notapprox_update[1])*latt_width_;
    const double max_x_real=max(max_vir_x_, limits_IS_notapprox_update[2])*latt_width_;
    const double max_y_real=max(max_vir_y_, limits_IS_notapprox_update[3])*latt_width_;
    
    const double center_rect_x= min_x_real + (max_x_real-min_x_real)/2.;
    const double center_rect_y= min_y_real + (max_y_real-min_y_real)/2.;
    double maxradius_tmp= sqrt(pow((max_y_real-min_y_real),2.)+pow((max_x_real-min_x_real),2.))/2. + 1;
 
    const double min_x_real_tot=min(min_vir_x_, limits_IS_notapprox_tot[0])*latt_width_;
    const double min_y_real_tot=min(min_vir_y_, limits_IS_notapprox_tot[1])*latt_width_;
    const double max_x_real_tot=max(max_vir_x_, limits_IS_notapprox_tot[2])*latt_width_;
    const double max_y_real_tot=max(max_vir_y_, limits_IS_notapprox_tot[3])*latt_width_;
    
    const double center_rect_x_tot= min_x_real_tot + (max_x_real_tot-min_x_real_tot)/2.;
    const double center_rect_y_tot= min_y_real_tot + (max_y_real_tot-min_y_real_tot)/2.;
    double maxradius_tot_tmp= sqrt(pow((max_y_real_tot-min_y_real_tot),2.)+pow((max_x_real_tot-min_x_real_tot),2.))/2. + 1;
 
    
    if (approx_fullIS){
        count_old_vir=vir_coords_.size();
        count_new_vir=0;
        
        idx_new_vir=vector<int>();
        idx_old_vir=vector<int>(vir_coords_.size());
        iota(idx_old_vir.begin(), idx_old_vir.end(), 0); // Fill with 0, 1, ..., vir_coords_.size() -1. In place. c++11
    }
    
    bool FFT_newvir=false;
    
    //if (count_new_vir > 30 && FFT_newvirs_gen) FFT_newvir=true; // FFT with too unbalanced and non uniform M and N works worse. Depends on c. Possibly better combining with far field  && count_new_vir > IS_coords_subset.size()/10000
    if (count_new_vir > 30 && FFT_newvirs_gen && count_new_vir*IS_update_coords_subset.size() >pow(10.,6.)) FFT_newvir=true; // FFT with too unbalanced and non uniform M and N works worse. Depends on c. Possibly better combining with far field  && count_new_vir > IS_coords_subset.size()/10000
   
    
    
    //scheme 2
    double eps_B=(0.25*recog_width_/maxradius_tmp)/4.; /**< parameter for kernel rescale rec_width   */
    eps_B=0.1;
    if((0.25*recog_width_/maxradius_tmp)<0.05) eps_B=0.;

    double eps_B_newvirs=(0.25*recog_width_/maxradius_tot_tmp)/4.; /**< parameter for kernel rescale rec_width   */
    eps_B_newvirs=0.1;
    if((0.25*recog_width_/maxradius_tot_tmp)<0.05) eps_B_newvirs=0.;
    
    const double maxradius= maxradius_tmp/(1. - (1/0.25)*eps_B/2.);
    const double maxradius_tot= maxradius_tot_tmp/(1. - (1/0.25)*eps_B_newvirs/2.);

        
    int d= DIM_; /**< number of dimensions    */
    //int N=IS_update_coords_.size(); /**< number of source nodes  */
    int N=IS_update_coords_subset.size(); /**< number of source nodes  */
    int N_newvirs=IS_coords_subset.size(); /**< number of source nodes  */
    //int M= count_old_vir; /**< number of target nodes  */
    int M= count_old_vir; /**< number of target nodes  */
    int M_newvirs= count_new_vir; /**< number of target nodes  */
    int m=5; /**< cut-off parameter       */
    double c=0.25*recog_width_/maxradius; /**< parameter for kernel rescale rec_width   */
    double c_tot=0.25*recog_width_/maxradius_tot; /**< parameter for kernel rescale rec_width   */
    //scheme 1
    //int n=min({2*int((100 + int(5/c))/2),  2*N/2, 500}) ; /**< expansion degree , tune it as function of c, N and M       */
    //int n_newvirs=min({2*int((100 + int(5/c_tot))/2), 2*N_newvirs/2, 500}); /**< expansion degree , tune it as function of c, N and M       */
    
    // scheme 2
    int n=min({2*int((100 + int(5/c))/2),  2*N/2, 300}) ; /**< expansion degree , tune it as function of c, N and M       */
    int n_newvirs=min({2*int((100 + int(5/c_tot))/2), 2*N_newvirs/2, 300}); /**< expansion degree , tune it as function of c, N and M       */
    
    
    
    
    double eps_I=(float(n)/N)*100./c; /**< parameter for near field  */
    double eps_I_max=min(0.01 + max(log10(100./c),0.)*0.01, 0.03); /**< parameter for near field  */
    if(eps_I < 0.01 && c>100) eps_I=0.;
    else if(eps_I < 0.01) eps_I=0.01;
    if(eps_I > eps_I_max) eps_I=eps_I_max;
    
    double eps_I_newvirs=(float(n_newvirs)/N_newvirs)*100./c_tot; /**< parameter for near field  */
    double eps_I_max_newvirs=min(0.01 + max(log10(100./c_tot),0.)*0.01, 0.03); /**< parameter for near field  */
    if(eps_I_newvirs  < 0.01 && c_tot>100) eps_I_newvirs=0.;
    else if(eps_I_newvirs < 0.01) eps_I_newvirs=0.01;
    if(eps_I_newvirs > eps_I_max_newvirs) eps_I_newvirs=eps_I_max_newvirs;
    
    //int p=10; /**< expansion degree , tune it as function of c, N and M       */
    
    int p=2*n*eps_I; /**< expansion degree , tune it as function of c, N and M       */
    if(p<2) p=2;
    if(p>12) p=12; 
    
    int p_newvirs=2*n_newvirs*eps_I_newvirs; /**< expansion degree , tune it as function of c, N and M       */
    if(p_newvirs<2) p_newvirs=2; 
    if(p_newvirs>12) p_newvirs=12; 
    
    fastsum_plan my_fastsum_plan; /**< plan for fast summation */
    
    fastsum_plan my_fastsum_plan_newvirs; /**< plan for fast summation */
    
    kernel k_nfft;
    
	if (kernel_=="exponential"){
	    
	    k_nfft=kern_nfft_exp;
	}
	else{
		cout << "FFT kernel not implemented yet, for now just exponential!" << endl;
		exit(1);
	    
	}
    
    
    if(N==0){
        
        vector<double> new_fitn_fft (vir_coords_.size(), 0.); 
        
        for(int j=0; j < vir_coords_.size(); j++){ 
            
            int x1 =  vir_coords_[j][0] - start_x_;
            int y1 =  vir_coords_[j][1] - start_y_;
            
            double fitness_prec=  world_[x1][y1].get_fitness();
            
            new_fitn_fft[j] = fitness_prec;
               
            
        } 
        
        return new_fitn_fft;
    }
    
    
    //~ cout <<   " FFT PARAMETERS:  N "   <<  N  << " M "   <<  M   << " n "   <<  n    << " c "   <<  c   << " p "   <<  p   << " eps_I "   <<  eps_I    << " eps_B "   <<  eps_B   << endl;
    //~ cout <<   " newvirs:  N_newvirs "   <<  N_newvirs  << " M_newvirs "   <<  M_newvirs   << " n_newvirs "   <<  n_newvirs    << " c "   <<  c_tot  << " p "   <<  p_newvirs   << " eps_I_newvirs "   <<  eps_I_newvirs    << " eps_B_newvirs "   <<  eps_B_newvirs   << endl;
    
    //~ outstream_param << fixed  << setprecision(0)<<setw(19) << N << setprecision(0)<<setw(19) << M <<setprecision(0)<<setw(19) << n  <<setprecision(0)<<setw(19) << m <<setprecision(10)<<setw(19) << c<<setprecision(0)<<setw(19) << p  <<setprecision(10)<<setw(19) << eps_I <<setprecision(10)<<setw(19) << eps_B ;
    
    //~ if (FFT_newvir)     outstream_param << fixed  << setprecision(0)<<setw(19) << N_newvirs << setprecision(0)<<setw(19) << M_newvirs <<setprecision(0)<<setw(19) << n_newvirs  <<setprecision(0)<<setw(19) << m <<setprecision(10)<<setw(19) << c_tot <<setprecision(0)<<setw(19) << p_newvirs <<setprecision(10)<<setw(19) << eps_I_newvirs <<setprecision(10)<<setw(19) << eps_B_newvirs ;
    //~ else    outstream_param << fixed  << setprecision(0)<<setw(19) << -1 << setprecision(0)<<setw(19) << -1 <<setprecision(0)<<setw(19) << -1  <<setprecision(0)<<setw(19) << -1 <<setprecision(10)<<setw(19) << -1 <<setprecision(0)<<setw(19) << -1  <<setprecision(0)<<setw(19) << -1 <<setprecision(10)<<setw(19) << -1 ;
    
    

    if(swap_tot_upd_FFT){
        
        if(FFT_newvir)
        {
            cout << "error, NFFT swap_tot_upd_FFT FFT_newvir must be false! " << endl;
            exit(1); 
        }
        
        if(count_new_vir >500)
        {
            cout << "error, NFFT swap_tot_upd_FFT count_new_vir too big! " << endl;
            exit(1); 
        }
        
        if(count_old_vir < 20 || count_old_vir*IS_update_coords_subset.size() < pow(10.,6.))
        {
            cout << "error, NFFT swap_tot_upd_FFT newvirs not worth FFTing! " << endl;
            exit(1); 
        }
        if( 0.1*count_new_vir*IS_coords_subset.size() > count_old_vir*IS_update_coords_subset.size() )
        {
            cout << "error, NFFT swap_tot_upd_FFT newvirs too little fraction of total! " << endl;
            exit(1); 
        }
        
    }    
    
    
    
    fastsum_init_guru(&my_fastsum_plan, d, N, M, k_nfft, &c, 0, n, m, p, eps_I, eps_B);
  
    //fastsum_init_guru(&my_fastsum_plan, d, N, M,  c, kernel_, 0, n, m);
    if (FFT_newvir)
        fastsum_init_guru(&my_fastsum_plan_newvirs, d, N_newvirs, M_newvirs, k_nfft,  &c_tot, 0, n_newvirs, m, p_newvirs, eps_I_newvirs, eps_B_newvirs);
        //fastsum_init_guru(&my_fastsum_plan_newvirs, d, N_newvirs, M_newvirs,  c_tot, kernel_, 0, n_newvirs, m);
 
 
 
    //~ cout <<  " initialize vir coords rescaled "   << endl;
    

    int count_old_vir_curr=0;
    int count_new_vir_curr=0;
    
    
    for(int j=0; j < vir_coords_.size(); j++){ 
        
        
        double r_max = 0.25;
        double r2 = 0.0;
    
        int x1 =  vir_coords_[j][0] - start_x_;
        int y1 =  vir_coords_[j][1] - start_y_;
        
        double fitness_prec=  world_[x1][y1].get_fitness();
           
        
        double x1_real_resc= 0.25*(idx_to_real_x(x1) - center_rect_x)/maxradius;
        double y1_real_resc= 0.25*(idx_to_real_y(y1) - center_rect_y)/maxradius;
           
        
        double x1_real_resc_tot= 0.25*(idx_to_real_x(x1) - center_rect_x_tot)/maxradius_tot;
        double y1_real_resc_tot= 0.25*(idx_to_real_y(y1) - center_rect_y_tot)/maxradius_tot;
        
        if(swap_tot_upd_FFT){
            if (fitness_prec !=0) fitness_prec=0;
            else fitness_prec=-5;
        }
           
                     
        if (fitness_prec !=0 || approx_fullIS){ // virus was already there, or I don't care
            
            
       
           
            my_fastsum_plan.y[count_old_vir_curr * d + 0] = x1_real_resc;
            my_fastsum_plan.y[count_old_vir_curr * d + 1] = y1_real_resc;
        
            r2 += pow(my_fastsum_plan.y[count_old_vir_curr * d + 0],2.)+pow(my_fastsum_plan.y[count_old_vir_curr * d + 1],2.);
              
            count_old_vir_curr++;
       
        }
        
        else if (FFT_newvir){

           
            my_fastsum_plan_newvirs.y[count_new_vir_curr * d + 0] = x1_real_resc_tot;
            my_fastsum_plan_newvirs.y[count_new_vir_curr * d + 1] = y1_real_resc_tot;
        
            r2 += pow(my_fastsum_plan_newvirs.y[count_new_vir_curr * d + 0],2.)+pow(my_fastsum_plan_newvirs.y[count_new_vir_curr * d + 1],2.);
            
            count_new_vir_curr++;

        }
        
        
        
        
        if (r2 > r_max * r_max)
        {
            cout << "wrong rescaling of vir coords! "<< r2 << " x2_real_resc " << x1_real_resc << " y2_real_resc " << y1_real_resc << " x2 " << idx_to_real_x(x1) << " y2 " << idx_to_real_y(y1)<< " max_x_real " << max_x_real << " max_y_real " << max_y_real << " min_x_real " << min_x_real << " min_y_real " << min_y_real  << " fitness prec " << fitness_prec << endl;
            exit(1);
        } 
    } 

    //~ cout <<  " initialize IS coords rescaled "   << endl;

    if (FFT_newvir){
        for(int i=0; i < IS_coords_subset.size(); i++){ 
            
            
            double r_max = 0.25;
            double r2 = 0.0;
         
            long int x2 =  IS_coords_subset[i][0] - start_x_;
            long int y2 =  IS_coords_subset[i][1] - start_y_;
            
            long int num_IS =  world_[x2][y2].get_num_IS();
            
            
            my_fastsum_plan_newvirs.alpha[i] = (num_IS);///float(people_number_)
            
            
            double x2_real_resc= 0.25*(idx_to_real_x(x2) - center_rect_x_tot)/maxradius_tot;
            double y2_real_resc= 0.25*(idx_to_real_y(y2) - center_rect_y_tot)/maxradius_tot;
                
            my_fastsum_plan_newvirs.x[i * d + 0] = x2_real_resc;
            my_fastsum_plan_newvirs.x[i * d + 1] = y2_real_resc;
        
            r2 += pow(my_fastsum_plan_newvirs.x[i * d + 0],2.)+pow(my_fastsum_plan_newvirs.x[i * d + 1],2.);
              
            if (r2 > r_max * r_max)
            {
                cout << "wrong rescaling of IS  coords! "<< r2 << " x2_real_resc " << x2_real_resc << " y2_real_resc " << y2_real_resc << " x2 " << idx_to_real_x(x2) << " y2 " << idx_to_real_y(y2)<< " max_x_real " << max_x_real << " max_y_real " << max_y_real << " min_x_real " << min_x_real << " min_y_real " << min_y_real  << endl;
                exit(1);
            } 
            
        }
    }
 

    //~ cout <<  " initialize IS update rescaled "   << endl;


    for(int i=0; i < IS_update_coords_subset.size(); i++){ 
     
           
        
        double r_max = 0.25;
        double r2 = 0.0;
     
        long int x2 =  IS_update_coords_subset[i][0] - start_x_;
        long int y2 =  IS_update_coords_subset[i][1] - start_y_;
    
        long int num_upd;
        if(approx_fullIS || swap_tot_upd_FFT)  num_upd = world_[x2][y2].get_num_IS();
        else num_upd =  IS_update_coords_subset[i][2];
    
        
        
        my_fastsum_plan.alpha[i] = (num_upd);// /float(people_number_)
        
        
        double x2_real_resc= 0.25*(idx_to_real_x(x2) - center_rect_x)/maxradius;
        double y2_real_resc= 0.25*(idx_to_real_y(y2) - center_rect_y)/maxradius;
            
        my_fastsum_plan.x[i * d + 0] = x2_real_resc;
        my_fastsum_plan.x[i * d + 1] = y2_real_resc;
    
        r2 += pow(my_fastsum_plan.x[i * d + 0],2.)+pow(my_fastsum_plan.x[i * d + 1],2.);
          
        if (r2 > r_max * r_max)
        {
            cout << "wrong rescaling of IS update coords! "<< r2 << " x2_real_resc " << x2_real_resc << " y2_real_resc " << y2_real_resc << " x2 " << idx_to_real_x(x2) << " y2 " << idx_to_real_y(y2)<< " max_x_real " << max_x_real << " max_y_real " << max_y_real << " min_x_real " << min_x_real << " min_y_real " << min_y_real << endl;
            cout << "IS update "<< i << " IS_update_coords_subset.size() " << IS_update_coords_subset.size() << " x2 " << x2 << " y2 " << y2 << " num_upd " << num_upd  << endl;
            cout << start_x_ << "   " << start_y_  << endl;
            exit(1);
        } 
         
    }
    
 
    //~ chrono::steady_clock::time_point end = std::chrono::steady_clock::now();  

  
    //~ cout << "Time difference = " << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;


 

    //~ /** precomputation */
    //~ cout << "pre-computation:    " << endl;

    begin = std::chrono::steady_clock::now();

    //t0 = getticks();
    fastsum_precompute(&my_fastsum_plan);
    if (FFT_newvir)
        fastsum_precompute(&my_fastsum_plan_newvirs);
    //~ end = std::chrono::steady_clock::now();
    
    //~ cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;

 
 
     //~ cout << "fast computation:   " << endl;
     
     begin = std::chrono::steady_clock::now();
     
     fastsum_trafo(&my_fastsum_plan);
     if (FFT_newvir)
         fastsum_trafo(&my_fastsum_plan_newvirs);
         
    //~ end = std::chrono::steady_clock::now();
    
    //~ cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;

     
     
     //~ cout << "compute fitness vector:   " << endl;
      
      begin = std::chrono::steady_clock::now();

     vector<double> new_fitn_fft (vir_coords_.size(), 0.); // in theory if preserve order I do not need to transform back in real space. BUT  I DO NEED TO KNOW WHAT COORDINATES FALL INTO WHAT SET OF VIRUSES
     //I could directly set fitness here
 
 
    const double center_rect_x_vir= (min_vir_x_ + (max_vir_x_-min_vir_x_)/2.)*latt_width_;
    const double center_rect_y_vir= (min_vir_y_ + (max_vir_y_-min_vir_y_)/2.)*latt_width_;
  
    const double center_rect_x_IS= (min_IS_x_ + (max_IS_x_-min_IS_x_)/2.)*latt_width_;
    const double center_rect_y_IS= (min_IS_y_ + (max_IS_y_-min_IS_y_)/2.)*latt_width_;
     
    const double x2_ref_center_IS_notnorm = center_rect_x_IS - center_rect_x_vir; //center rect IS direction, unit cector
    const double y2_ref_center_IS_notnorm = center_rect_y_IS - center_rect_y_vir;
    const double x2_ref_center_IS = x2_ref_center_IS_notnorm/(sqrt(pow(x2_ref_center_IS_notnorm,2.) + pow(y2_ref_center_IS_notnorm,2.)));
    const double y2_ref_center_IS = y2_ref_center_IS_notnorm/(sqrt(pow(x2_ref_center_IS_notnorm,2.) + pow(y2_ref_center_IS_notnorm,2.)));
    
    const double theta_IS_ref=atan2(y2_ref_center_IS,x2_ref_center_IS);   

 
 
    for (int j = 0; j < my_fastsum_plan.M_total; j++)
    {
        double conv_FFT=abs(my_fastsum_plan.f[j]);
        
        int idx_coord = idx_old_vir[j];
        
        if( idx_new_vir.size()==0 && idx_coord!=j){
            cout << "wrong idx_old_vir!" <<     j   <<   "  "   <<  idx_coord<< endl;
            exit(1);
        }
        
        
        int x1 =  vir_coords_[idx_coord][0] - start_x_;
        int y1 =  vir_coords_[idx_coord][1] - start_y_;
        
        
        double fitness_prec=  world_[x1][y1].get_fitness();
        
        double conv_tot_FFT= conv_FFT;            
            
        if (combined_prep && thetabins){
   
            
            for(int i=0; i < conv_precomp_update.size(); i++){ 
                double theta_bin=conv_precomp_update[i][0];
                
                if (theta_bin!=100){
            
                    double x1_ref_center_rv = idx_to_real_x(x1) - center_rect_x_vir; //specific IS direction
                    double y1_ref_center_rv = idx_to_real_y(y1) - center_rect_y_vir;
                    
                    double dist_vir_center_rv= distance_real(x1_ref_center_rv, y1_ref_center_rv, 0., 0.);
                   
                    //theta w.r.t. IS ref
                    // project change basis to ref and find theta in that basis
                    double x2_ref_center_basis_ISc = x1_ref_center_rv* x2_ref_center_IS + y1_ref_center_rv *y2_ref_center_IS ; 
                    double y2_ref_center_basis_ISc = -x1_ref_center_rv* y2_ref_center_IS + y1_ref_center_rv *x2_ref_center_IS ;
                    
                    double theta;
                    if(y2_ref_center_basis_ISc==0 && x2_ref_center_basis_ISc==0) theta=0; 
                    
                    else theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
    
                    //double theta=atan2(y1_ref_center_rv,x1_ref_center_rv);
                    
                    double x2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*cos(theta_bin);//ptojection on bin average direction, in the rotated basis, but only needed for distance, basis choice does not matter (orthonormal)
                    double y2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*sin(theta_bin);
                    
                    int sign_proj= 1;
                    if (cos(theta - theta_bin)<0) sign_proj= -1;
                    
                    double dist_vir_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.); // parallel bin
                    
                    if (isnan(conv_precomp_update[i][2])){
                        cout << "nan precomp update!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                        exit(1);
                        
                    }
                    
                    //cout << "nan precomp update! " <<     conv_precomp_update[i][2] <<endl;
                  
                    conv_tot_FFT += conv_precomp_update[i][2]*exp(sign_proj*dist_vir_center/recog_width_);
                                
                    if (isnan(conv_tot_FFT)){
                        cout << "nan conv precomp update!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                        exit(1);
                        
                    }
                }

            }
        }
        else if (combined_prep && !thetabins) conv_tot_FFT += conv_precomp_update[0][0];
        
        conv_tot_FFT=conv_tot_FFT/float(people_number_);
        
        if(!approx_fullIS && !swap_tot_upd_FFT) conv_tot_FFT+= fitness_inv_fct(fitness_prec, mem_points_, F0_);

  
        double fitn_FFT = fitness_fct(conv_tot_FFT, mem_points_, F0_);
        
        if (fitn_FFT<-1) fitn_FFT = -1;
        
        new_fitn_fft[idx_coord] += fitn_FFT;
        
        if(upd_fitn) world_[x1][y1].set_fitness(fitn_FFT);
        
    }
             
    if (FFT_newvir){
        for (int j = 0; j < my_fastsum_plan_newvirs.M_total; j++)
        {
            double conv_FFT=abs(my_fastsum_plan_newvirs.f[j]);
            
            int idx_coord = idx_new_vir[j];
            
            int x1 =  vir_coords_[idx_coord][0] - start_x_;
            int y1 =  vir_coords_[idx_coord][1] - start_y_;
            
            
            if(combined_prep_newvirs && thetabins){
                
                for(int i=0; i < conv_precomp_tot.size(); i++){ 
                    double theta_bin=conv_precomp_tot[i][0];
                        if (theta_bin!=100){
                
                        double x1_ref_center_rv = idx_to_real_x(x1) - center_rect_x_vir; //specific IS direction
                        double y1_ref_center_rv = idx_to_real_y(y1) - center_rect_y_vir;
                        
                        double dist_vir_center_rv= distance_real(x1_ref_center_rv, y1_ref_center_rv, 0., 0.);
                        
                        //theta w.r.t. IS ref
                        // project change basis to ref and find theta in that basis
                        double x2_ref_center_basis_ISc =  x1_ref_center_rv* x2_ref_center_IS + y1_ref_center_rv *y2_ref_center_IS; 
                        double y2_ref_center_basis_ISc = -x1_ref_center_rv* y2_ref_center_IS + y1_ref_center_rv *x2_ref_center_IS;
                        
                        double theta;
                        if(y2_ref_center_basis_ISc==0 && x2_ref_center_basis_ISc==0) theta=0; 
                        
                        else theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
    
                        //double theta=atan2(y1_ref_center_rv,x1_ref_center_rv);
                        
                        double x2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*cos(theta_bin);//ptojection on bin average direction
                        double y2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*sin(theta_bin);
                        
                        int sign_proj= 1;
                        if (cos(theta - theta_bin)<0) sign_proj= -1;
                        
                        double dist_vir_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.); // parallel bin
              
                        if (isnan(conv_precomp_tot[i][2])){
                            cout << "nan precomp tot!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                            exit(1);
                            
                        }
                        
                        //cout << "nan precomp tot! " <<     conv_precomp_tot[i][2] <<endl;
                        
                        conv_FFT += conv_precomp_tot[i][2]*exp(sign_proj*dist_vir_center/recog_width_);
                        
                        if (isnan(conv_FFT)){
                            cout << "nan conv precomp tot!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                            exit(1);
                            
                        }
                    }
                }        
            }
            else if (combined_prep_newvirs && !thetabins) conv_FFT += conv_precomp_tot[0][0];
            
            conv_FFT=conv_FFT/float(people_number_);

            double fitn_FFT = fitness_fct(conv_FFT, mem_points_, F0_);
            
            if (fitn_FFT<-1) fitn_FFT = -1;
            
            new_fitn_fft[idx_coord] += fitn_FFT;
            
            
            
            if(upd_fitn) world_[x1][y1].set_fitness(fitn_FFT);
            
        }
    }
    
    else{
            
     
        for(int j=0; j < idx_new_vir.size(); j++){  // if approx_fullIS empty vector
            
            double conv=0.;
            double fitn=0.;
        
            int idx_coord = idx_new_vir[j];
            
            int x1 =  vir_coords_[idx_coord][0] - start_x_;
            int y1 =  vir_coords_[idx_coord][1] - start_y_;
            
            double fitness_prec=  world_[x1][y1].get_fitness();
        
            
            for(int i=0; i < IS_coords_subset.size(); i++){ 
             
                int x2 =  IS_coords_subset[i][0] - start_x_;
                int y2 =  IS_coords_subset[i][1] - start_y_;
                
                long int num_upd;
                if(swap_tot_upd_FFT) num_upd =  IS_coords_subset[i][2];
                else  num_upd = world_[x2][y2].get_num_IS();

                conv += (num_upd)*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);
                
            }
            
            if(combined_prep_newvirs && thetabins){
                
                for(int i=0; i < conv_precomp_tot.size(); i++){ 
                    double theta_bin=conv_precomp_tot[i][0];
                        if (theta_bin!=100){
                
                        double x1_ref_center_rv = idx_to_real_x(x1) - center_rect_x_vir; //specific IS direction
                        double y1_ref_center_rv = idx_to_real_y(y1) - center_rect_y_vir;
                        
                        double dist_vir_center_rv= distance_real(x1_ref_center_rv, y1_ref_center_rv, 0., 0.);
                        
                        //theta w.r.t. IS ref
                        // project change basis to ref and find theta in that basis
                        double x2_ref_center_basis_ISc =  x1_ref_center_rv* x2_ref_center_IS + y1_ref_center_rv *y2_ref_center_IS; 
                        double y2_ref_center_basis_ISc = -x1_ref_center_rv* y2_ref_center_IS + y1_ref_center_rv *x2_ref_center_IS;
                        
                        double theta;
                        if(y2_ref_center_basis_ISc==0 && x2_ref_center_basis_ISc==0) theta=0; 
                        
                        else theta=atan2(y2_ref_center_basis_ISc,x2_ref_center_basis_ISc); // it's out of vir rect, therefore can't be 0,0
    
                        //double theta=atan2(y1_ref_center_rv,x1_ref_center_rv);
                        
                        double x2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*cos(theta_bin);//ptojection on bin average direction
                        double y2_ref_center = dist_vir_center_rv*cos(theta - theta_bin)*sin(theta_bin);
                        
                        int sign_proj= 1;
                        if (cos(theta - theta_bin)<0) sign_proj= -1;
                        
                        double dist_vir_center= distance_real(x2_ref_center, y2_ref_center, 0., 0.); // parallel bin
              
                        if (isnan(conv_precomp_tot[i][2])){
                            cout << "nan precomp tot!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                            exit(1);
                            
                        }
                        
                        //cout << "nan precomp tot! " <<     conv_precomp_tot[i][2] <<endl;
                        
                        conv += conv_precomp_tot[i][2]*exp(sign_proj*dist_vir_center/recog_width_);
                        
                        if (isnan(conv)){
                            cout << "nan conv precomp tot!" <<     i   <<   "  "   <<  theta_bin  <<   "  "   <<  theta <<   "  "   <<  dist_vir_center  << endl;
                            exit(1);
                            
                        }
                    }
                }        
            }
            else if (combined_prep_newvirs && !thetabins) conv += conv_precomp_tot[0][0];
            conv= conv/float(people_number_);
            
            if(swap_tot_upd_FFT) conv+= fitness_inv_fct(fitness_prec, mem_points_, F0_);
            
            fitn = fitness_fct(conv, mem_points_, F0_);

            
            
            if (fitn<-1) fitn = -1;
            
            new_fitn_fft[idx_coord] += fitn;
            
            
            if(upd_fitn) world_[x1][y1].set_fitness(fitn);
            
        } 
               
        
    }
             
             

    
     /** finalise the plan */
     fastsum_finalize(&my_fastsum_plan);
     if (FFT_newvir)
         fastsum_finalize(&my_fastsum_plan_newvirs);
         
    return new_fitn_fft;

};


// use NFFT after far field pre-computation
vector<double> antigenic_space::update_fitness_IS_update_combine_farfield_NFFT(bool upd_fitn, bool thetabins, bool combined_prep_newvirs, bool FFT_newvirs_gen, bool approx_fullIS){
    
    vector<vector<long int>> IS_update_coords_notapprox; 
    vector<vector<int>> IS_coords_notapprox;
    double conv_precomp_update=0.;
    double conv_precomp_tot=0.;
    vector<vector<double>> conv_precomp_update_theta;
    vector<vector<double>> conv_precomp_tot_theta;
    vector<int> limits_IS_notapprox_update {0,0,0,0};
    vector<int> limits_IS_notapprox_tot {0,0,0,0};
    int count_old_vir=-1;
    int count_new_vir=-1;
    vector<int> idx_old_vir;
    vector<int> idx_new_vir;
    
    vector<double> new_fitn;
    
    if(thetabins){
            
        update_fitness_IS_update_prep_farfield_thetabins(combined_prep_newvirs, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update_theta, conv_precomp_tot_theta, limits_IS_notapprox_update, limits_IS_notapprox_tot, count_old_vir, count_new_vir, idx_old_vir, idx_new_vir, approx_fullIS);
        
    }
    else{
       
        update_fitness_IS_update_prep_farfield(IS_update_coords_notapprox, conv_precomp_update,limits_IS_notapprox_update, count_old_vir, count_new_vir, idx_old_vir, idx_new_vir, approx_fullIS);
        // instead of this overload at this point it may just make sense to incorporate the two methods with extra bool thetabins
        
        combined_prep_newvirs=false;
        conv_precomp_update_theta= vector<vector<double>>(1,vector<double>(1, conv_precomp_update));
        conv_precomp_tot_theta=vector<vector<double>>(1,vector<double>(1, conv_precomp_tot));
        
    }

    new_fitn = update_fitness_IS_update_NFFT(upd_fitn, true, combined_prep_newvirs, thetabins, FFT_newvirs_gen, approx_fullIS, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update_theta, conv_precomp_tot_theta, limits_IS_notapprox_update, limits_IS_notapprox_tot, count_old_vir, count_new_vir, idx_old_vir, idx_new_vir); 
    // print error directly in NFFT method
    
    
    return new_fitn;
};


// choose convolution algorithm runtime, based on estimated complexity, tracks proxy for fitness approximation error
vector<double> antigenic_space::update_fitness_IS_update_switch(bool upd_fitn, double & maxfiterr_o_std){
    const bool approx_fullIS=false;
    const bool thetabins=true;
    const bool combined_prep = true;
    const bool combined_prep_newvirs =true;
    const bool FFT_newvirs_gen =true;

    int algo=0;
 
    int count_old_vir=0;
    int count_new_vir=0;

    for(int j=0; j < vir_coords_.size(); j++){ 
        
    
        int x1 =  vir_coords_[j][0] - start_x_;
        int y1 =  vir_coords_[j][1] - start_y_;
        
        double fitness_prec=  world_[x1][y1].get_fitness();
                     
        if (fitness_prec !=0 ){ //  
            count_old_vir+=1;
        }
        
        else{
            count_new_vir+=1;
        }
        
    } 
    long int time_compl_full= count_old_vir*IS_update_coords_.size() + count_new_vir*IS_coords_.size();
    
    int time_abs_thr=pow(10., 6.);
    
    vector<vector<long int>> IS_update_coords_notapprox; 
    vector<vector<int>> IS_coords_notapprox;
    //double conv_precomp_update=0.;
    //double conv_precomp_tot=0.;
    vector<vector<double>> conv_precomp_update_theta;
    vector<vector<double>> conv_precomp_tot_theta;
    vector<int> limits_IS_notapprox_update {0,0,0,0};
    vector<int> limits_IS_notapprox_tot {0,0,0,0};
    vector<int> idx_old_vir;
    vector<int> idx_new_vir;
    
    vector<double> new_fitn;
    
    cout << time_compl_full << ' ' << count_old_vir << ' ' << count_new_vir << ' ' << IS_coords_.size() << ' ' << IS_update_coords_.size()  << endl;
    
    if(time_compl_full<=time_abs_thr) new_fitn = update_fitness_fullconv_update(true);
    else{
            
        update_fitness_IS_update_prep_farfield_thetabins(combined_prep_newvirs, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update_theta, conv_precomp_tot_theta, limits_IS_notapprox_update, limits_IS_notapprox_tot, count_old_vir, count_new_vir, idx_old_vir, idx_new_vir, approx_fullIS);
        
        long int size_subs_upd=IS_update_coords_notapprox.size();
        long int size_subs_tot=IS_coords_notapprox.size();
    
        
        if (conv_precomp_update_theta.size()==0 || IS_update_coords_notapprox.size()==0){
            size_subs_upd=IS_update_coords_.size();
        }
        if (conv_precomp_tot_theta.size()==0 || IS_coords_notapprox.size()==0){
            size_subs_tot=IS_coords_.size();
        }
        long int time_compl_farf_conv= count_old_vir*size_subs_upd + count_new_vir*size_subs_tot;
        
        cout << time_compl_farf_conv << ' ' << count_old_vir << ' ' << count_new_vir << ' ' << size_subs_upd << ' ' << size_subs_tot << endl;
        
        int num_newvir_thr=30;
        
        //scheme 1
        //int num_oldvir_thr=500;
        
        //scheme 2
        int num_oldvir_thr=300;
        
        bool swap_tot_upd_FFT= false; // in this condition is convenient to do FFT only on newvirs
        
        if(time_compl_farf_conv<=time_abs_thr || count_old_vir < num_oldvir_thr){
            new_fitn = update_fitness_fullconv_update_precomp_thetabins(upd_fitn, combined_prep, combined_prep_newvirs, approx_fullIS, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update_theta, conv_precomp_tot_theta);  
            
            algo=1; 
        } 
        else{
            new_fitn = update_fitness_IS_update_NFFT(upd_fitn, combined_prep, combined_prep_newvirs, thetabins, FFT_newvirs_gen, approx_fullIS, IS_update_coords_notapprox, IS_coords_notapprox, conv_precomp_update_theta, conv_precomp_tot_theta, limits_IS_notapprox_update, limits_IS_notapprox_tot, count_old_vir, count_new_vir, idx_old_vir, idx_new_vir); 
            
            algo=2; 
       
        }   
        
    }
    
    
    
    vector<int> idxs_max= get_maxobs(); //num virs, IS, fitn
    
    vector<double> FFT_errors; //error fullconv to do vs actual. contains for each entry: fitn1, fitn2, abs(fitn1 - fitn2). Total 9 elements 
    
    for(int j=0; j < idxs_max.size(); j++){ 
        
        int idx = idxs_max[j];
        
        int x1 =  vir_coords_[idx][0] - start_x_;
        int y1 =  vir_coords_[idx][1] - start_y_;
        
        double fitness_FFT=  new_fitn[idx];
        
        FFT_errors.push_back(fitness_FFT);
        
    
        double conv=0.;
        double fitn=0.;
    
        
        for(int i=0; i < IS_coords_.size(); i++){ 
         
            int x2 =  IS_coords_[i][0] - start_x_;
            int y2 =  IS_coords_[i][1] - start_y_;
            
            conv += (world_[x2][y2].get_num_IS())*kern(distance(x1, y1, x2, y2), recog_width_, kernel_);
            
        }
        
        conv= conv/float(people_number_);

        fitn = fitness_fct(conv, mem_points_, F0_);

        
        
        if (fitn<-1) fitn = -1;
        
         FFT_errors.push_back(fitn);
         
         FFT_errors.push_back(abs(fitness_FFT- fitn));
        
       
    
        //~ outstream << fixed  << setprecision(10)<<setw(19) << fitness_FFT << setprecision(10)<<setw(19) << fitn <<setprecision(10)<<setw(19) << abs(fitness_FFT- fitn) ;

        
    }
    
    
    
	vector<double> fitn_stat = get_avg_fitness();    
	//double avg_fitness = fitn_stat[0];    
	double std_fitness = sqrt(fitn_stat[1]);    
    
    maxfiterr_o_std=max({FFT_errors.back(), FFT_errors[2], FFT_errors[5]})/std_fitness;
    
    
    //~ outstream  <<endl;
    
    return new_fitn;

};

// time and errors benchmark function comparing all convolution algorithms, in the end updates using the full convolution over all populated sites. Call this every 10000 cycles

double antigenic_space::benchmark_update_fitness_fullconv_update_NFFT_farfield(){
    cout << endl;
    cout <<  " update_fitness_fullconv update FFT benchmark "   << endl;
    
    bool approx_fullIS= false;
     
    
    string times_b_f ("times_errors_benchmark_conv.dat");
    

    string outfile= out_dir_+ times_b_f;

    fstream outstream;
    
    outstream.open(outfile, fstream::out | fstream::app); 
     
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();  
     
     vector<double> new_fitn_fft =  update_fitness_IS_update_NFFT(false, false, false, false , true, approx_fullIS); // in theory if preserve order I do not need to transform back in real space. BUT  I DO NEED TO KNOW WHAT COORDINATES FALL INTO WHAT SET OF VIRUSES
     // fifth bool for FFT newvirs
     
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();  
    cout << "TIME difference tot NFFT = " << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
    
    outstream << fixed  << setprecision(5)<<setw(19) << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  ;

     
    cout << endl; 
     

     begin = std::chrono::steady_clock::now();  

     vector<double> new_fitn_farfield =  update_fitness_IS_update_combine_farfield_fullconv(false, true, true, approx_fullIS); // second bool decides if theta bins, third if farfield only on update, fourth approx_fullIS
     

    end = std::chrono::steady_clock::now();  
    cout << "TIME difference tot FARFIELD = " << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
     
     outstream << fixed  << setprecision(5)<<setw(19) << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  ;
     
    cout << endl; 
     
     begin = std::chrono::steady_clock::now();       
          
     vector<double> new_fitn_fullconv_upd =  update_fitness_fullconv_update(false); // in theory if preserve order I do not need to transform back in real space. BUT  I DO NEED TO KNOW WHAT COORDINATES FALL INTO WHAT SET OF VIRUSES
 


    end = std::chrono::steady_clock::now();  
    cout << "TIME difference tot FULLCONV = " << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
    
    double fullconv_time=(chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0 ;
    cout << endl; 
    
    
     begin = std::chrono::steady_clock::now();       
     
     vector<double> new_fitn_farfield_nfft =  update_fitness_IS_update_combine_farfield_NFFT(false, true, true, true, approx_fullIS); // second bool decides if theta bins, third sets if all IS coords approximated by farfield theta or just the update (if not, thetabins automathically set to false), fourth if FFT_newvirs, fifth if approximate directly all IS coords
     


    end = std::chrono::steady_clock::now();  
    cout << "TIME difference tot FARFIELD NFFT = " << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
    double farf_nfft_time=(chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0 ;

    outstream << fixed  << setprecision(5)<<setw(19) << farf_nfft_time  ;
     
    cout << endl; 
    

     
     begin = std::chrono::steady_clock::now();       
     
     double maxfiterr_o_std=0;
          
     //vector<double> new_fitn_switch =  update_fitness_IS_update_switch(true); // in theory if preserve order I do not need to transform back in real space. BUT  I DO NEED TO KNOW WHAT COORDINATES FALL INTO WHAT SET OF VIRUSES
     vector<double> new_fitn_switch =  update_fitness_IS_update_switch(false, maxfiterr_o_std); // in theory if preserve order I do not need to transform back in real space. BUT  I DO NEED TO KNOW WHAT COORDINATES FALL INTO WHAT SET OF VIRUSES
 


    end = std::chrono::steady_clock::now();  
    cout << "TIME difference tot SWITCH  = " << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
    
    double switch_time=(chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0 ;
    

    cout << endl; 
    
     
    outstream << fixed  << setprecision(5)<<setw(19) << fullconv_time  ;
    
    begin = std::chrono::steady_clock::now(); 
    
    //new_fitn_fullconv_upd =  update_fitness_fullconv(false); // if update by approximation, errors must be calculated w.r.t. the true fitness
    new_fitn_fullconv_upd =  update_fitness_fullconv(true); // if update by approximation, errors must be calculated w.r.t. the true fitness
     
    end = std::chrono::steady_clock::now();  
    cout << "TIME difference tot FULLCONV FULL = " << (chrono::duration_cast<chrono::nanoseconds>(end - begin).count()) /1000000000.0  << "[s]" << endl;
     
     /** compute max error */
     double error = 0.; /**< for error computation   */
     double error_farfield = 0.; /**< for error computation   */
     double error_farfield_nfft = 0.; /**< for error computation   */
     double error_switch = 0.; /**< for error computation   */
           
     double conv_error = 0.; /**< for error computation   */
     double conv_error_farfield = 0.; /**< for error computation   */
     double conv_error_farfield_nfft = 0.; /**< for error computation   */
     double conv_error_switch = 0.; /**< for error computation   */
           
     double full_conv_avg        =0.;
     double full_conv_sq        =0.;
           
     double full_fitn_avg        =0.;
     double full_fitn_sq        =0.;
     
     double direct_tot_avg          =0.;
     double fast_tot_avg          =0.;
     double farfield_tot_avg      =0.;
     double farfield_nfft_tot_avg=0.;
     double switch_tot_avg=0.;
     
     double fast_err_avg          =0.;
     double farfield_err_avg      =0.;
     double farfield_nfft_err_avg=0.;
     double switch_err_avg=0.;
     
     double fast_err_conv_avg          =0.;
     double farfield_err_conv_avg      =0.;
     double farfield_nfft_err_conv_avg=0.;
     double switch_err_conv_avg=0.;
     
     double fast_err_conv_resc_avg          =0.;
     double farfield_err_conv_resc_avg      =0.;
     double farfield_nfft_err_conv_resc_avg=0.;
     double switch_err_conv_resc_avg=0.;
     
     long int tot_vir=0;
           
     double direct_tot=0.;
     double fast_tot=0.;
     double farfield_tot=0.;
     double farfield_nfft_tot=0.;
     double switch_tot=0.;
     for (int j = 0; j < vir_coords_.size(); j++)
     {
       
       double fitness_fullconv=  new_fitn_fullconv_upd[j];
       double fitness_FFT=  new_fitn_fft[j];
       double fitness_farfield=  new_fitn_farfield[j];
       double fitness_farfield_nfft=  new_fitn_farfield_nfft[j];
       double fitness_switch=  new_fitn_switch[j];
        
       double conv_fullconv     =  fitness_inv_fct(fitness_fullconv, mem_points_, F0_);     
       double conv_FFT          =  fitness_inv_fct(fitness_FFT, mem_points_, F0_);          
       double conv_farfield     =  fitness_inv_fct(fitness_farfield, mem_points_, F0_);     
       double conv_farfield_nfft=  fitness_inv_fct(fitness_farfield_nfft, mem_points_, F0_);
       double conv_switch=  fitness_inv_fct(fitness_switch, mem_points_, F0_);
        
        int x =  vir_coords_[j][0] - start_x_;
        int y =  vir_coords_[j][1] - start_y_;
        long int num_vir=world_[x][y].get_num_vir();
        
        tot_vir+= num_vir;

       
         
       if (abs(fitness_fullconv - fitness_FFT) / abs(fitness_fullconv) > error)
         error = abs(fitness_fullconv - fitness_FFT) / abs(fitness_fullconv);
       if (abs(fitness_fullconv - fitness_farfield) / abs(fitness_fullconv) > error_farfield)
         error_farfield = abs(fitness_fullconv - fitness_farfield) / abs(fitness_fullconv);
       if (abs(fitness_fullconv - fitness_switch) / abs(fitness_fullconv) > error_switch)
         error_switch = abs(fitness_fullconv - fitness_switch) / abs(fitness_fullconv);
       if (abs(fitness_fullconv - fitness_farfield_nfft) / abs(fitness_fullconv) > error_farfield_nfft)
         error_farfield_nfft = abs(fitness_fullconv - fitness_farfield_nfft) / abs(fitness_fullconv);
         
       if (abs(conv_fullconv - conv_FFT) / abs(conv_fullconv) > conv_error)
         conv_error = abs(conv_fullconv - conv_FFT) / abs(conv_fullconv);
       if (abs(conv_fullconv - conv_farfield) / abs(conv_fullconv) > conv_error_farfield)
         conv_error_farfield = abs(conv_fullconv - conv_farfield) / abs(conv_fullconv);
       if (abs(conv_fullconv - conv_switch) / abs(conv_fullconv) > conv_error_switch)
         conv_error_switch = abs(conv_fullconv - conv_switch) / abs(conv_fullconv);
       if (abs(conv_fullconv - conv_farfield_nfft) / abs(conv_fullconv) > conv_error_farfield_nfft)
         conv_error_farfield_nfft = abs(conv_fullconv - conv_farfield_nfft) / abs(conv_fullconv);
         
       direct_tot += abs(fitness_fullconv);
       fast_tot += abs(fitness_FFT);
       farfield_tot += abs(fitness_farfield);
       farfield_nfft_tot += abs(fitness_farfield_nfft);
       switch_tot += abs(fitness_switch);
         
       full_fitn_avg       += num_vir*fitness_fullconv;
       full_fitn_sq       += num_vir*pow(fitness_fullconv,2.);
         
       full_conv_avg       += num_vir*conv_fullconv;
       full_conv_sq       += num_vir*pow(conv_fullconv,2.);
       
       direct_tot_avg       += num_vir*abs(fitness_fullconv);
       fast_tot_avg         += num_vir*abs(fitness_FFT);
       farfield_tot_avg     += num_vir*abs(fitness_farfield);
       farfield_nfft_tot_avg+= num_vir*abs(fitness_farfield_nfft);
       switch_tot_avg+= num_vir*abs(fitness_switch);
         
       fast_err_avg         += num_vir*abs(fitness_fullconv - fitness_FFT);
       farfield_err_avg     += num_vir*abs(fitness_fullconv - fitness_farfield);
       farfield_nfft_err_avg+= num_vir*abs(fitness_fullconv - fitness_farfield_nfft);
       switch_err_avg+= num_vir*abs(fitness_fullconv - fitness_switch);
         
       fast_err_conv_avg         += num_vir*abs(conv_fullconv - conv_FFT);
       farfield_err_conv_avg     += num_vir*abs(conv_fullconv - conv_farfield);
       farfield_nfft_err_conv_avg+= num_vir*abs(conv_fullconv - conv_farfield_nfft);
       switch_err_conv_avg+= num_vir*abs(conv_fullconv - conv_switch);
         
       fast_err_conv_resc_avg         += num_vir*abs(conv_fullconv - conv_FFT)/(1.-conv_FFT/mem_points_);
       farfield_err_conv_resc_avg     += num_vir*abs(conv_fullconv - conv_farfield)/(1.-conv_FFT/mem_points_);
       farfield_nfft_err_conv_resc_avg+= num_vir*abs(conv_fullconv - conv_farfield_nfft)/(1.-conv_FFT/mem_points_);
       switch_err_conv_resc_avg+= num_vir*abs(conv_fullconv - conv_switch)/(1.-conv_FFT/mem_points_);
     }
     
     
     full_fitn_avg       /=tot_vir;
     full_fitn_sq        /=tot_vir;
     
     double full_fitn_stdev=sqrt(full_fitn_sq - pow(full_fitn_avg,2.) );
     
     full_conv_avg       /=tot_vir;
     full_conv_sq        /=tot_vir;
     
     double full_conv_stdev=sqrt(full_conv_sq - pow(full_conv_avg,2.) );
     
     direct_tot_avg       /=tot_vir;
     fast_tot_avg         /=tot_vir;
     farfield_tot_avg     /=tot_vir;
     farfield_nfft_tot_avg/=tot_vir;
     switch_tot_avg/=tot_vir;
   
     fast_err_avg         /=tot_vir;
     farfield_err_avg     /=tot_vir;
     farfield_nfft_err_avg/=tot_vir;
     switch_err_avg/=tot_vir;
   
     fast_err_conv_avg         /=tot_vir;
     farfield_err_conv_avg     /=tot_vir;
     farfield_nfft_err_conv_avg/=tot_vir;
     switch_err_conv_avg/=tot_vir;
   
   
     fast_err_conv_resc_avg         /=tot_vir;
     farfield_err_conv_resc_avg     /=tot_vir;
     farfield_nfft_err_conv_resc_avg/=tot_vir;
     switch_err_conv_resc_avg/=tot_vir;
   
   
        // I could reimplement on world with cutoff at like 7 rec_width_
    cout <<   " direct tot "   <<  direct_tot   << endl;
    cout <<  " fast tot "   <<    fast_tot   << endl;
    cout <<  " farfield tot "   <<    farfield_tot   << endl;
    cout <<  " farfield NFFT tot "   <<    farfield_nfft_tot   << endl;
    cout <<  " switch tot "   <<    switch_tot   << endl;
    cout <<  " max relative error FFT "   <<    error  << endl;
    cout <<  " max relative error farfield "   <<    error_farfield  << endl;
    cout <<  " max relative error farfield NFFT "   <<    error_farfield_nfft  << endl;
    cout <<  " max relative error switch "   <<    error_switch  << endl;
    cout <<  " max relative conv_error FFT "   <<    conv_error  << endl;
    cout <<  " max relative conv_error farfield "   <<    conv_error_farfield  << endl;
    cout <<  " max relative conv_error farfield NFFT "   <<    conv_error_farfield_nfft  << endl;
    cout <<  " max relative conv_error switch "   <<    conv_error_switch  << endl;
    
    outstream << fixed  << setprecision(5)<<setw(19) << direct_tot << setprecision(5)<<setw(19) << fast_tot << setprecision(5)<<setw(19) << farfield_tot << setprecision(5)<<setw(19) << farfield_nfft_tot << setprecision(5)<<setw(19) << error << setprecision(5)<<setw(19) << error_farfield << setprecision(5)<<setw(19) << error_farfield_nfft  ;
    
    outstream << fixed << setprecision(9)<<setw(19) << direct_tot_avg << setprecision(9)<<setw(19) << fast_tot_avg << setprecision(9)<<setw(19) << farfield_tot_avg << setprecision(9)<<setw(19) << farfield_nfft_tot_avg << setprecision(9)<<setw(19) << conv_error << setprecision(9)<<setw(19) << conv_error_farfield << setprecision(9)<<setw(19) << conv_error_farfield_nfft ;
    
    outstream << fixed << setprecision(9)<<setw(19) << fast_err_avg << setprecision(9)<<setw(19) << farfield_err_avg << setprecision(9)<<setw(19) << farfield_nfft_err_avg  << setprecision(9)<<setw(19) << fast_err_conv_avg << setprecision(9)<<setw(19) << farfield_err_conv_avg << setprecision(9)<<setw(19) << farfield_nfft_err_conv_avg << setprecision(9)<<setw(19) << fast_err_conv_resc_avg << setprecision(9)<<setw(19) << farfield_err_conv_resc_avg << setprecision(9)<<setw(19) << farfield_nfft_err_conv_resc_avg << setprecision(9)<<setw(19) << full_fitn_stdev << setprecision(9)<<setw(19) << full_conv_stdev  ;
    
    outstream << fixed  << setprecision(5)<<setw(19) << switch_time  ;

    outstream << fixed << setprecision(5)<<setw(19) << switch_tot << setprecision(5)<<setw(19) << error_switch<< setprecision(9)<<setw(19) << switch_tot_avg << setprecision(9)<<setw(19) << conv_error_switch  << setprecision(9)<<setw(19) << switch_err_avg   << setprecision(9)<<setw(19) << switch_err_conv_avg  << setprecision(9)<<setw(19) << switch_err_conv_resc_avg; 
      
    outstream  <<endl;

    cout << endl;
    cout << endl;
    
    return switch_err_avg/full_fitn_stdev;
    
};


void antigenic_space::update_IS(){
    cout <<  " update_IS "   << endl;
    
    
    
    int max_x_tot=max(max_vir_x_, max_IS_x_);
    int max_y_tot=max(max_vir_y_, max_IS_y_);
    int min_x_tot=min(min_vir_x_, min_IS_x_);
    int min_y_tot=min(min_vir_y_, min_IS_y_);

    
    //bool agent_based=false;
    bool agent_based=false;

    vector<vector<long int>>().swap( IS_update_coords_); // clear IS updates, and memory, to allocate the new ones
    
    long int tot_update=0;
    
    

    vector<long int> new_IS_vec (vir_coords_.size(), 0);
    vector<long int> prec_IS_vec (vir_coords_.size(), 0);
    //vector<vector<int>> new_IS_coord;
    
    for(int j=0; j < vir_coords_.size(); j++){ 
        
        int x =  vir_coords_[j][0] - start_x_;
        int y =  vir_coords_[j][1] - start_y_;
     
        long int num_virs= world_[x][y].get_num_vir();
        long int num_IS  = world_[x][y].get_num_IS();
        
        long int new_IS=num_virs;
        
        if (new_IS + num_IS > mem_points_*people_number_) new_IS = mem_points_*people_number_ - num_IS;
        
        if(new_IS<0){
            cout << "new_IS less than 0! "<< new_IS << endl;
            exit(1); 
        }
        new_IS_vec[j] =new_IS;
        prec_IS_vec[j] =num_IS;
        tot_update +=new_IS;
   
   
        
        long int toupdate=tot_update;
        
        
        
        int min_IS_x_chk= numeric_limits<int>::max();  
        int min_IS_y_chk= numeric_limits<int>::max();  
        int max_IS_x_chk= numeric_limits<int>::min();  
        int max_IS_y_chk= numeric_limits<int>::min();  
    
    
        
        for(int j=0; j < IS_coords_.size(); j++){
            
            //cout <<endl;
            //cout << j <<endl;

            int x =  IS_coords_[j][0] - start_x_;
            int y =  IS_coords_[j][1] - start_y_;
            
            long int num_vir  = world_[x][y].get_num_vir();
            long int num_IS  = world_[x][y].get_num_IS();
            
            double prob_upd=tot_update/float(mem_points_*people_number_);
            
            binomial_distribution<long int> bin_distribution(num_IS,prob_upd);
            
            long int num_upd=bin_distribution(mt);


            
            if (num_upd >toupdate) num_upd = toupdate;
            
            world_[x][y].set_num_IS(num_IS - num_upd);
            
            
            if (num_upd!=0 && num_vir ==0){ // If vir non zero I will add to coords later
            
                vector< long int> newcoords {x + start_x_, y + start_y_, - num_upd};
                
                IS_update_coords_.push_back(newcoords);
                
                
                //cout << "IS update check first loop "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << IS_update_coords_[IS_update_coords_.size()-1][2]  << endl;
                
               if (IS_update_coords_[IS_update_coords_.size()-1][0] > max_x_tot || IS_update_coords_[IS_update_coords_.size()-1][0] < min_x_tot || IS_update_coords_[IS_update_coords_.size()-1][1] > max_y_tot || IS_update_coords_[IS_update_coords_.size()-1][1] < min_y_tot )
                {
                    cout << "wrong IS update coords first loop!  x " << IS_update_coords_[IS_update_coords_.size()-1][0] << " y " << IS_update_coords_[IS_update_coords_.size()-1][1] << " max_x_tot " << max_x_tot << " max_y_tot " << max_y_tot << " min_x_tot " << min_x_tot << " min_y_tot " << min_y_tot << endl;
                    cout << "IS update "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << IS_update_coords_[IS_update_coords_.size()-1][2]  << endl;
                 // Generate an interrupt
                //~ std::raise(SIGINT);

                    exit(1);
                } 
                
            }
            
            
            
        
            if(num_IS==num_upd){// sent to zero the IS, now have to check for virs. I can remove it from IS list, the ones I have to add are in vir coords anyway. Here check if shrink world
                
                           
                if (IS_coords_.size() >= 2)
                {
                    swap(IS_coords_[j], *(IS_coords_.end() - 1));
                    IS_coords_.pop_back();
                    j--;
                }
                else
                {
                    vector<vector<int>>().swap(IS_coords_); //clears memory, including pointers
                }
                
                long int num_vir  = world_[x][y].get_num_vir();
                
                if (num_vir == 0){
                    //shrink space
                    
                    //~ cout <<  " check shrink_world IS loop"   << endl;

                    
                    //checkpad_and_shrink(x, y);
                    
                   
                   // check if moved min or max, and shrink world consequently
                           
                    if (x + start_x_ == min_IS_x_){
                        min_IS_x_ = check_new_border( x, y, 1, false) + start_x_;
                        
                        int min_tot= min(min_IS_x_ , min_vir_x_) - 1;
                        
                        if (num_x_ > pad_size_ && (min_tot - start_x_) >= pad_size_){
                                
                            int size_shrink=(min_tot - start_x_)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(3, size_shrink); // changes start_x_, careful what you do next
                        }
                    }
                    
                    else if (x + start_x_ == max_IS_x_){ 
                        max_IS_x_ = check_new_border( x, y, 2, false) + start_x_;
                        
                        int max_tot= max(max_IS_x_ , max_vir_x_) + 1;
                        
                        if (num_x_ > pad_size_ && (max_tot - start_x_) < num_x_ - pad_size_){
                                
                            int size_shrink=(num_x_ - (max_tot - start_x_) - 1)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(1, size_shrink);
                        }
                        
                    }
                    
                    if (y + start_y_ == min_IS_y_) {
                        min_IS_y_ = check_new_border(x, y, 3, false) + start_y_;
                        
                        int min_tot= min(min_IS_y_ , min_vir_y_) - 1;

                        if (num_y_ > pad_size_ && (min_tot - start_y_) >= pad_size_){
                                
                            int size_shrink=(min_tot - start_y_)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(2, size_shrink);
                        }
                    }
                    
                    
                
                    else if (y + start_y_ == max_IS_y_){ 
                        max_IS_y_ = check_new_border( x, y, 4, false) + start_y_;
                        
                        int max_tot= max(max_IS_y_ , max_vir_y_) + 1;
                        
                        if (num_y_ > pad_size_ && (max_tot - start_y_) < num_y_ - pad_size_){
                                
                            int size_shrink=(num_y_ - (max_tot - start_y_) - 1)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(0, size_shrink);
                
                        }
                        
                    }

                   //~ cout <<  " checked"   << endl;




                }
                
            }
           else{
            
                if(x + start_x_ < min_IS_x_chk) min_IS_x_chk = x + start_x_ ;
                if(y + start_y_ < min_IS_y_chk) min_IS_y_chk = y + start_y_ ;
                if(x + start_x_ > max_IS_x_chk) max_IS_x_chk = x + start_x_ ;
                if(y + start_y_ > max_IS_y_chk) max_IS_y_chk = y + start_y_ ;
                
            }             
            
            toupdate -= num_upd;
 
            //cout << toupdate <<endl;
           
            if (toupdate ==0) break;

            
        }; 
        
        if(toupdate>0){ // I can correct only if actually checked all coords
                       
            if(min_IS_x_chk!=min_IS_x_){ // if sth goes wron in minmax tracking, correct here. It should not
                //~ cout << "correcting min_IS_x_ "<< min_IS_x_ << " to " << min_IS_x_chk << endl;
                min_IS_x_=min_IS_x_chk;
            }
            
            
            if(max_IS_x_chk!=max_IS_x_){ // if sth goes wron in minmax tracking, correct here. It should not
                //~ cout << "correcting max_IS_x_ "<< max_IS_x_ << " to " << max_IS_x_chk << endl;
                max_IS_x_=max_IS_x_chk;
            }
            
            if(min_IS_y_chk!=min_IS_y_){
                //~ cout << "correcting min_IS_y_ "<< min_IS_y_ << " to " << min_IS_y_chk << endl;
                min_IS_y_=min_IS_y_chk;
            }
        
            if(max_IS_y_chk!=max_IS_y_){ // if sth goes wron in minmax tracking, correct here. It should not
                //~ cout << "correcting max_IS_y_ "<< max_IS_y_ << " to " << max_IS_y_chk << endl;
                max_IS_y_=max_IS_y_chk;
            }
            
            
            max_x_tot=max(max_vir_x_, max_IS_x_);
            max_y_tot=max(max_vir_y_, max_IS_y_);
            min_x_tot=min(min_vir_x_, min_IS_x_);
            min_y_tot=min(min_vir_y_, min_IS_y_);
        }
    
        
        
        
        while(toupdate >0){
            //~ cout <<  " draw j"   << endl;
               
            uniform_int_distribution<int> distr_int(0,IS_coords_.size()-1);

            //~ cout <<  " ready"   << endl;

             int j =  distr_int(mt); //segfault if virs > IS
             
             //~ cout <<  " drawn"   << endl;
             
            
            //~ cout << j <<endl;
    
            int x =  IS_coords_[j][0] - start_x_;
            int y =  IS_coords_[j][1] - start_y_;
            //~ cout << x <<endl;

            long int num_vir  = world_[x][y].get_num_vir();

            long int num_IS = world_[x][y].get_num_IS();  
            
            long int toupdate_tmp=toupdate;
            
            if (num_IS > toupdate){
                
                world_[x][y].set_num_IS(num_IS-toupdate);
               
                if (toupdate_tmp!=0  && num_vir ==0){ // I should check if this coord exists already, for now let's not
                
                    vector< long int> newcoords {x + start_x_, y + start_y_, - toupdate_tmp};
                    
                    IS_update_coords_.push_back(newcoords);
                    
                        
                        
                    
                    //cout << "IS update check first loop "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << IS_update_coords_[IS_update_coords_.size()-1][2]  << endl;
                    
                   if (IS_update_coords_[IS_update_coords_.size()-1][0] > max_x_tot || IS_update_coords_[IS_update_coords_.size()-1][0] < min_x_tot || IS_update_coords_[IS_update_coords_.size()-1][1] > max_y_tot || IS_update_coords_[IS_update_coords_.size()-1][1] < min_y_tot )
                    {
                        cout << "wrong IS update coords draw loop last upd!  x " << IS_update_coords_[IS_update_coords_.size()-1][0] << " y " << IS_update_coords_[IS_update_coords_.size()-1][1] << " max_x_tot " << max_x_tot << " max_y_tot " << max_y_tot << " min_x_tot " << min_x_tot << " min_y_tot " << min_y_tot << endl;
                        cout << "IS update "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << IS_update_coords_[IS_update_coords_.size()-1][2]  << endl;
                        exit(1);
                    } 
                    
                    
                }
                
                
                toupdate=0;
            }
            
                 
            else{// sent to zero the IS, now have to check for virs. I can remove it from IS list, the ones I have to add are in vir coords anyway. Here check if shrink world
                
                toupdate-=num_IS;
                world_[x][y].set_num_IS(0);

                if (num_IS!=0  && num_vir ==0 ){ // I should check if this coord exists already, for now let's not
                
                    vector< long int> newcoords {x + start_x_, y + start_y_, - num_IS};
                    
                    IS_update_coords_.push_back(newcoords);
                    
                    
                    
                    
                    
                    //cout << "IS update check first loop "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << newcoords[2]  << endl;
                    
                   if (IS_update_coords_[IS_update_coords_.size()-1][0] > max_x_tot || IS_update_coords_[IS_update_coords_.size()-1][0] < min_x_tot || IS_update_coords_[IS_update_coords_.size()-1][1] > max_y_tot || IS_update_coords_[IS_update_coords_.size()-1][1] < min_y_tot )
                    {
                        cout << "wrong IS update coords draw loop while!  x " << IS_update_coords_[IS_update_coords_.size()-1][0] << " y " << IS_update_coords_[IS_update_coords_.size()-1][1] << " max_x_tot " << max_x_tot << " max_y_tot " << max_y_tot << " min_x_tot " << min_x_tot << " min_y_tot " << min_y_tot << endl;
                        cout << "IS update "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << IS_update_coords_[IS_update_coords_.size()-1][2]  << endl;
                        exit(1);
                    } 
                    
                    
                }
                            
                           
                if (IS_coords_.size() >= 2)
                {
                    swap(IS_coords_[j], *(IS_coords_.end() - 1));
                    IS_coords_.pop_back();
                    j--;
                }
                else
                {
                    vector<vector<int>>().swap(IS_coords_); //clears memory, including pointers
                }
                
                long int num_vir  = world_[x][y].get_num_vir();
                
                if (num_vir == 0){
   
                 
                    if (x + start_x_ == min_IS_x_){
                        min_IS_x_ = check_new_border( x, y, 1, false) + start_x_;
                        
                        int min_tot= min(min_IS_x_ , min_vir_x_) - 1;
                        
                        if (num_x_ > pad_size_ && (min_tot - start_x_) >= pad_size_){
                                
                            int size_shrink=(min_tot - start_x_)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(3, size_shrink); // changes start_x!
                        }
                    }

                    else if (x + start_x_ == max_IS_x_){ 
                        max_IS_x_ = check_new_border( x, y, 2, false) + start_x_;
                        
                        int max_tot= max(max_IS_x_ , max_vir_x_) + 1;
                        
                        if (num_x_ > pad_size_ && (max_tot - start_x_) < num_x_ - pad_size_){
                                
                            int size_shrink=(num_x_ - (max_tot - start_x_) - 1)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(1, size_shrink);
                        }
                        
                    }   
                                     
                    if (y + start_y_ == min_IS_y_) {
                        min_IS_y_ = check_new_border(x, y, 3, false) + start_y_;
                        
                        int min_tot= min(min_IS_y_ , min_vir_y_) - 1;

                        if (num_y_ > pad_size_ && (min_tot - start_y_) >= pad_size_){
                                
                            int size_shrink=(min_tot - start_y_)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(2, size_shrink);
                        }
                    }
                
                    else if (y + start_y_ == max_IS_y_){ 
                        max_IS_y_ = check_new_border( x, y, 4, false) + start_y_;
                        
                        int max_tot= max(max_IS_y_ , max_vir_y_) + 1;
                        
                        if (num_y_ > pad_size_ && (max_tot - start_y_) < num_y_ - pad_size_){
                                
                            int size_shrink=(num_y_ - (max_tot - start_y_) - 1)/pad_size_;
                            size_shrink*=pad_size_; // remove the closest multiple of pad_size
                            shrink_world(0, size_shrink);
                        }
                        
                    }

                    
                    
                   //~ cout <<  " checked"   << endl;
                }
                
            }                
           //~ cout << toupdate <<endl; 
        }
        //~ cout << toupdate <<endl;
                
//    }
        
        
        
	//~ tot_num_vir_ch = get_tot_vir();    
    //~ if (tot_num_vir_ch!=numvirtot_){
        //~ cout << "IS upd 2 inconsistent numvirtot_! "<< numvirtot_ << " vs " << tot_num_vir_ch << endl;
        //~ exit(1); 
    //~ }
    
    //~ cout <<  " add to vir coord "   << endl;

    for(int j=0; j < vir_coords_.size(); j++){ 
        
        int x =  vir_coords_[j][0] - start_x_;
        int y =  vir_coords_[j][1] - start_y_;
     
        long int num_virs= world_[x][y].get_num_vir();
        
        long int num_IS  = world_[x][y].get_num_IS(); //now after removal
        
        long int old_IS=prec_IS_vec[j]; //before update
        long int new_IS=new_IS_vec[j]; // to be added to the actual one due to viruses
        
        if (new_IS + num_IS > mem_points_*people_number_) {
            cout << "new_IS inconsistent! "<< new_IS << endl;
            exit(1); 
        }
        
        if(new_IS<0){
            cout << "new_IS less than 0! "<< new_IS << endl;
            exit(1); 
        }
       
        world_[x][y].set_num_IS(new_IS + num_IS);

        if (new_IS + num_IS - old_IS!=0){ //  
        
            vector< long int> newcoords {x + start_x_, y + start_y_, new_IS + num_IS - old_IS};
            
            IS_update_coords_.push_back(newcoords);
            
    
    
            
            //cout << "IS update check first loop "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << IS_update_coords_[IS_update_coords_.size()-1][2]  << endl;
            
           if (IS_update_coords_[IS_update_coords_.size()-1][0] > max_x_tot || IS_update_coords_[IS_update_coords_.size()-1][0] < min_x_tot || IS_update_coords_[IS_update_coords_.size()-1][1] > max_y_tot || IS_update_coords_[IS_update_coords_.size()-1][1] < min_y_tot )
            {
                cout << "wrong IS update coords vir loop!  x " << IS_update_coords_[IS_update_coords_.size()-1][0] << " y " << IS_update_coords_[IS_update_coords_.size()-1][1] << " max_x_tot " << max_x_tot << " max_y_tot " << max_y_tot << " min_x_tot " << min_x_tot << " min_y_tot " << min_y_tot << endl;
                cout << "IS update "<< j << " IS_update_coords_.size() " << IS_update_coords_.size()  << " num_upd " << IS_update_coords_[IS_update_coords_.size()-1][2]  << endl;
                exit(1);
            } 
            
        }       
        
        if (num_IS==0 && new_IS>0){// populated a new coord with IS, which had viruses
            
            IS_coords_.push_back(vir_coords_[j]);
            
        
            // check min max x y
            
            if (x + start_x_ < min_IS_x_) min_IS_x_ = x + start_x_;
            if (y + start_y_ < min_IS_y_) min_IS_y_ = y + start_y_;
            if (x + start_x_ > max_IS_x_) max_IS_x_ = x + start_x_;
            if (y + start_y_ > max_IS_y_) max_IS_y_ = y + start_y_;
            
        }
        
    }
    
    
    
	    
 
};




void antigenic_space::growth_virs(){
 
        cout <<  " growth_virs "   << endl;
        cout <<  vir_coords_.size()  << endl;



    
    
    int min_vir_x_chk= numeric_limits<int>::max();  
    int min_vir_y_chk= numeric_limits<int>::max();  
    int max_vir_x_chk= numeric_limits<int>::min();  
    int max_vir_y_chk= numeric_limits<int>::min();  


    
    for(int j=0; j < vir_coords_.size(); j++){
        //cout <<  " cycle"   << endl;
        //
        //cout <<  j  << endl;
        
        int x =  vir_coords_[j][0] - start_x_;
        int y =  vir_coords_[j][1] - start_y_;
     
        long int next_inf;
        long int num_vir_prec= (world_[x][y].get_num_vir());
        
        double vir_gr = ((world_[x][y].get_fitness()) + 1)*num_vir_prec; // linear approximation
        
        //cout <<  num_vir_prec  << endl;
        //cout <<  vir_gr  << endl;
        
        if(vir_gr<0){
            cout << "fitness less than -1! "<< (world_[x][y].get_fitness()) << endl;
            exit(1); 
        }
        
        //POISSON DISTRIBUTION, with average value vir_gr
        
        poisson_distribution<long int> poisson(vir_gr);

        //cout <<  " growth ready"   << endl;
        
        next_inf = poisson(mt);
  

        //cout <<  " drawn"   << endl;
        //
        //cout <<  next_inf  << endl;
        
        //if (next_inf < - num_vir_prec) next_inf = - num_vir_prec;
        //vir_coords_[j].update_virs(next_inf);
        
        if (next_inf <= 0){
            next_inf = 0;
            world_[x][y].set_fitness(0.);
        }
        world_[x][y].set_num_virs(next_inf);
        numvirtot_+=next_inf - num_vir_prec;
        
        if(next_inf==0){// sent to zero the virs, now have to check for IS
            
                       
            if (vir_coords_.size() >= 2)
            {
                swap(vir_coords_[j], *(vir_coords_.end() - 1));
                vir_coords_.pop_back();
                j--;
            }
            else
            {
                vector<vector<int>>().swap(vir_coords_); //clears memory, including pointers

                //infected.clear(); // doesn't clear memory
            }
            
            // check min max x y
            
            if (x + start_x_ == min_vir_x_) min_vir_x_ = check_new_border( x, y, 1, true) + start_x_;
            if (y + start_y_ == min_vir_y_) min_vir_y_ = check_new_border( x, y, 3, true) + start_y_;
            
            if (x + start_x_ == max_vir_x_) max_vir_x_ = check_new_border( x, y, 2, true) + start_x_;
            if (y + start_y_ == max_vir_y_) max_vir_y_ = check_new_border( x, y, 4, true) + start_y_;
            
        }  
        else{
        
            if(x + start_x_ < min_vir_x_chk) min_vir_x_chk = x + start_x_ ;
            if(y + start_y_ < min_vir_y_chk) min_vir_y_chk = y + start_y_ ;
            if(x + start_x_ > max_vir_x_chk) max_vir_x_chk = x + start_x_ ;
            if(y + start_y_ > max_vir_y_chk) max_vir_y_chk = y + start_y_ ;
            
        } 
        
        //cout <<  j  << endl;
        //cout <<  vir_coords_.size()  << endl;

    }
    
    
    if(min_vir_x_chk!=min_vir_x_){ // if sth goes wron in minmax tracking, correct here. It should not
        //~ cout << "correcting min_vir_x_ "<< min_vir_x_ << " to " << min_vir_x_chk << endl;
        min_vir_x_=min_vir_x_chk;
    }
    
    
    if(max_vir_x_chk!=max_vir_x_){ // if sth goes wron in minmax tracking, correct here. It should not
        //~ cout << "correcting max_vir_x_ "<< max_vir_x_ << " to " << max_vir_x_chk << endl;
        max_vir_x_=max_vir_x_chk;
    }
    
    if(min_vir_y_chk!=min_vir_y_){
        //~ cout << "correcting min_vir_y_ "<< min_vir_y_ << " to " << min_vir_y_chk << endl;
        min_vir_y_=min_vir_y_chk;
    }

    if(max_vir_y_chk!=max_vir_y_){ // if sth goes wron in minmax tracking, correct here. It should not
        //~ cout << "correcting max_vir_y_ "<< max_vir_y_ << " to " << max_vir_y_chk << endl;
        max_vir_y_=max_vir_y_chk;
    }
        
	//~ long int tot_num_vir = get_tot_vir();    
    //~ if (tot_num_vir!=numvirtot_){
        //~ cout << "vir gr inconsistent numvirtot_! "<< numvirtot_ << " vs " << tot_num_vir << endl;
        //~ exit(1); 
    //~ }
    
}


void antigenic_space::diffusion_virs(){
    cout <<  " diffusion_virs "   << endl;
    cout << " vir_coords_ " << vir_coords_.size() << endl;
    
    if (jump_size_>1) cout <<  " GAMMA "   << endl;
    else cout <<  " RW "   << endl;



    
    discrete_distribution<int> direction_distr {1,1,1,1};
    //poisson_distribution<long int> poisson(mu_);
    
    vector<double> poiss_probs;
    
    int ev=1;
    double prob= pow(mu_,float(ev))*exp(- mu_ )/tgamma(ev+1);
    poiss_probs.push_back(prob);
    
    while (prob > pow(10., -9.) && poiss_probs.size() < 10 ){
        ev++;
        prob= pow(mu_,float(ev))*exp(- mu_ )/tgamma(ev+1);
        poiss_probs.push_back(prob);
    }
    
    discrete_distribution<int> poiss_positive (poiss_probs.begin(), poiss_probs.end());
    uniform_real_distribution<double> distr_theta(0.0, M_PI); 
    gamma_distribution<double> gamma_distr(20.0,jump_size_/20.0); // if change k=10 remember to change initial conditions!

    
    
    
    double p_mut= 1. - exp(- mu_); // time step 1
    
    int maxnummut=numeric_limits<int>::min();  
    
    vector<long int> num_vir_prec (vir_coords_.size(), 0);
    
    for(int j=0; j < vir_coords_.size(); j++){ // store the previous state 
     
        int x =  vir_coords_[j][0] - start_x_;
        int y =  vir_coords_[j][1] - start_y_;
        
        long int num_virs= world_[x][y].get_num_vir();
        num_vir_prec[j] =num_virs;
        
    }
    
    vector<vector<int>> vir_coords_new;

    for(int j=0; j < vir_coords_.size(); j++){
    
        int x =  vir_coords_[j][0] - start_x_;
        int y =  vir_coords_[j][1] - start_y_;
 
        long int num_virs_p=num_vir_prec[j];
        long int num_virs= world_[x][y].get_num_vir();
        
        binomial_distribution<long int> bin_distribution(num_virs_p,p_mut);
        
        //long int num_mut=bin_distribution(mt);
        long int num_mut=bin_distribution(mt);


       
        world_[x][y].set_num_virs(num_virs - num_mut);
        
        if (num_virs== num_mut){ // Sent virs to 0
            
            world_[x][y].set_fitness(0.);
                              
            if (vir_coords_.size() >= 2)
            {
                swap(vir_coords_[j], *(vir_coords_.end() - 1));
                vir_coords_.pop_back();
                swap(num_vir_prec[j], *(num_vir_prec.end() - 1));
                num_vir_prec.pop_back();
                j--;
            }
            else
            {
                vector<vector<int>>().swap(vir_coords_); //clears memory, including pointers
                vector<long int>().swap(num_vir_prec); //clears memory, including pointers

                //infected.clear(); // doesn't clear memory
            }
           
            // check min max x y
            
            if (x + start_x_ == min_vir_x_) min_vir_x_ = check_new_border( x, y, 1, true) + start_x_;
            if (y + start_y_ == min_vir_y_) min_vir_y_ = check_new_border( x, y, 3, true) + start_y_;
            
            if (x + start_x_ == max_vir_x_) max_vir_x_ = check_new_border( x, y, 2, true) + start_x_;
            if (y + start_y_ == max_vir_y_) max_vir_y_ = check_new_border( x, y, 4, true) + start_y_;
             
        }
        
        int x_ref = x + start_x_;
        int y_ref = y + start_y_;

        for(int k=0; k < num_mut; k++){// super slow with many virs, if becomes bottleneck I can sample dir from GSL multinomial
                
                
            int mult_mut= 1;
            
            if(poiss_probs.size() > 1 ) mult_mut= poiss_positive(mt) + 1; // +1 because draws vect index from 0
            
            if (mult_mut <=0 || mult_mut>10){
                cout << "wrong number of double mutations! "<< mult_mut << " vs " << poiss_probs.size() << endl;
                // Generate an interrupt, for debugger
                //~ std::raise(SIGINT);
        
                exit(1); 
            }
            
            if(mult_mut > maxnummut) maxnummut = mult_mut ;
           
            int x_nn= x_ref - start_x_; // starting point on the grid, correct even if start changes
            int y_nn= y_ref - start_y_;
            
            
            // gamma distributed jumps
            
            if (jump_size_>1){
    
                double theta=distr_theta(mt);
                theta*=2.0;
    
                double dist;
                
                int x_displ; // on the grid directly
                int y_displ;
                
                do{
        
                    dist=gamma_distr(mt)/latt_width_;
                    
                    x_displ= round(cos(theta)*dist);
                    y_displ= round(sin(theta)*dist);
                    
                } while(x_displ==0 && y_displ==0); // choose jump distr so that probability of this happening is low
                
                
                while(mult_mut>0){
                    x_nn += x_displ;
                    y_nn += y_displ;
                        
                    if (y_nn > num_y_ -1) y_nn= num_y_ -1;
                    if (y_nn < 0)  y_nn=0;
                    if (x_nn > num_x_ -1) x_nn= num_x_ - 1;
                    if (x_nn < 0)  x_nn=0;
                    
                    theta=distr_theta(mt);
                    theta*=2.0;
                    dist=gamma_distr(mt)/latt_width_;
                    x_displ= round(cos(theta)*dist);
                    y_displ= round(sin(theta)*dist);  
                           
                    mult_mut--;
                }
                
            }
            // random walk, not used in paper
            else{
                int dir=direction_distr(mt); // from 0 to 3
                
                while(mult_mut>0){
                     
                    if (dir==0) y_nn++;
                    if (dir==2) y_nn--;
                    if (dir==1) x_nn++;
                    if (dir==3) x_nn--;
                    dir=direction_distr(mt);
                    mult_mut--;
                    
                    if (y_nn > num_y_ -1  || y_nn < 0) {
                        cout << "viral coord already at y border! "<< y_nn << " founder " << y_ref - start_y_ << " num_y_ " << num_y_ << endl;
                        // Generate an interrupt
                        //~ std::raise(SIGINT);
                
                        exit(1); 
                    }
                    
                    if (x_nn > num_x_ -1  || x_nn < 0) {
                        cout << "viral coord already at x border! "<< x_nn << " founder " << x_ref - start_x_ << " num_x_ " << num_x_ << endl;
                        // Generate an interrupt
                        //~ std::raise(SIGINT);
                
                        exit(1); 
                    }

                }
            }
        
            double displ_fin_x = (x_nn - (x_ref - start_x_)) * latt_width_;
            double displ_fin_y = (y_nn - (y_ref - start_y_)) * latt_width_;
            
            displ_x_    +=displ_fin_x;
            displ_y_    +=displ_fin_y;
            displ_x_sq_ +=pow(displ_fin_x, 2.);
            displ_y_sq_ +=pow(displ_fin_y, 2.);
            
            displ_tot_ += sqrt(pow(displ_fin_x, 2.) + pow(displ_fin_y, 2.));
            
            count_displ_++;
            
            //antigenic_coordinate nn = vir_coords_[j]. get_nn(dir);
            
            long int num_virs_nn= world_[x_nn][y_nn].get_num_vir();
    
//           
//            cout <<  " diffusion in "   << endl;
//            cout << x_nn << " "  << y_nn << endl;
            
            world_[x_nn][y_nn].set_num_virs(num_virs_nn +1);
            
            if (num_virs_nn==0){ // created new virus! expand grid at need
                vector<int> newcoords {x_nn + start_x_, y_nn + start_y_};
                
                vir_coords_new.push_back(newcoords);
                
                // check min max x y
                
                if (x_nn + start_x_ < min_vir_x_) min_vir_x_ = x_nn + start_x_;
                if (y_nn + start_y_ < min_vir_y_) min_vir_y_ = y_nn + start_y_;
                if (x_nn + start_x_ > max_vir_x_) max_vir_x_ = x_nn + start_x_;
                if (y_nn + start_y_ > max_vir_y_) max_vir_y_ = y_nn + start_y_;
                
            }
        
            if (y_nn >= num_y_ -40) expand_world(0);
            if (y_nn < 40) expand_world(2);
            
            
            if (x_nn >= num_x_ -40) expand_world(1);
            if (x_nn < 40) expand_world(3);
            
                
        } 
     
    }    
    // append new coord to virs
    
    //~ cout << "maxnummut "<< maxnummut << endl;
    //~ cout << "vir_coords_new "<< vir_coords_new.size() << " vir_coords_ " << vir_coords_.size() << endl;
    
    vir_coords_.insert( vir_coords_.end(), vir_coords_new.begin(), vir_coords_new.end() );

    
	long int tot_num_vir = get_tot_vir();    
    if (tot_num_vir!=numvirtot_){
        cout << "diff virs inconsistent numvirtot_! "<< numvirtot_ << " vs " << tot_num_vir << endl;
        cout << "vir_coords_new "<< vir_coords_new.size() << " vir_coords_ " << vir_coords_.size() << endl;
        // Generate an interrupt
        //~ std::raise(SIGINT);

        exit(1); 
    }
}


void antigenic_space::expand_world(int dir){

    cout <<  " expand_world "   << endl;
    cout << dir << endl;
    
    
    antigenic_coordinate empty;
    
    int padsize=pad_size_;
    
    if (dir==1){//right
    
        //deque<antigenic_coordinate> empty_col(world_[0].size(), empty);
        vector<antigenic_coordinate> empty_col(world_[0].size(), empty);
        
        for(int i=0; i < padsize; i++){
            world_.push_back(empty_col);
        }
        
        set_num_x(num_x_ + padsize);
       
        
    }
    else if (dir==3){//left
    
        //deque<antigenic_coordinate> empty_col(world_[0].size(), empty);
        vector<antigenic_coordinate> empty_col(world_[0].size(), empty);
        
        for(int i=0; i < padsize; i++){
            // for deque
            //world_.push_front(empty_col);
            //for vector
            world_.push_back(empty_col);
            rotate(world_.rbegin(), world_.rbegin() + 1, world_.rend());
        }
        //.insert(world_[i].begin(), empty_col_slice.begin(), empty_col_slice.end());  
        
        set_num_x(num_x_ + padsize);
        set_start_x(start_x_ - padsize);
    }
    else if (dir==0){//up
    
        //deque<antigenic_coordinate> empty_col_slice(padsize, empty);
        vector<antigenic_coordinate> empty_col_slice(padsize, empty);
        
        for(int i=0; i < world_.size(); i++){
            world_[i].insert(world_[i].end(), empty_col_slice.begin(), empty_col_slice.end());  
        }
        
        set_num_y(num_y_ + padsize);
    }
    else if (dir==2){//down
    
        //deque<antigenic_coordinate> empty_col_slice(padsize, empty);
        vector<antigenic_coordinate> empty_col_slice(padsize, empty);
        
        for(int i=0; i < world_.size(); i++){
            world_[i].insert(world_[i].begin(), empty_col_slice.begin(), empty_col_slice.end());  
        }
        
        set_num_y(num_y_ + padsize);
        set_start_y(start_y_ - padsize);
    }
    else{
        cout<<"Error, direction not implemented!"<<endl;
        exit(1);
    }    
    
//	long int tot_num_vir = get_tot_vir();    
//    if (tot_num_vir!=numvirtot_){
//        cout << "exp world inconsistent numvirtot_! "<< numvirtot_ << " vs " << tot_num_vir << endl;
//        exit(1); 
//    }
    
}

int antigenic_space::check_new_border( int x, int y, int bord, bool vir){    
 

    if (bord==1){ // check right for new xmin
        
        int xc=x; 
        while(xc < world_.size()){
            int yc=0; 
            while(yc < world_[xc].size()){
                
                //if (yc!=y || xc !=x){
                
                long int num_tot_c  = 0;
                if (vir) num_tot_c  = world_[xc][yc].get_num_vir();
                else num_tot_c  = world_[xc][yc].get_num_IS();
                
                if (num_tot_c>0) break;
                    
                yc++;
            }
            
            if (yc < world_[xc].size()) break;
            xc++;
        }
        if (xc == world_.size()){ // he is the "last", now calculate how much to shrink
            
            
            if(numvirtot_>0){
                cout<<"Error, new minx not found!"<<endl;
                 cout << x << " "  << y << endl;
                cout << num_x_ << " "  << num_y_ << endl;
                cout << min_IS_x_ << " "  << min_vir_x_ << " "  << min_IS_x_ - start_x_ << " "  << min_IS_y_ - start_y_ << " "  << start_x_ << " "  <<  start_y_ << endl;
                cout << world_.size() << " "  << world_[0].size() << " "  << world_[x].size() << endl;
                // Generate an interrupt
                //~ std::raise(SIGINT);

               exit(1);
            }
            
        }
        
        return xc;
        
        
    }
    else if (bord==2){ // check left for new xmax
    

        int xc=x; 
        while(xc >=0){
            int yc=0; 
            while(yc < world_[xc].size()){
                
                //if (yc!=y || xc !=x){

                long int num_tot_c  = 0;
                if (vir) num_tot_c  = world_[xc][yc].get_num_vir();
                else num_tot_c  = world_[xc][yc].get_num_IS();

                if (num_tot_c>0) break;

                yc++;
            }
            
            if (yc < world_[xc].size()) break;
            xc--;
        }
        if (xc == -1){ // he is the "last", now calculate how much to shrink
            
            if(numvirtot_>0){
                cout<<"Error, new maxx not found!"<<endl;
                 cout << x << " "  << y << endl;
                cout << num_x_ << " "  << num_y_ << endl;
                cout << max_IS_x_ << " "  << max_vir_x_ << " "  << max_IS_x_ - start_x_ << " "  << max_IS_y_ - start_y_ << " "  << start_x_ << " "  <<  start_y_ << endl;
                cout << world_.size() << " "  << world_[0].size() << " "  << world_[x].size() << endl;
                          // Generate an interrupt
                //~ std::raise(SIGINT);

                exit(1);
            }
            
            
        }
        
        return xc;

    }
    else if (bord==3){ // check top for new ymin
        

        int yc=y; 
        //while(yc < world_[x].size()){
        while(yc < world_[0].size()){
            int xc=0; 
            while(xc < world_.size()){
                
                //if (yc!=y || xc !=x){

                long int num_tot_c  = 0;
                if (vir) num_tot_c  = world_[xc][yc].get_num_vir();
                else num_tot_c  = world_[xc][yc].get_num_IS();
                
                if (num_tot_c>0) break;
                    
                xc++;
            }
            
            if (xc < world_.size()) break;
            yc++;
        }
        //if (yc == world_[x].size()){ // he is the "last", now calculate how much to shrink
        if (yc == world_[0].size()){ // he is the "last", now calculate how much to shrink
            
            if(numvirtot_>0){
                cout<<"Error, new miny not found!"<<endl;
                cout << x << " "  << y << endl;
                cout << num_x_ << " "  << num_y_ << endl;
                cout << min_IS_y_ << " "  << min_vir_y_ << " "  << min_IS_x_ - start_x_ << " "  << min_IS_y_ - start_y_ << " "  << start_x_ << " "  <<  start_y_ << endl;
                cout << world_.size() << " "  << world_[0].size() << " "  << world_[x].size() << endl;
                // Generate an interrupt
                //~ std::raise(SIGINT);

                exit(1);
            }
            
        }
        return yc;
    }
    else if (bord==4){ // check down for new ymax

        int yc=y; 
        while(yc >=0){
            int xc=0; 
            while(xc < world_.size()){
                
                //if (yc!=y || xc !=x){
                
                long int num_tot_c  = 0;
                if (vir) num_tot_c  = world_[xc][yc].get_num_vir();
                else num_tot_c  = world_[xc][yc].get_num_IS();
                
                if (num_tot_c>0) break;
                    
                xc++;
            }
            
            if (xc < world_.size()) break;
            yc--;
        }
        if (yc == -1){ // he is the "last", now calculate how much to shrink
            
            if(numvirtot_>0){
                cout<<"Error, new msxy not found!"<<endl;
                 cout << x << " "  << y << endl;
                cout << num_x_ << " "  << num_y_ << endl;
                cout << max_IS_y_ << " "  << max_vir_y_ << " "  << max_IS_x_ - start_x_ << " "  << max_IS_y_ - start_y_ << " "  << start_x_ << " "  <<  start_y_ << endl;
                cout << world_.size() << " "  << world_[0].size() << " "  << world_[x].size() << endl;
                // Generate an interrupt, send same signal to program as control C
                //~ std::raise(SIGINT);
               exit(1);
            }
            
        } 
        return yc;                      
    }
        
    else{
        cout<<"Error, border not implemented!"<<endl;
        exit(1);
    }
    
    
//	long int tot_num_vir = get_tot_vir();    
//    if (tot_num_vir!=numvirtot_){
//        cout << "check bord inconsistent numvirtot_! "<< numvirtot_ << " vs " << tot_num_vir << endl;
//        exit(1); 
//    }    
// 
    
}


void antigenic_space::checkpad_and_shrink(int x, int y){ 
    
    int padsize=pad_size_;
    
    if (num_x_ > padsize){

        if (x >= num_x_ - padsize){ //right
            
            int xc=x; 
            while(xc < world_.size()){
                int yc=0; 
                while(yc < world_[xc].size()){
                    
                    //if (yc!=y || xc !=x){
                    
                    long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                    
                    if (num_tot_c>0) break;
                        
                    yc++;
                }
                
                if (yc < world_[xc].size()) break;
                xc++;
            }
            if (xc == world_.size()){ // he is the "last", now calculate how much to shrink
                
                int xc=x-1; 
                while(xc >=0){
                    int yc=0; 
                    while(yc < world_[xc].size()){
                        
                        long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                        
                        if (num_tot_c>0) break;
                            
                        yc++;
                    }
                    
                    if (yc < world_[xc].size()) break;
                    xc--;
                }                       
                
                int size_shrink=(num_x_ - xc - 1)/padsize;
                size_shrink*=padsize; // remove the closest multiple of padsize
                shrink_world(1, size_shrink);
                
            }
            
            
        }
        if (x < padsize) { // left
        

            int xc=x; 
            while(xc >=0){
                int yc=0; 
                while(yc < world_[xc].size()){
                    
                    //if (yc!=y || xc !=x){
                    
                    long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                    
                    if (num_tot_c>0) break;
                        
                    yc++;
                }
                
                if (yc < world_[xc].size()) break;
                xc--;
            }
            if (xc == -1){ // he is the "last", now calculate how much to shrink
                
                int xc=x+1; 
                while(xc< world_.size()){
                    int yc=0; 
                    while(yc < world_[xc].size()){
                        
                        long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                        
                        if (num_tot_c>0) break;
                            
                        yc++;
                    }
                    
                    if (yc < world_[xc].size()) break;
                    xc++;
                }                       
                
                int size_shrink=(xc)/padsize;
                size_shrink*=padsize; // remove the closest multiple of padsize
                shrink_world(3, size_shrink);
                
            }
        }
        
    }
    if (num_y_ > padsize){
            
        if (y >= num_y_ - padsize) { //top
            

            int yc=y; 
            while(yc < world_[x].size()){
                int xc=0; 
                while(xc < world_.size()){
                    
                    //if (yc!=y || xc !=x){
                    
                    long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                    
                    if (num_tot_c>0) break;
                        
                    xc++;
                }
                
                if (xc < world_.size()) break;
                yc++;
            }
            if (yc == world_[x].size()){ // he is the "last", now calculate how much to shrink
                
                int yc=y-1; 
                while(yc >=0){
                    int xc=0; 
                    while(xc < world_.size()){
                        
                        long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                        
                        if (num_tot_c>0) break;
                            
                        xc++;
                    }
                    
                    if (xc < world_.size()) break;
                    yc--;
                }                       
                
                int size_shrink=(num_y_ - yc - 1)/padsize;
                size_shrink*=padsize; // remove the closest multiple of padsize
                shrink_world(0, size_shrink);
                
            }
        }
        if (y < padsize) { // down

            int yc=y; 
            while(yc >=0){
                int xc=0; 
                while(xc < world_.size()){
                    
                    //if (yc!=y || xc !=x){
                    
                    long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                    
                    if (num_tot_c>0) break;
                        
                    xc++;
                }
                
                if (xc < world_.size()) break;
                yc--;
            }
            if (yc == 0){ // he is the "last", now calculate how much to shrink
                
                int yc=y+1; 
                while(yc < world_[x].size()){
                    int xc=0; 
                    while(xc < world_.size()){
                        
                        long int num_tot_c  = world_[xc][yc].get_num_vir() + world_[xc][yc].get_num_IS();
                        
                        if (num_tot_c>0) break;
                        xc++;
                    }
                    
                    if (xc < world_.size()) break;
                    yc++;
                }                       
                
                int size_shrink= (yc)/padsize;
                size_shrink*=padsize; // remove the closest multiple of padsize
                shrink_world(2, size_shrink);
                
            }                       
        }
    }
    
    
}

void antigenic_space::shrink_world(int dir, int shrinksize){ // it's important that the check if I'm not cutting out stuff is done before

    cout <<  " shrink_world "   << endl;
    cout << dir << " "  << shrinksize << endl;
    cout << world_.size() << " "  << world_[0].size() << endl;
    
    
    if (dir==1){ //right
        
        world_.erase(world_.end() - shrinksize ,world_.end());
        world_.shrink_to_fit();
        
        set_num_x(num_x_ - shrinksize );
    }
    else if (dir==3){//left
        
        world_.erase(world_.begin(), world_.begin() + shrinksize);
        world_.shrink_to_fit();

        
        set_num_x(num_x_ - shrinksize);
        set_start_x(start_x_ + shrinksize);
    }
    else if (dir==0){//up
    
        
        for(int i=0; i < world_.size(); i++){
            world_[i].erase(world_[i].end() - shrinksize ,world_[i].end()); 
            world_[i].shrink_to_fit();

        }
        
        set_num_y(num_y_ - shrinksize);
    }
    else if (dir==2){//down
        
        for(int i=0; i < world_.size(); i++){
            world_[i].erase(world_[i].begin(), world_[i].begin() + shrinksize); 
            world_[i].shrink_to_fit();
        }
        
        set_num_y(num_y_ - shrinksize);
        set_start_y(start_y_ + shrinksize);
    }
    else{
        cout<<"Error, direction not implemented!"<<endl;
        exit(1);
    }    
     
     //~ cout <<  " shrinked "   << endl;
   
    //~ cout << world_.size() << " "  << world_[0].size() << endl;

}

//deque<antigenic_coordinate>& antigenic_space::operator[] (const int i)
vector<antigenic_coordinate>& antigenic_space::operator[] (const int i)
{
    if( world_.size()<=i ){
        cout<<"ERROR TRYING TO ACCESS NON EXISTING COORDINATES"<<endl;
    }
    else return world_[i];
};

//const deque<antigenic_coordinate>& antigenic_space::operator[] (const int i) const
const vector<antigenic_coordinate>& antigenic_space::operator[] (const int i) const
{
    if( world_.size()<=i ){
        cout<<"ERROR TRYING TO ACCESS NON EXISTING COORDINATES "<<endl;
    }
    else return world_[i];
};


//perform the infection cycle  and check for explosions and extinctions, saving appropriate files upon these events for downstream analysis. Decide if update the fitness using the exact convolution on updated coordinates, based on accumulated approximation error

void antigenic_space::infection_cycle(int & count_full_upd, bool benchmark_update_fitness) {
    bool approx_fullIS=false;
    
    
	update_IS(); //
	growth_virs(); //
    if (numvirtot_>people_number_/2.){
        cout << "too many viruses, aborting simulation! "<< numvirtot_ << " vs " << people_number_ << endl; 
        
  	    string expl_f ("expl_file.txt");
  	    
  	    fstream expl_file;
  	    expl_file.open(out_dir_ + expl_f, fstream::out); 
  	    
  	    expl_file<<"true"<<endl;
  	    expl_file.close(); 
  	    
        exit(1);
    }
    if (numvirtot_<=1){
        cout << "viruses extincts, aborting simulation! "<< numvirtot_ << endl; 
        
  	    string expl_f ("extinct_file.txt");
  	    
  	    fstream expl_file;
  	    expl_file.open(out_dir_ + expl_f, fstream::out); 
  	    
  	    expl_file<<"true"<<endl;
  	    expl_file.close(); 
  	    
        exit(1);
    }
	diffusion_virs(); //
    
    double fitn_avg_o_std=0;
    
    double maxfiterr_o_std=0;
    
	//update_fitness_fullconv(); //
	//update_fitness_fullconv_update(); //
	//update_fitness_IS_update_NFFT(true); //
    //update_fitness_IS_update_combine_farfield_NFFT(false, true, true); // second bool decides if theta bins, third sets if all is coords approximated by farfield theta (if no theta automathically set to false)
    if(benchmark_update_fitness){
         fitn_avg_o_std=benchmark_update_fitness_fullconv_update_NFFT_farfield();
         //if(fitn_avg_o_std> 0.1) count_full_upd =100; // benchmark updates fullconv, no need to check afterwards
         
         cout << "UPDATE FITNESS CHECK "<< fitn_avg_o_std << " " << count_full_upd << endl; 
         
    } //
    else if(count_full_upd > 0){ 
        update_fitness_fullconv_update(true);
        count_full_upd--;
        cout << "UPDATE FITNESS CHECK FULLUPD "<< count_full_upd << endl; 
    }
    else{
        update_fitness_IS_update_switch(true, maxfiterr_o_std);
        
         if(maxfiterr_o_std> 0.1) count_full_upd =10;
         
         cout << "UPDATE FITNESS CHECK "<< maxfiterr_o_std << " " << count_full_upd << endl; 

    }
    //else update_fitness_IS_update_combine_farfield_NFFT(true, true, true, true, approx_fullIS);
    //else update_fitness_fullconv_update(true);
    //else update_fitness_IS_update_combine_farfield_fullconv(true, true, true);

};


// creates the world, putting IS and virs in good initial conditions as explained in Methods
antigenic_space create_world(int DIM, double recog_width_real, double latt_width, double jump_size, long int people_number, int mem_points, string kernel, double mu, string out_dir, double F_0){ // THIS ASSUME S LATT WIDTH 1. IF DIFFERENT PASS recog_width/latt_width 
    
    bool anysot_IC=false;
 
    double recog_width = recog_width_real / latt_width  ; //;/jump_size
    
    double diff_const = pow(jump_size, 2.)*mu; // discrete variance  /(DIM  * 2.)
 
 
    int ret_inf= 1000; // timescale for population memory
    
    //ret_inf= int(pow(F_0,2/3.) * pow(recog_width,4/3.)/(pow(mu,2/3.)* pow(24*log(pow(mu,1/3.) * pow(F_0/recog_width,2/3.) * people_number/ret_inf ),1/3.)) ); // timescale for population memory
    ret_inf= int(pow(recog_width,4/3.)/(pow(diff_const,2/3.)* pow(mem_points*(exp(F_0/mem_points) -1.),4/3.) * pow(24*log(pow(diff_const,1/3.) * pow(mem_points*(exp(F_0/mem_points) -1.)/recog_width,2/3.) * people_number/ret_inf ),1/3.)) ) ; // timescale for population memory  /(jump_size *jump_size ) 
    
    long int num_vir_IC = people_number/ret_inf;//number of viruses in total
    
    double radius_vir_IC= pow(diff_const, 1/3.)*pow(mem_points*(exp(F_0/mem_points) -1.)/recog_width, -1/3.)*pow(24*log(num_vir_IC*pow(diff_const, 1/3.)*pow(mem_points*(exp(F_0/mem_points) -1.)/recog_width, 2/3.)), 1/6.); //jump_size * 
    
    
    
    
    int fittest_x= int(max(pow(diff_const, 1/3.)*pow(mem_points*(exp(F_0/mem_points) -1.)/recog_width, -1/3.)*pow(24*log(num_vir_IC*pow(diff_const, 1/3.)*pow(mem_points*(exp(F_0/mem_points) -1.)/recog_width, 2/3.)), 2/3.)/4., 2.));// jump_size *
    
    long int num_vir_fittest_line= (num_vir_IC/1000)*fittest_x;
    //int num_vir_fittest_line= fittest_x;


    long int num_vir_IC_per_pos= num_vir_IC - num_vir_fittest_line; //number of viruses in the gaussian
    
    
    
    double y_IS_IC=ceil(max(2*radius_vir_IC, 1.));
    
    int num_y_IS_IC= 2* y_IS_IC;
    
    

    
    double tau=mem_points*ret_inf;
    double s=(mem_points*(exp(F_0/mem_points) -1.))/recog_width;
    
    double v_tau= mem_points/s; // jump_size * 
 
 

    int memtrace=  v_tau * log(mem_points*people_number/(v_tau* num_y_IS_IC)); // this is the minimum x that IS will cover, starting from x=0
     
    int pad_size = 50;
    
    int num_x_IC = memtrace + max(1.*y_IS_IC, 1.*fittest_x) + 2*pad_size;//grid dimensions, number of entries
    //int num_y_IC = 6*radius_vir_IC + 100;
    int num_y_IC = num_y_IS_IC + 2*pad_size;
    
    int start_x_IC = -memtrace - pad_size; // coordinates corresponding to low left corner of current grid, in grid number units, leave 50 on each side IC 
    //int start_y_IC = -2*radius_vir_IC - 50;
    //int start_y_IC = -3*radius_vir_IC - 50;
    int start_y_IC = -y_IS_IC - pad_size;
    
    if(anysot_IC){
        
        num_x_IC = memtrace + max(1.*y_IS_IC, 1.*fittest_x) + 2*pad_size;
        num_y_IC = memtrace + max(1.*y_IS_IC, 1.*fittest_x) + 2*pad_size;
        start_x_IC = -memtrace - pad_size;
        start_y_IC = -memtrace - pad_size;
        
        
    }
    
    
    
    
    string IC_trav_wave_f ("IC_trav_wave.txt");

    string outfile= out_dir+ IC_trav_wave_f;

    fstream outstream;
    
    outstream.open(outfile, fstream::out); 
    outstream<<"# 1 Num vir"<<setw(30)<<" 2 tau "<<setw(30)<<" 3 vel"<<setw(30)<<"4 s"<<setw(30)<<"5 sigma space" <<setw(30)<<"6 u_c space" <<setw(30)<<"7 R pers_l" <<endl;
     
    
    outstream << fixed  << setprecision(0)<<setw(19) << num_vir_IC <<setprecision(5)<<setw(19) << tau <<setprecision(5)<<setw(19) << mem_points/(s*tau) <<setprecision(5)<<setw(19) <<  s <<setprecision(5)<<setw(19) <<  radius_vir_IC <<setprecision(5)<<setw(19) <<  fittest_x <<setprecision(5)<<setw(19) <<  mem_points*recog_width/(s*recog_width + mem_points) <<endl;

       
    
    
    cout << "IC PARAMETERS, num_vir_IC "<< num_vir_IC << " ret_inf " << ret_inf << " memtrace " << memtrace << " radius_vir_IC " << radius_vir_IC << " fittest_x " << fittest_x << " fittest_x " << fittest_x  << " max gauss " << int(num_vir_IC_per_pos /(2*M_PI*pow(radius_vir_IC, 2.))) << endl;
        
    long int numvireff=0;
    long int extra_IS=0;
    long int num_IS_eff=0;
 
 
    int min_vir_x= numeric_limits<int>::max();        
    int min_vir_y= numeric_limits<int>::max();        
    int max_vir_x= numeric_limits<int>::min();        
    int max_vir_y= numeric_limits<int>::min();        
    int min_IS_x = numeric_limits<int>::max();        
    int min_IS_y = numeric_limits<int>::max();        
    int max_IS_x = numeric_limits<int>::min();        
    int max_IS_y = numeric_limits<int>::min();        
            
        
        
 
 
   
    //deque<deque<antigenic_coordinate>> w;
    vector<vector<antigenic_coordinate>> w;
    vector<vector<int>> vir_coords;
    vector<vector<int>> IS_coords;
    vector<vector<long int>> IS_coords_update;
    
    double rot=-M_PI*0./8;
    //num_y_IS_IC-=2;
    y_IS_IC*=cos(rot) - sin(rot) ;
    //double  dens_fact=2*cos(rot);
    //double  dens_fact=num_y_IS_IC/(num_y_IS_IC +  ( num_y_IS_IC - 1.) );
    double  dens_fact=1.;
 
    
    // grid initialization
    
    antigenic_coordinate empty;
    
    //deque<antigenic_coordinate> empty_col(num_y_IC, empty);
    vector<antigenic_coordinate> empty_col(num_y_IC, empty);

    for(int i=0; i < num_x_IC; i++){
        w.push_back(empty_col);
    }    
    //for(int i=0; i < w.size(); i++){ 
    for(int i=w.size() -1; i >=0; i--){ 
        for(int j= w[i].size() -1 ; j >=0; j--){ 
            double dist_centr= sqrt( pow( i + start_x_IC, 2.) + pow( j + start_y_IC, 2.));
            double dist_centr_L1= sqrt( pow( i + start_x_IC, 2.)) + 10* sqrt(pow( j + start_y_IC, 2.));
            
            
            double erf_x=0;
            if (i>0) erf_x=((erf((i + start_x_IC)/(sqrt(2.)*radius_vir_IC)) - erf((i-1 + start_x_IC)/(sqrt(2.)*radius_vir_IC)) )/2) ;
            //if (i + start_x_IC == 0 || i + start_x_IC == 1)
            double erf_y=0;
            if (j>0) erf_y=((erf((j + start_y_IC)/(sqrt(2.)*radius_vir_IC)) - erf((j-1 + start_y_IC)/(sqrt(2.)*radius_vir_IC)) )/2) ;
            
            long int numvirsite=(long)(num_vir_IC_per_pos * erf_x * erf_y);
            
            
            double coord_1=i+ start_x_IC ;
            double coord_2=j+ start_y_IC ;
            
            if(anysot_IC) {
                coord_1= cos(rot)*(i+ start_x_IC) - sin(rot)*(j+ start_y_IC);
                coord_2= sin(rot)*(i+ start_x_IC) + cos(rot)*(j+ start_y_IC);
            }
            
            if (coord_1 <= fittest_x && coord_1 > 0 && coord_2  ==0){ 
                //numvirsite+=1;
                numvirsite+=(num_vir_IC/1000);
            }
            
            if (numvirsite>0){ 
                w[i][j].set_num_virs(numvirsite);
                vector<int> coord {i + start_x_IC, j + start_y_IC};
                vir_coords.push_back(coord);   
                numvireff+= numvirsite; 
                
                if(i + start_x_IC < min_vir_x) min_vir_x = i + start_x_IC ;
                if(j + start_y_IC < min_vir_y) min_vir_y = j + start_y_IC ;
                if(i + start_x_IC > max_vir_x) max_vir_x = i + start_x_IC ;
                if(j + start_y_IC > max_vir_y) max_vir_y = j + start_y_IC ;
                
            }  
                  
                  
            
            double bin_wrt0= 1 - (coord_1) ;
            
            long int numISsite=0;
            //long int numISsite_tmp=dens_fact* people_number*mem_points* exp(- bin_wrt0/v_tau)*(exp((1.)/v_tau) -1.)/num_y_IS_IC;
            long int numISsite_tmp=dens_fact* people_number*mem_points* exp(- bin_wrt0/v_tau)*(exp((cos(rot))/v_tau) -1.)/num_y_IS_IC;
            
            if (coord_1 > -memtrace && coord_1 <=0 && coord_2 > - y_IS_IC && coord_2 <= y_IS_IC && num_IS_eff + numISsite_tmp<= mem_points*people_number){ //  && (i+ start_x_IC + j + start_y_IC)%2 == 0  && fmod(coord_1/(cos(rot) + sin(rot)), 1) == 0
                numISsite += numISsite_tmp;
                
                //prob_x= mem_points* exp(coord_1/v_tau)/v_tau;
                //prob_y= 1/(2*y_IS_IC);
                
            
            
            }
            
            

            if (numISsite>0 ){ 
                w[i][j].set_num_IS(numISsite);
                vector<int> coord {i + start_x_IC, j + start_y_IC};
                IS_coords.push_back(coord);
                num_IS_eff += numISsite;
                

                if(i + start_x_IC < min_IS_x) min_IS_x = i + start_x_IC ;
                if(j + start_y_IC < min_IS_y) min_IS_y = j + start_y_IC ;
                if(i + start_x_IC > max_IS_x) max_IS_x = i + start_x_IC ;
                if(j + start_y_IC > max_IS_y) max_IS_y = j + start_y_IC ;
                
            }
            
        }
    }
    
    extra_IS = mem_points*people_number - num_IS_eff;
    
    cout << "IC PARAMETERS, numvireff "<< numvireff << endl;
    cout << "IC PARAMETERS, num_IS_eff "<< num_IS_eff << " extra_IS "<< extra_IS << endl;
    
    if(extra_IS >0 && anysot_IC){
        
        int extra_IS_tot = extra_IS;
        
        for(int j=0; j < IS_coords.size(); j++){
            
            //cout <<endl;
            //cout << j <<endl;

            int x =  IS_coords[j][0] - start_x_IC;
            int y =  IS_coords[j][1] - start_y_IC;
            
            long int num_IS  = w[x][y].get_num_IS();
            
            //double prob_upd=num_IS/float(mem_points*people_number);
            double prob_upd=num_IS/float(num_IS_eff);
            
            binomial_distribution<long int> bin_distribution(extra_IS_tot,prob_upd);
            
            long int num_upd=bin_distribution(mt);


            
            if (num_upd >extra_IS) num_upd = extra_IS;
            
            w[x][y].set_num_IS(num_IS + num_upd);
            
            extra_IS -= num_upd;
            
            if(extra_IS==0) break;
            
        }
        
    }
    
    cout << "IC PARAMETERS, num_IS_eff "<< num_IS_eff << " extra_IS "<< extra_IS << endl;

    if(extra_IS >0 ){
        int extra_dump_x = - start_x_IC - memtrace ;
        int extra_dump_y = - start_y_IC ;
        
        long int num_prec_extra_coord= w[extra_dump_x][extra_dump_y].get_num_IS();
        
        w[extra_dump_x ][extra_dump_y].set_num_IS(num_prec_extra_coord + extra_IS);
        
        if (num_prec_extra_coord == 0){
        
            vector<int> coord {extra_dump_x  + start_x_IC, extra_dump_y + start_y_IC};
            IS_coords.push_back(coord);
    
            if(extra_dump_x + start_x_IC < min_IS_x) min_IS_x = extra_dump_x + start_x_IC ;
            if(extra_dump_y + start_y_IC < min_IS_y) min_IS_y = extra_dump_y + start_y_IC ;
            if(extra_dump_x + start_x_IC > max_IS_x) max_IS_x = extra_dump_x + start_x_IC ;
            if(extra_dump_y + start_y_IC > max_IS_y) max_IS_y = extra_dump_y + start_y_IC ;
            
        }
    }
    
	antigenic_space world(vir_coords,  IS_coords, IS_coords_update, DIM, recog_width_real,  latt_width,  jump_size,  people_number,  mem_points,  kernel,  mu,  out_dir,   num_x_IC,   num_y_IC,   start_x_IC,   start_y_IC, w, F_0, pad_size, min_vir_x, min_vir_y, max_vir_x, max_vir_y, min_IS_x, min_IS_y, max_IS_x, max_IS_y);// 
    
    return world;

}; 



double fitness_fct(double conv, int M, double F_0){ // fitness gicen the convolution
    
    //double F_0=1.;
    //double F_0=0.1;
    //
    //double B= F_0*(F_0 + 1);
    //
    //return F_0 - B*conv;
    
    //double B= F_0*(F_0 + 1);
    
    return F_0 + M*log(1 - conv/M);
    
     // double fitn= 1.*exp(-10.*conv) -1 ;

};



double fitness_inv_fct(double fitn, int M, double F_0){ // return the convolution from fitness
    
    //double F_0=1.;
    //double F_0=0.1;
    //
    //double B= F_0*(F_0 + 1);
    //
    //return F_0 - B*conv;
    
    //double B= F_0*(F_0 + 1);
    
    return  M*(1 -exp((fitn - F_0)/M));
    
     // double fitn= 1.*exp(-10.*conv) -1 ;

};




vector<double> intersections_rect(double x1, double y1, double h, double w) // given a vector, and h and w of a rectangle , finds two vectors (the origin is the rect center) on the first vector direction, intersecting perimeter.
{
    //rescale vector
    
    double x1_resc=x1/sqrt( pow( x1, 2.) + pow( y1, 2.));
    double y1_resc=y1/sqrt( pow( x1, 2.) + pow( y1, 2.));
    
    
    double lambda_x;
    double lambda_y;
    
    if (x1_resc==0) lambda_x = numeric_limits<double>::infinity();
    else lambda_x = w/(2*fabs(x1_resc));
    
    if (y1_resc==0) lambda_y = numeric_limits<double>::infinity();
    else lambda_y = h/(2*fabs(y1_resc));
    
    double lambda = min(lambda_x, lambda_y);
    
    double x_inters=lambda* x1_resc;
    double y_inters=lambda* y1_resc;
    double x_opposite=-lambda* x1_resc;
    double y_opposite=-lambda* y1_resc;
    
    vector<double> intersections { x_inters, y_inters, x_opposite, y_opposite};

    return intersections;   
    
};


    
    
double distance_real(double x1, double y1, double x2, double y2)
{
    
    
    double dist =sqrt( pow( x1 - x2, 2.) + pow( y1 - y2, 2.));
    

    return dist;   
    
};


    

    
