#ifndef class_h
#define class_h
#include<vector>
#include<deque>
#include<cmath>
#include <cstdlib>
#include<random>
#include <algorithm>
#include <utility>
#include <limits>
#include <iomanip>
#include <string>
#include <fstream>
#include <memory>
#include <iostream>
#include <sstream>
#include <chrono>
#include  <utility>
#include <sys/types.h>
#include <sys/stat.h>
#include <csignal>
//#include "fastsum_NFFT_simple.h"

// classes header file

extern int dtorcalls;
using namespace std;

extern mt19937 mt;

class antigenic_coordinate{ // class for an antigenic coordinate.

    
    
    long int num_virs_;// number of   virs at that coord
    long int num_IS_;// number of  immune receptors at that coord
    
    double fitness_; // fitness of the viruses, if no viruses this is set to 0

    
    

public:
	antigenic_coordinate(); 
	antigenic_coordinate(long int num_virs, long int num_IS, double fitness=0.); 
   	antigenic_coordinate(const antigenic_coordinate & c);

	~antigenic_coordinate() { 
        //dtorcalls++;        
        //std::cout << "coord dtor @" << this  << " time "  << dtorcalls << std::endl; 
        //std::cout << "Virus destroy count_infected_ @" << count_infected_ <<" "<< *count_infected_ <<" "<< count_infected_.get() <<" "<< count_infected_.use_count() << std::endl;
    }; //need delete only if points to the heap, that is if initialized with new
    


	long int get_num_vir() const {return num_virs_;};
	long int get_num_IS() const {return num_IS_;};
	
    double get_fitness() const {return fitness_;};
    
    

	void set_num_virs(long int num_virs) {
        if(num_virs<0){
            cout << "num_virs less than 0! "<< num_virs << endl;
            exit(1); 
        }        
        num_virs_= num_virs;
    };
	void set_num_IS(long int num_IS) {
        if(num_IS<0){
            cout << "num_IS less than 0! "<< num_IS << endl;
            exit(1); 
        }
        num_IS_= num_IS;
    };
    
    
    void set_fitness(double fitness) {
        if(fitness<-1){
            cout << "Fitness less than -1! "<< fitness << endl;
            exit(1); 
        }
        fitness_= fitness;};
    
    
        antigenic_coordinate & operator= (const antigenic_coordinate & c);


};  



//antigenic_space class, represents an antigenic_space lattice,  with all the parameters of viral dnamics
class antigenic_space{
    
  
    // <mu> <recognition width> <lattice space> <people number>  <out_dir>  <simulation time> <full_frame save time> <n_reals> <real> <initial_condition> <save_final_configuration> <kernel> <save_phylo> <fake_in_cond_str> <mem_points>
    
    int DIM_; // space dimension
    double recog_width_; // recog_width
    double latt_width_; // lattice spacing
    double jump_size_; // average mutation effect
    long int people_number_; // number of hosts
    int mem_points_; // number of memory receptors per host
    
    string kernel_;
    string out_dir_;
    
    double mu_; // viral mutation rate
    double F0_;

  
    vector<vector<int>> vir_coords_; // coordinates with viruses, in the "original grid", to acces vir_x, do vir_coords[0] - start_x
    vector<vector<int>> IS_coords_;  // coordinates with immune receptors
    vector<vector<long int>> IS_update_coords_;  // coordinates with immune receptors that have been updated each cycle. Holds also the number, instead of embedding it in the coordinates themselves. Therefore 3 columns
 
    long int numvirtot_;
    
    
    int num_x_;//grid dimensions, number of entries
    int num_y_;
    
    int start_x_; // coordinates corresponding to low left corner of current grid, in grid number units
    int start_y_;
    
    int min_vir_x_; // extremes of viral cloud, to be checked and changed whenever a site is emptied or colonized
    int min_vir_y_;
    int max_vir_x_;
    int max_vir_y_;
    int min_IS_x_;
    int min_IS_y_;
    int max_IS_x_;
    int max_IS_y_;
    
    int pad_size_; // pad size around border occupied sites
    
    vector<vector<antigenic_coordinate>> world_; // lattice
    
    
    // some mutation statistics records
    double displ_tot_    ;
    double displ_x_    ;
    double displ_y_    ;
    double displ_x_sq_ ;
    double displ_y_sq_ ;
    int count_displ_;
 
 
public:

	
	antigenic_space ();// 
    
	antigenic_space (vector<vector<int>> & vir_coords,  vector<vector<int>> & IS_coords,  vector<vector<long int>> & IS_update_coords, int DIM, double recog_width, double latt_width, double jump_size, long int people_number, int mem_points, string kernell, double mu, string out_dir,  int num_x,  int num_y,  int start_x,  int start_y, vector<vector<antigenic_coordinate>> & world, double F0, int pad_size, int min_vir_x, int min_vir_y, int max_vir_x, int max_vir_y, int min_IS_x, int min_IS_y, int max_IS_x, int max_IS_y);//  
	antigenic_space(const antigenic_space & h);
	~antigenic_space() {//std::cout << "Host dtor @" << this << std::endl; 
	};
	

	int get_DIM() const {return DIM_;};
	long int  get_people_number() const {return people_number_;};
	int get_mem_points() const {return mem_points_;};
	double get_F0() const {return F0_;};
    
	double get_recog_width() const {return recog_width_;};
	double get_mu() const {return mu_;};
	
    string  get_kernel()  const {return kernel_; };
    string  get_out_dir() const {return out_dir_;};
    
	int get_num_vir_coords() const {return vir_coords_.size();};    
	int get_num_IS_coords()  const {return IS_coords_.size(); };    
	int get_num_IS_update_coords()  const {return IS_update_coords_.size(); };    
	
    
	long int get_numvirtot()  const {return numvirtot_; };    
    
    
	int get_num_x()  const {return  num_x_; };    
	int get_num_y()  const {return  num_y_; };    
	int get_start_x()  const {return  start_x_; };    
	int get_start_y()  const {return  start_y_; };    
	
    int get_min_vir_x()  const {return  min_vir_x_; };    
    int get_min_vir_y()  const {return  min_vir_y_; };    
    int get_max_vir_x()  const {return  max_vir_x_; };    
    int get_max_vir_y()  const {return  max_vir_y_; };    
    int get_min_IS_x()   const {return  min_IS_x_; };    
    int get_min_IS_y()   const {return  min_IS_y_; };    
    int get_max_IS_x()   const {return  max_IS_x_; };    
    int get_max_IS_y()   const {return  max_IS_y_; };    
    
    int get_pad_size()  const {return pad_size_; };    
    
    double idx_to_real_x(int i1) const; //translare world_idx to real coord
    double idx_to_real_y(int j1) const;
    double distance(int i1, int j1, int i2, int j2) const;
    
    
	long int get_tot_IS()  const;    
	long int get_tot_IS_update()  const;    
	long int get_tot_vir() const;    
	vector<double> get_avg_fitness() const;    
	vector<double> get_avg_vir_coord() const;    
    
	int get_maxobs(int obs) const;    // return the index of vir_coords_where some observable is maximum- 0 num_virs, 1 num_is 2 fitn
	vector<int> get_maxobs() const;    // overload, return vector of indexes of vir_coords_where observables are maximum- 0 num_virs, 1 num_is 2 fitn
    
    
    
	void set_num_x  (int num_x  ) 
    {
        if(num_x!=world_.size()){
            cout << "num_x different than grid size! "<< num_x << " world " << world_.size() << endl;
            exit(1); 
        }         
        num_x_  = num_x  ; 
    };    
	void set_num_y  (int num_y  ) 
    {
        if(num_y!=world_[0].size()){
            cout << "num_y different than grid size! "<< num_y << " world " << world_[0].size() << endl;
            exit(1); 
        }    
        num_y_  = num_y  ; 
    };    
	void set_start_x(int start_x) {start_x_= start_x; };    
	void set_start_y(int start_y) {start_y_= start_y; };    
	
    void set_min_vir_x(int min_vir_x) {min_vir_x_= min_vir_x; };    
    void set_min_vir_y(int min_vir_y) {min_vir_y_= min_vir_y; };    
    void set_max_vir_x(int max_vir_x) {max_vir_x_= max_vir_x; };    
    void set_max_vir_y(int max_vir_y) {max_vir_y_= max_vir_y; };    
    void set_min_IS_x(int min_IS_x)   {min_IS_x_ = min_IS_x; };    
    void set_min_IS_y(int min_IS_y)   {min_IS_y_ = min_IS_y; };    
    void set_max_IS_x(int max_IS_x)   {max_IS_x_ = max_IS_x; };    
    void set_max_IS_y(int max_IS_y)   {max_IS_y_ = max_IS_y; };    

    
	void set_out_dir(string out_dir) {out_dir_ = out_dir;}; 
    
    // convolution update

	vector<double> update_fitness_fullconv(bool upd_fitn = true); //
    
    
    // different approximation algorithms to speed up convolution update
	
    vector<double> update_fitness_fullconv_update(bool upd_fitn, bool combined_prep = false, bool combined_prep_newvirs = false, bool approx_fullIS = false, vector<vector<long int>> IS_update_coords_subset = vector<vector<long int>>(), vector<vector<int>> IS_coords_subset= vector<vector<int>>(), double conv_precomp_update=0., double conv_precomp_tot=0.); //
    vector<double> update_fitness_fullconv_update_precomp_thetabins(bool upd_fitn, bool combined_prep = false, bool combined_prep_newvirs = false, bool approx_fullIS = false, vector<vector<long int>> IS_update_coords_subset = vector<vector<long int>>(), vector<vector<int>> IS_coords_subset= vector<vector<int>>(), vector<vector<double>> conv_precomp_update=vector<vector<double>>(), vector<vector<double>>  conv_precomp_tot=vector<vector<double>>()); // overload for thetabins
	
    vector<double> update_fitness_IS_update_NFFT(bool upd_fitn, bool combined_prep = false, bool combined_prep_newvirs = false, bool thetabins = false, bool FFT_newvirs_gen = true, bool approx_fullIS = false, vector<vector<long int>> IS_update_coords_subset = vector<vector<long int>>(), vector<vector<int>> IS_coords_subset= vector<vector<int>>(), vector<vector<double>> conv_precomp_update=vector<vector<double>>(), vector<vector<double>>  conv_precomp_tot=vector<vector<double>>(), vector<int> limits_IS_notapprox_update=vector<int>(), vector<int> limits_IS_notapprox_tot=vector<int>(), int count_old_vir=0, int count_new_vir=0, vector<int> idx_old_vir=vector<int>(), vector<int> idx_new_vir=vector<int>(), bool swap_tot_upd_FFT = false); //
    
    
    void update_fitness_IS_update_prep_farfield(vector<vector<long int>> & IS_update_coords_notapprox, double & conv_precomp_update, vector<int> & limits_IS_notapprox_update, int & count_old_vir, int & count_new_vir, vector<int> & idx_old_vir, vector<int> & idx_new_vir, bool approx_fullIS = false); //
    
    void update_fitness_IS_update_prep_farfield_thetabins(bool combined_prep_newvirs, vector<vector<long int>> & IS_update_coords_notapprox, vector<vector<int>> & IS_coords_notapprox, vector<vector<double>> & conv_precomp_update, vector<vector<double>> & conv_precomp_tot, vector<int> & limits_IS_notapprox_update, vector<int> & limits_IS_notapprox_tot, int & count_old_vir, int & count_new_vir, vector<int> & idx_old_vir, vector<int> & idx_new_vir, bool approx_fullIS = false); //
    
    vector<double> update_fitness_IS_update_combine_farfield_fullconv(bool upd_fitn, bool thetabins=false,  bool combined_prep_newvirs=false, bool approx_fullIS = false); //

    vector<double> update_fitness_IS_update_combine_farfield_NFFT(bool upd_fitn, bool thetabins=false,  bool combined_prep_newvirs=false,  bool FFT_newvirs_gen=true, bool approx_fullIS = false); //
    
    vector<double> update_fitness_IS_update_switch(bool upd_fitn, double & maxfiterr_o_std); // choose convolution algorithm runtime, based on estimated complexity
    
	
    double benchmark_update_fitness_fullconv_update_NFFT_farfield(); //
    
    
    // cycles steps
	void update_IS(); //
	void growth_virs(); //
	void diffusion_virs(); //
	//~ void growth_plus_diffusion_virs(); //
    
    
    
    // handle lattice cropping it around viruses
    void expand_world(int dir); //
    void shrink_world(int dir, int shrinksize); //
    
    void checkpad_and_shrink( int x, int y); //
    
    
    int check_new_border( int x, int y, int bord, bool vir); // if remove vir or IS coord, and it's a border, {xmin, xmax, ymin,ymax} need to check for the new border






	void print_to_file(double t) const; // viruses snapshot
    void print_to_file_IS_update(double t) const; // IS subsample, those updated this cycle
	
    void print_to_file_fullstatus(bool first_save, double t) const; // checkpoint
    
    
    
	void print_avg_to_file(fstream& outstream, double t); // time series averages
    
    
    antigenic_space & operator= (const antigenic_space& h);



    vector<antigenic_coordinate>& operator[] (const int i); // gets the coordiates vector i . 
    const vector<antigenic_coordinate>& operator[] (const int i) const; //  


    void infection_cycle(int & count_full_upd, bool benchmark_update_fitness=false); // perform the infection cycle 
    
};


antigenic_space create_world(int DIM, double recog_width, double latt_width, double jump_size, long int people_number, int mem_points, string kernel, double mu, string out_dir, double F0); // populates the lattice according to some initial condition

double fitness_fct(double conv, int M, double F0);
double fitness_inv_fct(double fitn, int M, double F0);

double kern( double distance, double rec, string ker); // kernell evaluation

double distance_real(double x1, double y1, double x2, double y2); // distance in real world GIVEN IDX
vector<double> intersections_rect(double x1, double y1, double h, double w);





template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random, int curr_untor_t_inf_end=-1) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
        uniform_int_distribution<int> distr_int(0,left-1);
        bidiiter last = end - 1;
        bidiiter r = last;
	int ind;
	//bidiiter new= last;
	//do {
	    ind=distr_int(mt);
	    //advance(new, -ind);
	//}while((last - ind)->get_t_inf_end()>=curr_untor_t_inf_end);
        advance(r, -ind);
        swap(*last, *r);
        --end;
        --left;
    }
    return end;
}

inline const char * const BoolToString(bool b)
{
  return b ? "true" : "false";
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

#endif





