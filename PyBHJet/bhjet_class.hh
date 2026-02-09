#ifndef PYBHJET_CLASS
#define PYBHJET_CLASS

#include "jet_output.hh" 
#include <unordered_map>
#include <vector>
#include <string>

class BhJetClass {
public:
    BhJetClass(); 

    void load_params(const std::string& file);
    void print_parameters() const; 
    void run();
    void run_singlezone();
    const JetOutput& get_output() const;

    //Accessing parameters by name in python 
    double get_parameter(const std::string& name) const;
    void set_parameter(const std::string& name, double value);

    int cutoff_type = 0;
    int get_cutoff_type() const {return cutoff_type; }
    void set_cutoff_type(int t) {cutoff_type = t; }

    // expose parameter names to Python
    std::vector<std::string> get_parameter_names() const;

    double Mbh, Eddlum, Rg, theta, dist, redsh, jetrat, zmin, r_0, h, z_acc, z_diss, z_max, t_e;
    double f_nth, f_pl, pspec, f_heat, f_beta, f_sc, p_beta, sig_acc, l_disk, r_in, r_out;
    double compar1, compar2, compar3, compsw, velsw;
    int infosw, EBLsw;
    std::vector<double> energy_grid,total_flux_vals; 

private:
    int npar, ne;
    double emin, emax;
    bool params_loaded = false; //Checking if parameters were loaded first before running code 

    std::vector<double> params;
    std::vector<std::pair<std::string, std::string>> param_units; // Add units map
    std::unordered_map<std::string, size_t> param_name_to_index; // Map of parameter names to indices
    
    void initialize_parameter_map();
    void initialize_parameter_units();

    JetOutput output;  // JetOutput instance to store results, maybe should change to be more detailed per output type: 

    void update_internal_parameters();
};

#endif // PYBHJET_CLASS