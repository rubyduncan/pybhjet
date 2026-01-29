#include "bhjet_class.hh"
#include "jetmain.hh"

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <memory>

using namespace std;

BhJetClass::BhJetClass() : params(28, 0.0) {
    //initializing a vector with 28 elements 
    initialize_parameter_map(); //setting up mapping between the parameter names and their indicies in the file 
    initialize_parameter_units();
}

// read parameters from a file (like .dat, normal) from some file 
void BhJetClass::load_params(const std::string& file) {
    std::ifstream inFile(file); //opening file 
    if (!inFile) {
        throw std::runtime_error("Cannot open parameter file: " + file);
    }

    std::string line;
    int line_nb = 0;
    while (std::getline(inFile, line)) {
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char c) { return !std::isspace(c); }));
        if (line.empty() || line[0] == '#') {
            continue;
        }
        params[line_nb] = std::atof(line.c_str()); //stores new values in the param vector
        line_nb++;
    }
    inFile.close();

    update_internal_parameters(); //now updates the member variable w/ new param value 
    params_loaded = true; //for the run function 
}


void BhJetClass::initialize_parameter_map() {  //how the params can be accessed by their name instead of indicies
    param_name_to_index = {
        {"Mbh", 0},
        {"theta", 1},
        {"dist", 2},
        {"redsh", 3},
        {"jetrat", 4},
        {"r_0", 5},
        {"z_diss", 6},
        {"z_acc", 7},
        {"z_max", 8},
        {"t_e", 9},
        {"f_nth", 10},
        {"f_pl", 11},
        {"pspec", 12},
        {"f_heat", 13},
        {"f_beta", 14},
        {"f_sc", 15},
        {"p_beta", 16},
        {"sig_acc", 17},
        {"l_disk", 18},
        {"r_in", 19},
        {"r_out", 20},
        {"compar1", 21},
        {"compar2", 22},
        {"compar3", 23},
        {"compsw", 24},
        {"velsw", 25},
        {"infosw", 26},
        {"EBLsw", 27}
    };
}

void BhJetClass::initialize_parameter_units() {
    param_units = {
        {"Mbh", "M_sun"},
        {"theta", "degrees"},
        {"dist", "kpc"},
        {"redsh", ""},
        {"jetrat", "Jet Power, L_edd"},
        {"r_0", "Jet base radius, R_g"},
        {"z_diss", "nT particle injection region, R_g"},
        {"z_acc", "Jet acceleration region, R_g"},
        {"z_max", "Maximum Jet Distance, R_g"},
        {"t_e", "Temperature of rel. electrons in nozzle, keV"},
        {"f_nth", "Percentage of particles acc. into nT tail "},
        {"f_pl", "Dissipation Parameter"},
        {"pspec", "Slope of nT lepton distribution"},
        {"f_heat", "Shock Heating"},
        {"f_beta", "Adiabatic Cooling Timescale"},
        {"f_sc", " < 0.1, Max. E of nT particles"},
        {"p_beta", "Plasma Beta [Regulates Pair Content in BHJet for velsw > 1] "},
        {"sig_acc", "Magnetization at acceleration region"},
        {"l_disk", "Disk Luminosity, L_edd"},
        {"r_in", "Inner Disk Radius, R_g"},
        {"r_out", "Outer Disk Radius, R_g"},
        {"compar1", "Compsw=1: Temp of BB [K], Compsw=2: % of disk photons reprocessed by BLR"},
        {"compar2", "Compsw=1: Lum of BB [ergs], Compsw=2: % of disk photons reprocessed by torus"},
        {"compar3", "Compsw=1: BB Energy Density"},
        {"compsw", "Adds Ext. Photon Field, = 1: BB, = 2: BLR & Torus"},
        {"velsw", "=0, 1: agnjet(pressure driven), > 1: bljet(magnetic driven)"},
        {"infosw", ""},
        {"EBLsw", "Account for EBL attenuation"}
    };
}

void BhJetClass::print_parameters() const {
    std::cout << "BHJet Parameters:" << std::endl;
    std::cout << "-----------------" << std::endl;

    // Determine the maximum width for parameter names
    size_t max_name_length = 0;
    for (const auto& pair : param_units) {
        max_name_length = std::max(max_name_length, pair.first.length());
    }

    // Print the parameters with alignment
    for (const auto& pair : param_units) {
        const auto& name = pair.first;
        const auto& unit = pair.second;
        auto it = param_name_to_index.find(name);
        if (it != param_name_to_index.end()) {
            size_t index = it->second;

            // Align the parameter name and value
            std::cout << std::left << std::setw(max_name_length + 2) << name
                      << ": " << std::setw(12) << std::left << params[index]
                      << unit << std::endl;
        }
    }
}

// assigning the names and values of parameters  
void BhJetClass::update_internal_parameters() {
    Mbh = params[0];
    theta = params[1];
    dist = params[2]; // in kpc;
    redsh = params[3];
    jetrat = params[4]; // in edd lum 
    r_0 = params[5]; // in  Rg;
    z_diss = params[6]; // in  Rg;
    z_acc = params[7]; // in  Rg;
    z_max = params[8]; // in  Rg;
    t_e = params[9];
    f_nth = params[10];
    f_pl = params[11];
    pspec = params[12];
    f_heat = params[13];
    f_beta = params[14];
    f_sc = params[15];
    p_beta = params[16];
    sig_acc = params[17];
    l_disk = params[18];
    r_in = params[19]; // in rg
    r_out = params[20]; // in rg;
    compar1 = params[21]; 
    compar2 = params[22];
    compar3 = params[23];
    compsw = params[24];
    velsw = params[25];
    infosw = static_cast<int>(params[26]);
    EBLsw = static_cast<int>(params[27]);
}

double BhJetClass::get_parameter(const std::string& name) const {
    auto it = param_name_to_index.find(name); //using the param map created above now 
    if (it != param_name_to_index.end()) {
        size_t index = it->second;
        return params[index];
    } else {
        throw std::invalid_argument("Parameter name not found: " + name);
    }
}

// this takes a parameter name and a new value to assign to it 
void BhJetClass::set_parameter(const std::string& name, double value) {
    auto it = param_name_to_index.find(name); //param map 
    if (it != param_name_to_index.end()) {
        size_t index = it->second;
        params[index] = value;
        update_internal_parameters();  // Update dependent variables
    } else {
        throw std::invalid_argument("Parameter name not found: " + name);
    }
}

std::vector<std::string> BhJetClass::get_parameter_names() const {
    std::vector<std::string> names;
    names.reserve(param_name_to_index.size());
    for (const auto& kv : param_name_to_index) {
        names.push_back(kv.first);
    }
    return names;
}

//kinda weird, I think this works, for getting the populated jet output to python  -- 
const JetOutput& BhJetClass::get_output() const {
    return output;
}

void BhJetClass::run() {
    // if (!params_loaded) {
    //     throw std::runtime_error("Parameters have not been loaded. Please call load_params() first.");
    // }

    // auto start_time = std::chrono::high_resolution_clock::now();

    //this is what was used in the bhwrap file for running bhjet alone ---- 
    int npar = 28;
    int ne = 201;
    double emin	= -9;
    double emax	= 3;
    double einc	= (emax-emin)/ne;

    auto ebins = std::make_unique<double[]>(ne); //energy bins 
    auto spec = std::make_unique<double[]>(ne - 1); // photon energies 
    auto dumarr = std::make_unique<double[]>(ne - 1); // dummy photon spectrum, if needed 

    //fill the energy bins 
    for (int i = 0; i < ne; ++i) { 
        ebins[i] = std::pow(10, emin + i * einc);
    }

    output.clear();

    // run the jetmain function: 
    jetmain(*this, ebins.get(), ne - 1, spec.get(), dumarr.get(), output);

    // Stop the timer
    // auto end_time = std::chrono::high_resolution_clock::now();

    // // Calculate and print elapsed time
    // std::chrono::duration<double> elapsed_time = end_time - start_time;
    // std::cout << "run() execution time: " << elapsed_time.count() << " seconds" << std::endl;

}


void BhJetClass::run_singlezone(){

    if (!params_loaded) {
        throw std::runtime_error("Parameters have not been loaded. Please call load_params() first.");
    }

    output.clear();

    singlezone_jetmain(*this,output); 
}


