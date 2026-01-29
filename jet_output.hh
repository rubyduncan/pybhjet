#ifndef JET_OUTPUT_HH
#define JET_OUTPUT_HH

// #include "jetmain.hh"
#include <vector>
#include <string>

// Define DataPoint struct for 2 value outputs (emission components, zones )
struct DataPoint {
    double energy; // nu [Hz]
    double flux;   // flux [mJy]
};

struct NumDenPoint {
    double momentum;
    double gamma;
    double n_p;
    double n_g;
};


class JetOutput {
public:
    
    // infosw = 1 ---------
    //all output here is in units of (nu [Hz], flux [mJy])
    std::vector<DataPoint> presyn;
    std::vector<DataPoint> postsyn;
    std::vector<DataPoint> precom;
    std::vector<DataPoint> postcom;
    std::vector<DataPoint> disk;
    std::vector<DataPoint> bb;
    std::vector<DataPoint> total;

    // For infosw >= 2 ------
    
    //output is (p [g cm s-1], g [], n(p) [# cm^-3 p^-1], n(g) [# cm^-3 g^-1])
    std::vector<NumDenPoint> numdens;

    //for infosw >=5 ---- 
    struct JetProfile {
        std::vector<double> z_rg; 
        std::vector<double> zone_rg; 
        std::vector<double> zone_bfield; 
        std::vector<double> zone_lepdens; 
        std::vector<double> zone_gamma; 
        std::vector<double> zone_eltemp; 

        void clear() {
            z_rg.clear();
            zone_rg.clear();
            zone_bfield.clear();
            zone_lepdens.clear();
            zone_gamma.clear();
            zone_eltemp.clear();
        }
    };
    JetProfile jetprofile;  
    
    // output for these is in units of (nu [Hz], flux [mJy])
    std::vector<DataPoint> cyclosyn_zones;
    std::vector<DataPoint> compton_zones;

    // For infosw >= 3 --------
    struct SpectralProperties {
        std::vector<double> disk_lum;  //Observed 0.3-5 keV disk luminosity:
        std::vector<double> IC_lum;   //Observed 0.3-300 keV Inverse Compton luminosity
        std::vector<double> xray_lum; //Observed 1-10 keV total luminosity
        std::vector<double> radio_lum;  //Observed 4-6 GHz luminosity
        std::vector<double> xray_index; //X-ray 10-100 keV photon index estimate
        std::vector<double> radio_index; //Radio 10-100 GHz spectral index estimate
        std::vector<double> jetbase_compactness; //Jet base compactness

        void clear() {
            disk_lum.clear(); 
            IC_lum.clear(); 
            xray_lum.clear();
            radio_lum.clear();
            xray_index.clear();
            radio_index.clear(); 
            jetbase_compactness.clear();
        }
    
    };
    SpectralProperties spectral_properties;

    struct JetBaseProperties {
        std::vector<double> pair_content;
        std::vector<double> init_mag;
        std::vector<double> particle_avg_lorentz_factor;
        std::vector<double> jet_nozzle_end;
        std::vector<double> jet_nozzle_optical_depth;

        void clear() {
            pair_content.clear();
            init_mag.clear();
            particle_avg_lorentz_factor.clear();
            jet_nozzle_end.clear();
            jet_nozzle_optical_depth.clear();
        }
    };
    JetBaseProperties jet_base_properties;

    // For infosw >= 5 ---- 
    //this returns values for each zone, as the code loops over each segement of the jet: 
    struct JetZoneProperties {
        std::vector<double> jet_bfield;
        std::vector<double> lepton_ndens; 
        std::vector<double> speed_gamma; 
        std::vector<double> delta; 
        std::vector<double> tshift; 
        std::vector<double> temp_kev; 
        std::vector<double> grid_r;  
        std::vector<double> delz; 
        std::vector<double> dist_z;  
        std::vector<double> z_delz; 
        std::vector<double> equpar_check; 
        std::vector<double> ue_ub; 

        // Method to clear all vectors
        void clear() {
            jet_bfield.clear();
            lepton_ndens.clear();
            speed_gamma.clear();
            delta.clear();
            tshift.clear();
            temp_kev.clear();
            grid_r.clear();
            delz.clear();
            dist_z.clear();
            z_delz.clear();
            equpar_check.clear();
            ue_ub.clear();
        }
    };
    JetZoneProperties jet_zone_properties;

    // Constructor
    JetOutput() = default;

    // Clear method
    void clear() {
        // Clear all vectors
        presyn.clear();
        postsyn.clear();
        precom.clear();
        postcom.clear();
        disk.clear();
        bb.clear();
        total.clear();
        cyclosyn_zones.clear();
        compton_zones.clear();
        jetprofile.clear();
        numdens.clear();
        spectral_properties.clear();
        jet_base_properties.clear();
        jet_zone_properties.clear();

    }
};

#endif // JET_OUTPUT_H