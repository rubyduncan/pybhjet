#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <cmath>

#define cee             GSL_CONST_CGSM_SPEED_OF_LIGHT
#define emgm            GSL_CONST_CGSM_MASS_ELECTRON
#define pmgm            GSL_CONST_CGSM_MASS_PROTON
#define kboltz          GSL_CONST_CGSM_BOLTZMANN
#define kboltz_kev2erg  1.6022e-9
#define gr_to_kev       5.6095883571872e+29
#define me_kev          511.
#define emerg           (GSL_CONST_CGSM_MASS_ELECTRON*pow(GSL_CONST_CGSM_SPEED_OF_LIGHT,2.))
#define pi              M_PI
#define charg           4.8e-10
#define sigtom          GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define herg            GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define hkev            (GSL_CONST_CGSM_PLANCKS_CONSTANT_H*6.2415e8)
#define mjy             1.e-26
#define re0             2.81794e-13
#define gconst          GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT
#define sbconst         GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT
#define aconst          7.56e-15

//Template class for particle distributions
//This class contains members and methods that are used for thermal, non-thermal and mixed distributions

//Structures used for GSL integration
typedef struct pl_params { 
    double s; 
    double n;
} pl_params;

typedef struct bkn_params { 
    double s1;
    double s2;
    double brk;
    double max; 
    double m;
    int cutoff_type;
} bkn_params;

typedef struct th_params { 
    double t; 
    double n;
    double m; 
} th_params;

typedef struct k_params { 
    double t;
    double k; 
} k_params;

typedef struct injection_mixed_params {
    double s; 
    double t; 
    double nth;
    double npl;
    double m; 
    double min;
    double max;
    double cutoff;
    int cutoff_type; 
} injection_mixed_params;

typedef struct injection_kappa_params {
    double t; 
    double k;
    double n; 
    double m;
} injection_kappa_params;

typedef struct injection_pl_params {
    double s;
    double n; 
    double m; 
    double max;
    int cutoff_type;
} injection_pl_params;

typedef struct injection_bkn_params {
    double s1;
    double s2;
    double brk;
    double max; 
    double m;
    double n;
    int cutoff_type;
} injection_bkn_params;

class Particles {
    protected:
        int size;
        
        double mass_gr;         //particle mass in grams
        double mass_kev;        //same as above but in keV, using electrons as "reference"
        
        double *p;				//array of particle momenta
        double *ndens;			//array of number density per unit volume, per unit momentum
        double *gamma;			//array of particle kinetic energies for each momentum
        double *gdens;		    //array of number density per unit volume, per unit gamma
        double *gdens_diff;     //array with differential of number density for radiation calculation

        //for the cutoff function
        int cutoff_type = 0;    // setting cutoff type within powerlaw/mixed 
						        
    public:
        ~Particles();
        
        static inline double cutoff_factor(double x, int type); //static helper for cutoff definitions ??
        int get_cutoff_type()const              { return cutoff_type; } // cutoff getter 
        void set_cutoff_type(int t)             { cutoff_type = t; }  //setter 
        
        void set_mass(double m);
        void initialize_gdens();	
        void initialize_pdens();
        void gdens_differentiate();

        const double *get_p()const              { return p; }
        const double *get_pdens()const          { return ndens; }
        const double *get_gamma()const          { return gamma; }
        const double *get_gdens()const          { return gdens; }
        const double *get_gdens_diff()const     { return gdens_diff; }
			
        double beta(int i);							

        double count_particles();
        double count_particles_energy();
        double av_p();
        double av_gamma();
        double av_psq();
        double av_gammasq();

        void test_arrays();	
        
};

//need this to be accessible by both mixed and powerlaw, so here goes: 
// this has three options for the behavior of the cutoff in the electron spectrum: 
// 1. the traditional exponential cutoff exp(e/ecut)
// 2. this is from Comisso 2021, exp[(e/ecut)^2] - magnetic turbulence PIC sim 
// 3. also comisso 2021, sech[(e/ecut)^2]

// in principle, meant to act like this: (in set_ndens or injection_mixed_int), in mixed.cpp or powerlaw.cpp: 
// C = cutoff_factor(p[i]/pmax_pl, cutoff_type); then ndens = npl * mom_int thing * C

// Three cutoff types, with x = p/p_cut (momentum)
//cutoff type 0 = exp(-x)
// cutoff type 1 = exp(-x^2)

// cutoff type 2 = sech(x)^2 = 1/cosh(x)^2 
// cosh(x) = ( e^x + e^-x ) / 2
// 
inline double Particles::cutoff_factor(double x, int type) {
    if (type == 0) return std::exp(-x); // classic exponential cutoff 
    if (type == 1) return std::exp(-(x * x)); // 
    if (type == 2) {    // sech^2 cutoff
        const double ax = std::fabs(x); 
        if (ax < 30.0) { //cosh(30) ~ 5e12 
            const double c = std::cosh(ax); //directly calculate sech(x)
            return 1.0 / (c * c);
        } else {
            const double t = std::exp(-2.0 * ax); //approx
            const double denom = 1.0 + t;
            return 4.0 * t / (denom * denom);
        }
    }
    return std::exp(-x); 
}

#endif
