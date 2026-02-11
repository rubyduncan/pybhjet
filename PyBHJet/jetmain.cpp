#include "jetmain.hh"

#include <fstream>
#include <stdarg.h>

using namespace std;

void jetmain(BhJetClass& bhjet, double* ear, int ne, double* photeng, double* photspec, JetOutput& output) {

    bool IsShock = false;						//flag to set shock heating

    int nz = 100;								//total number of zones
    int nel = 70;
    int syn_res = 10;							//number of bins per decade in synch frequency;
    int com_res = 6;							//number of bins per decade in compton frequency;
    int nsyn,ncom;								//number of bins in synch/compton frequency;
    int npsw = 1;								//switch to define number of protons calculations in agnjet

    // Use named variables directly 
    double Mbh = bhjet.Mbh;
    double Eddlum = 1.25e38 * Mbh;
    double Rg = gconst * Mbh * msun / (cee * cee);
    double theta = bhjet.theta;
    double dist = bhjet.dist * kpc; 
    double redsh = bhjet.redsh;
    double jetrat = bhjet.jetrat * Eddlum; 
    double r_0 = bhjet.r_0 * Rg; 
    double z_diss = bhjet.z_diss * Rg; 
    double z_acc = bhjet.z_acc * Rg;
    double z_max = bhjet.z_max * Rg;
    double t_e = bhjet.t_e;
    double f_nth = bhjet.f_nth;
    double f_pl = bhjet.f_pl;
    double pspec = bhjet.pspec;
    double f_heat = bhjet.f_heat;
    double f_beta = bhjet.f_beta;
    double f_sc = bhjet.f_sc;
    double p_beta = bhjet.p_beta;
    double sig_acc = bhjet.sig_acc;
    double l_disk = bhjet.l_disk;
    double r_in = bhjet.r_in * Rg; 
    double r_out = bhjet.r_out * Rg; 
    double compar1 = bhjet.compar1;
    double compar2 = bhjet.compar2;
    double compar3 = bhjet.compar3;
    double compsw = bhjet.compsw;
    double velsw = bhjet.velsw;
    int infosw = bhjet.infosw;
    int EBLsw = bhjet.EBLsw;
    double zmin = 2.0 * Rg;

    double z;									//distance along the jet axis
    double tshift;								//temperature shift from initial value due to ad. cooling	
    double Urad;								//estimate of total radiation energy density in each zone
    double Ubb1,Ubb2;							//estimate of comoving energy density of external photons
    double gmin,gmax;							//minimum/maximum Lorentz factors over which to integrate
    double syn_min,syn_max;						//interval for synchrotron calculation in each zone
    double com_min,com_max;						//interval for inverse Compton calculation in each zone

    double *syn_en;								//sychrotron energy array for jet+counterjet summed
    double *syn_lum;							//synchrotron luminosity array for jet+counterjet summed
    double *com_en;								//compton energy array for jet+counterjet summed
    double *com_lum;							//compton luminosity array for jet+counterjet summed

    double *tot_en = new double[ne];			//energy arrray for sum of all zones and/or components
    double *tot_syn_pre = new double[ne];		//specific synchrotron luminosity arrays for all zones 
    double *tot_syn_post = new double[ne];		//pre/post particle	acceleration
    double *tot_com_pre = new double[ne];		//same as above but for the inverse Compton part
    double *tot_com_post = new double[ne];
    double *tot_lum = new double[ne];			//specific luminosity array for sum of all components	

    //these are all defined in the header file: 
    grid_pars grid;								//structure with grid parameters
    jet_dynpars jet_dyn;						//structure with jet dynamical parameters
    jet_enpars nozzle_ener;						//structure with jet energetic parameters
    zone_pars zone;								//strucutre with parameters of each individual zone
    com_pars agn_com;							//structure with parameters for inverse Compton fields in AGN
    
    //External photon object declarations - from Kariba 
    ShSDisk Disk;
    BBody BLR;
    BBody Torus;
    BBody BlackBody;

    //splines for jet acceleration
    gsl_interp_accel *acc_speed = gsl_interp_accel_alloc();
    gsl_spline *spline_speed = gsl_spline_alloc(gsl_interp_steffen,54);

    //splines for electron distribution
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);

    gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
    gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  


    //initialize total energy/luminosity arrays
    for(int i=0;i<ne;i++){
       	tot_en[i] = (ear[i] + (ear[i+1]-ear[i])/2.)*herg/hkev;	
       	tot_syn_pre[i] = 1.;
       	tot_syn_post[i] = 1.;
       	tot_com_pre[i] = 1.;
       	tot_com_post[i] = 1.;
       	tot_lum[i] = 1.;	
    }

//STEP 3: DISK/EXTERNAL PHOTON CALCULATIONS

    if(r_in<r_out){
    	Disk.set_mbh(Mbh);
	    Disk.set_rin(r_in);
	    Disk.set_rout(r_out);
	    Disk.set_luminosity(abs(l_disk));		
	    Disk.set_inclination(theta);
	    Disk.disk_spectrum();
	    if (compsw != 2 && l_disk > 0) {
	        sum_ext(50,ne,Disk.get_energy_obs(),Disk.get_nphot_obs(),tot_en,tot_lum);   	    
	    }			
        if(infosw >= 4) {
            Disk.test();
            cout << endl;
        }	
    }

    //Depending on the value of compsw, we either include a) an extra homogeneous black body in every zone or
    //b) the radiation field of the broad line region/torus of a bright AGN. For bright AGN, the reprocessed
    //fraction of disk luminosity is removed from the observed disk luminosity
    if (compsw==1){
        BlackBody.set_temp_k(compar1);
        BlackBody.set_lum(compar2);
        Ubb1 = compar3;	
        BlackBody.bb_spectrum();
        sum_ext(40,ne,BlackBody.get_energy_obs(),BlackBody.get_nphot_obs(),tot_en,tot_lum);
    }
    else if (compsw==2 && r_in<r_out){
        agn_photons_init(Disk.total_luminosity(),compar1,compar2,agn_com);

        BLR.set_temp_kev(agn_com.tblr);
        BLR.set_lum(agn_com.lblr);
        BLR.bb_spectrum();			

        Torus.set_temp_k(agn_com.tdt);
        Torus.set_lum(agn_com.ldt);	
        Torus.bb_spectrum();
        
        Disk.cover_disk(compar1+compar2);	

        if(infosw>=3){
            cout << "BLR radius in Rg: " << agn_com.rblr/Rg << " and in cm: " << agn_com.rblr<< endl;
            cout << "DT radius in Rg: " << agn_com.rdt/Rg  << " and in cm: " << agn_com.rdt<< endl;
        }
        sum_ext(40,ne,Torus.get_energy_obs(),Torus.get_nphot_obs(),tot_en,tot_lum);
        if (l_disk > 0) {
            sum_ext(50,ne,Disk.get_energy_obs(),Disk.get_nphot_obs(),tot_en,tot_lum);    
        }          		
    }	

//STEP 4: JET BASE EQUIPARTITION CALCULATIONS AND SETUP
    //Dummy particle distribution, needed for average lorentz factor in equipartition function
    //The number density is just set to unity, the normalisation is not needed to calculate the average Lorenz factor
    //of the thermal distribution anyway

    Thermal dummy_elec(nel);
    dummy_elec.set_temp_kev(t_e);
    dummy_elec.set_p();	
    dummy_elec.set_norm(1.);	
    dummy_elec.set_ndens();

    grid.nz = nz;
    grid.cut = 0;
    grid.zcut = 1.e3*Rg;

    jet_dyn.min = zmin;
    jet_dyn.max = z_max;
    jet_dyn.h0 = 2.*r_0+zmin;
    jet_dyn.r0 = r_0;
    jet_dyn.acc = z_acc;
    jet_dyn.beta0 = sqrt(4./3.*(4./3.-1.)/(4./3.+1.));	//set initial jet speed for relativistic fluid, g=4/3 
    jet_dyn.gam0 = 1./sqrt(1.-(pow(jet_dyn.beta0,2.)));	//set corresponding lorentz factor
    jet_dyn.gamf = velsw;
    jet_dyn.Rg = Rg;  

    nozzle_ener.pbeta = p_beta;
    nozzle_ener.Nj = jetrat;
    nozzle_ener.sig_acc = sig_acc; 
    nozzle_ener.av_gamma = dummy_elec.av_gamma();

    //set up jet velocity profile depending on choice of adiabatic,isothermal,magnetically dominated jet   
    //note: the adiabatic jet only runs correctly if the final temperature is above ~1kev, which means the
    //initial temperature has to be ~10^4 kev to avoid numerical issues
    if(velsw==0){   
        velprof_ad(spline_speed);
        equipartition(npsw,jet_dyn,nozzle_ener);	
    } else if(velsw==1){
      	velprof_iso(spline_speed);
      	equipartition(npsw,jet_dyn,nozzle_ener);
    } else {
        velprof_mag(jet_dyn,spline_speed);
        equipartition(jetrat,jet_dyn,nozzle_ener);
    }	

    //check that the pair content is not negative, and also if running bljet that it's not too high
    if(nozzle_ener.eta<1){
        cout << "Unphysical pair content: " << nozzle_ener.eta << " pairs per proton. Check the value of " <<
                "plasma beta!" << endl;
    } else if (velsw>1 && dummy_elec.av_gamma()*nozzle_ener.eta >= 3e2){
           cout << "Pair content or temperature too high for  for bljet! " << endl;
           cout << "Pair content: " << nozzle_ener.eta << " pairs per proton" << endl;
           cout << "Average lepton Lorenz factor: " << dummy_elec.av_gamma() << endl;
           cout << "Check the value of Te and/or plasma beta!" << endl;
    }

    if(infosw>=3){
        // cout << "Jet base parameters: " << endl;
        // cout << "Pair content (ne/np): " << nozzle_ener.eta << endl;
        // cout << "Initial magnetization: " << nozzle_ener.sig0 << endl;
        // cout << "Particle average Lorenz factor: " << dummy_elec.av_gamma() << endl;
        // cout << "Jet nozzle ends at: " << jet_dyn.h0/Rg << " Rg" << endl ;
        // cout << "Jet nozzle optical depth: " << jet_dyn.r0*nozzle_ener.lepdens*sigtom << endl << endl;


        // want to save everything out, rather than having to tell bhjet to do so also 
        output.jet_base_properties.pair_content.push_back(nozzle_ener.eta);
        output.jet_base_properties.init_mag.push_back(nozzle_ener.sig0);
        output.jet_base_properties.particle_avg_lorentz_factor.push_back(dummy_elec.av_gamma());
        output.jet_base_properties.jet_nozzle_end.push_back(jet_dyn.h0 / Rg);
        output.jet_base_properties.jet_nozzle_optical_depth.push_back(jet_dyn.r0 * nozzle_ener.lepdens * sigtom);
        }



//STEP 5: TOTAL JET CALCULATIONS, LOOPING OVER EACH SEGMENT OF THE JET

    for(int i=0;i<nz;i++){

        //calculate dynamics/energetics in each zone
        jetgrid(i,grid,jet_dyn,zone.r,zone.delz,z);	
        if(velsw==0){
            adjetpars(z,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);	 
        } else if (velsw==1){
            isojetpars(z,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);
        } else {
            bljetpars(z,f_beta,jet_dyn,nozzle_ener,tshift,zone,spline_speed,acc_speed);
        }
        zone.delta = 1./(zone.gamma*(1.-zone.beta*cos(theta*pi/180.)));

        //This is to avoid crashes due to low (sub 1 kev) particle temperatures
        if (z < z_diss) {
            zone.eltemp = max(tshift*t_e,1.);    
        } else {
            zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
        }
                
        //This is to evolve the fraction of non thermal particles along the jet, and change the distribution 
        //of non thermal particles appropriately
        if (z < z_diss) {
            zone.nth_frac = 0.;
        } else {
            zone.nth_frac = f_nth*pow(log10(z_diss)/log10(z),f_pl);           
        }

        if (r_in < r_out) {
            double Rdisk = pow(r_in,2.)+pow(z,2.);
            double delta_disk, theta_disk;
            theta_disk = pi-atan(r_in/z);
            delta_disk = 1./(zone.gamma-zone.beta*cos(theta_disk));
            Urad = pow(delta_disk,2.)*l_disk*Eddlum/(4.*pi*Rdisk*cee);
        } else {
            Urad = 0.;
        }
    
        if(compsw==1){			
            Urad = Urad + pow(zone.delta,2.)*Ubb1;
        } else if (compsw==2 && r_in<r_out){
            zone_agn_phfields(z,zone,Ubb1,Ubb2,agn_com);
            Urad = Urad + agn_com.urad_total;
        }

        if(zone.nth_frac == 0.){		
            Thermal th_lep(nel);
            th_lep.set_temp_kev(zone.eltemp);
            th_lep.set_p();
            th_lep.set_norm(zone.lepdens);
            th_lep.set_ndens();

            gmin = th_lep.get_gamma()[0];
            gmax = th_lep.get_gamma()[nel-1];

            zone.avgammasq = pow(th_lep.av_gamma(),2.);

            gsl_spline_init(spline_eldis,th_lep.get_gamma(),th_lep.get_gdens(),nel);
            gsl_spline_init(spline_deriv,th_lep.get_gamma(),th_lep.get_gdens_diff(),nel); 

        if (infosw >=2){
        //saving the number densities 
            store_numdens(nel, th_lep.get_p(),th_lep.get_gamma(),th_lep.get_pdens(),th_lep.get_gdens(), output.numdens); 
        }

        } else if (zone.nth_frac < 0.5){ 
            if (IsShock==false){
                t_e = f_heat*t_e;
                IsShock = true;
                zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
            }					
            Mixed acc_lep(nel);
            acc_lep.set_cutoff_type(bhjet.get_cutoff_type());
            acc_lep.set_temp_kev(zone.eltemp);
            acc_lep.set_pspec(pspec);
            acc_lep.set_plfrac(zone.nth_frac);
            // std::cout << "gmax before:" << acc_lep.get_gamma()[nel-1] << std::endl;
            
            //if f_sc < 10 it's the acceleration efficiency, else it's the desired maximum lorentz factor
            if (f_sc<10.){
                acc_lep.set_p(Urad,zone.bfield,f_beta,zone.r,f_sc);
            } else{
                acc_lep.set_p(f_sc);
            }	

            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,f_beta);
            //Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel-1];
            // std::cout << "gmax (emax) calc from egrid:" << gmax << std::endl;
            
            zone.avgammasq = pow(acc_lep.av_gamma(),2.);

            gsl_spline_init(spline_eldis,acc_lep.get_gamma(),acc_lep.get_gdens(),nel);
            gsl_spline_init(spline_deriv,acc_lep.get_gamma(),acc_lep.get_gdens_diff(),nel); 
                
            if(infosw>=2){
                store_numdens(nel, acc_lep.get_p(),acc_lep.get_gamma(),acc_lep.get_pdens(),acc_lep.get_gdens(), output.numdens);
            }
        } else if (zone.nth_frac < 1.) {
            if (IsShock==false){
                t_e = f_heat*t_e;
                zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
                IsShock = true;
            }
            Thermal dummy_elec(nel);
            dummy_elec.set_temp_kev(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pbrk = dummy_elec.av_p();
            
            Bknpower acc_lep(nel);
            acc_lep.set_cutoff_type(bhjet.get_cutoff_type());
            acc_lep.set_pspec1(-2.);
            acc_lep.set_pspec2(pspec);
                        
            if (f_sc<10.){
                acc_lep.set_p(0.1*pbrk,pbrk,Urad,zone.bfield,f_beta,zone.r,f_sc);
            } else{
                acc_lep.set_p(0.1*pbrk,pbrk,f_sc);
            }	
            
            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,f_beta);
            //Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel-1];

            zone.avgammasq = pow(acc_lep.av_gamma(),2.);

            gsl_spline_init(spline_eldis,acc_lep.get_gamma(),acc_lep.get_gdens(),nel);
            gsl_spline_init(spline_deriv,acc_lep.get_gamma(),acc_lep.get_gdens_diff(),nel); 
                
            if(infosw>=2){
                store_numdens(nel, acc_lep.get_p(),acc_lep.get_gamma(),acc_lep.get_pdens(),acc_lep.get_gdens(), output.numdens);
            }
        } else if (zone.nth_frac == 1.) {
            if (IsShock==false){
                t_e = f_heat*t_e;
                zone.eltemp = max(tshift*t_e*pow(log10(z_diss)/log10(z),f_pl),1.);
                IsShock = true;
            }
            Thermal dummy_elec(nel);
            dummy_elec.set_temp_kev(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pmin = dummy_elec.av_p();
            
            Powerlaw acc_lep(nel);
            acc_lep.set_cutoff_type(bhjet.get_cutoff_type());
            acc_lep.set_pspec(pspec);
                        
            if (f_sc<10.){
                acc_lep.set_p(pmin,Urad,zone.bfield,f_beta,zone.r,f_sc);
            } else{
                acc_lep.set_p(pmin,f_sc);
            }	
            
            acc_lep.set_norm(zone.lepdens);	
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad,zone.lepdens,zone.bfield,zone.r,f_beta);
            //Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel-1];

            zone.avgammasq = pow(acc_lep.av_gamma(),2.);

            gsl_spline_init(spline_eldis,acc_lep.get_gamma(),acc_lep.get_gdens(),nel);
            gsl_spline_init(spline_deriv,acc_lep.get_gamma(),acc_lep.get_gdens_diff(),nel); 
             	
            if(infosw>=2){
                store_numdens(nel, acc_lep.get_p(),acc_lep.get_gamma(),acc_lep.get_pdens(),acc_lep.get_gdens(), output.numdens);
            }
        }

        if (infosw>=3){
            double Up,Ue,Ub;
            Ue = sqrt(zone.avgammasq)*zone.lepdens*emerg;
            Up = (zone.lepdens/nozzle_ener.eta)*pmgm*pow(cee,2.);
            Ub = pow(zone.bfield,2.)/(8.*pi);

            output.jet_zone_properties.jet_bfield.push_back(zone.bfield); 
            output.jet_zone_properties.lepton_ndens.push_back(zone.lepdens); 
            output.jet_zone_properties.speed_gamma.push_back(zone.gamma); 
            output.jet_zone_properties.delta.push_back(zone.delta); 
            output.jet_zone_properties.tshift.push_back(tshift);
            output.jet_zone_properties.temp_kev.push_back(zone.eltemp);  
            output.jet_zone_properties.grid_r.push_back(zone.r/Rg); 
            output.jet_zone_properties.delz.push_back(zone.delz/Rg);
            output.jet_zone_properties.dist_z.push_back(z/Rg);  
            output.jet_zone_properties.z_delz.push_back((zone.delz+z)/Rg); 
            output.jet_zone_properties.delta.push_back(zone.delta); 
            output.jet_zone_properties.equpar_check.push_back(2.*Ub/Up);
            output.jet_zone_properties.ue_ub.push_back(Ue/Ub); 

            // saving the profiles: 
            output.jetprofile.z_rg.push_back(z/Rg); 
            output.jetprofile.zone_rg.push_back(zone.r/Rg); 
            output.jetprofile.zone_bfield.push_back(zone.bfield); 
            output.jetprofile.zone_lepdens.push_back(zone.lepdens); 
            output.jetprofile.zone_gamma.push_back(zone.gamma); 
            output.jetprofile.zone_eltemp.push_back(zone.eltemp); 
        }

        //calculate emission of each zone		
        //note: the syn_en array is used for the seed photon fields in the IC part, so it needs to include
        //both the black body and disk part. This is why the maximum frequency is taken as the maximum of the
        //two scale frequencies.

        syn_min = 0.1*pow(gmin,2.)*charg*zone.bfield/(2.*pi*emgm*cee);

        if(r_in<r_out){
            syn_max = max(50.*pow(gmax,2.)*charg*zone.bfield/(2.*pi*emgm*cee),20.*Disk.tin()*kboltz/herg);
        }
        else {
            syn_max = 50.*pow(gmax,2.)*charg*zone.bfield/(2.*pi*emgm*cee);
        }

        nsyn = int(log10(syn_max)-log10(syn_min))*syn_res;
        syn_en = new double[nsyn];			
        syn_lum = new double[nsyn];	
        Cyclosyn Syncro(nsyn);
        Syncro.set_frequency(syn_min,syn_max);	

        com_min = 0.1*Syncro.nu_syn();
        com_max = ear[ne-1]/hkev;
        ncom = int(log10(com_max)-log10(com_min))*com_res;
        com_en = new double[ncom];			
        com_lum = new double[ncom];		
        Compton InvCompton(ncom,nsyn);
        InvCompton.set_frequency(com_min,com_max);

        if(infosw>1) {
            for (int k=0;k<ncom;k++){
                com_lum[k] = 0;
                com_en[k] = InvCompton.get_energy()[k];
            }
        }

        //calculate cyclosynchrotron spectrum
        //Set up the calculation by reading in magnetic field,beaming,volume,counterjet presence	
        Syncro.set_bfield(zone.bfield);
        Syncro.set_beaming(theta,zone.beta,zone.delta);
        Syncro.set_geometry("cylinder",zone.r,zone.delz);
        Syncro.set_counterjet(true);
        Syncro.cycsyn_spectrum(gmin,gmax,spline_eldis,acc_eldis,spline_deriv,acc_deriv);
        sum_counterjet(nsyn,Syncro.get_energy_obs(),Syncro.get_nphot_obs(),syn_en,syn_lum);	        

        if (infosw>=4){
            Syncro.test();
        }  

        //Include zone's emission to the pre/post particle acceleration spectrum
        if(z<z_diss){
            sum_zones(nsyn,ne,syn_en,syn_lum,tot_en,tot_syn_pre);
        } else {
            sum_zones(nsyn,ne,syn_en,syn_lum,tot_en,tot_syn_post);
        }

        //calculate inverse Compton spectrum, if it's expected to be bright enough	
        if (Compton_check(IsShock,i,Mbh,jetrat,Urad,velsw,zone) == true){

            InvCompton.set_beaming(theta,zone.beta,zone.delta);
            InvCompton.set_geometry("cylinder",zone.r,zone.delz);
            InvCompton.set_counterjet(true);	
            InvCompton.set_tau(zone.lepdens,zone.eltemp);

            if(InvCompton.get_ypar() > 1.e-2 && InvCompton.get_tau() > 5.e-2){
                InvCompton.set_niter(15);		
            }

            //Cyclosynchrotron photons are always considered in the scattering						
            InvCompton.cyclosyn_seed(Syncro.get_energy(),Syncro.get_nphot());
            
            //Disk photons are included only if the disk is present
            if(r_in<r_out){
                InvCompton.shsdisk_seed(Syncro.get_energy(),Disk.tin(),r_in,r_out,Disk.hdisk(),z+zone.delz/2.);
            }
            //Black body photons included only if compsw==1
            if(compsw==1){						
                InvCompton.bb_seed_k(Syncro.get_energy(),Ubb1,zone.delta*BlackBody.temp_k());
            }
            //AGN photon fields photons are considered only if disk is present and compsw==2
            if(compsw==2 && r_in<r_out){
                InvCompton.bb_seed_k(Syncro.get_energy(),Ubb1,zone.delta*BLR.temp_k());
                InvCompton.bb_seed_k(Syncro.get_energy(),Ubb2,zone.delta*Torus.temp_k());
            }

            // Calculate the spectrum with whichever fields have been invoked		
            InvCompton.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
            sum_counterjet(ncom,InvCompton.get_energy_obs(),InvCompton.get_nphot_obs(),com_en,com_lum);  

            if (infosw>=4){
                InvCompton.test();
            }

            if(z<z_diss){
                sum_zones(ncom,ne,com_en,com_lum,tot_en,tot_com_pre);
            } else {
                sum_zones(ncom,ne,com_en,com_lum,tot_en,tot_com_post);
            }
        } else if (infosw>=5){
            cout << "Out of the Comptonization region" << endl;
        }
    
        if(infosw>=2){
            store_output(nsyn,syn_en,syn_lum,output.cyclosyn_zones, dist, redsh); 
            store_output(ncom,com_en,com_lum, output.compton_zones, dist,redsh);
        }			
        delete[] syn_en,delete[] syn_lum;
        delete[] com_en,delete[] com_lum;

    } //closing bracket for jet loop 


//FINAL STEP: SUM JET COMPONENTS TO TOTAL OUTPUT, WRITE/CLOSE PLOT FILES, FREE MEMORY

    for(int k=0;k<ne;k++){
        tot_lum[k] = (tot_lum[k]+tot_syn_pre[k]+tot_syn_post[k]+tot_com_pre[k]+tot_com_post[k]);
        photeng[k] = log10(tot_en[k]/herg);
    }
    
    if(infosw==0){
        cout << "Infosw is turned off, nothing will be written to the output directories and arrays" << endl;
    }

    if(infosw>=1){
        
        store_output(ne,tot_en,tot_syn_pre,output.presyn,dist,redsh); 
        store_output(ne,tot_en,tot_syn_post,output.postsyn,dist,redsh); 
        store_output(ne,tot_en,tot_com_pre,output.precom,dist,redsh); 
        store_output(ne,tot_en,tot_com_post,output.postcom,dist,redsh);
        store_output(50,Disk.get_energy_obs(),Disk.get_nphot_obs(),output.disk,dist,redsh); 

        if(compsw==2){
            store_output(40,Torus.get_energy_obs(),Torus.get_nphot_obs(),output.bb,dist,redsh); 
        } else {
            store_output(40,BlackBody.get_energy_obs(),BlackBody.get_nphot_obs(),output.bb,dist,redsh); 
        }			
        store_output(ne,tot_en,tot_lum,output.total,dist,redsh); 
    }

    if (infosw >=3) {
        double disk_lum,IC_lum,Xray_lum,Radio_lum,Xray_index,Radio_index,compactness;
        disk_lum = integrate_lum(50,0.3*2.41e17,5.*2.41e17,Disk.get_energy_obs(),Disk.get_nphot_obs());
        IC_lum = integrate_lum(ne,0.3*2.41e17,300.*2.41e17,tot_en,tot_com_pre);
        Xray_lum = integrate_lum(ne,1.*2.41e17,10.*2.41e17,tot_en,tot_lum);
        Radio_lum = integrate_lum(ne,4e9,6e9,tot_en,tot_lum);
        Xray_index = photon_index(ne,10.*2.41e17,100.*2.41e17,tot_en,tot_lum);
        Radio_index = 1.+photon_index(ne,1e10,1e11,tot_en,tot_lum);
        compactness = integrate_lum(ne,0.1*2.41e17,300.*2.41e17,tot_en,tot_com_pre)*sigtom/(r_0*emerg*cee);   

        output.spectral_properties.disk_lum.push_back(disk_lum); 
        output.spectral_properties.IC_lum.push_back(IC_lum); 
        output.spectral_properties.xray_lum.push_back(Xray_lum); 
        output.spectral_properties.radio_lum.push_back(Radio_lum); 
        output.spectral_properties.xray_index.push_back(Xray_index); 
        output.spectral_properties.radio_index.push_back(Radio_index); 
        output.spectral_properties.jetbase_compactness.push_back(compactness); 
    }

//Warning for compactness --- 
    // if (compactness >= 10.*(param[9]/511.)*exp(511./param[9])) {
    //         cout << "Possible pair production in the jet base!" << endl; 
    //         cout << "Lower limit on allowed compactness: " << 10.*(param[9]/511.)*exp(511./param[9]) << endl;
    //         cout << "Note: this is for a slab, a cylinder allows higher l by a factor of ~10" << std::endl;}
    // }	

    gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);
    gsl_spline_free(spline_deriv), gsl_interp_accel_free(acc_deriv);	
    gsl_spline_free(spline_speed), gsl_interp_accel_free(acc_speed);
    delete[] tot_en,delete[] tot_lum;
    delete[] tot_syn_pre,delete[] tot_syn_post;
    delete[] tot_com_pre,delete[] tot_com_post;

}




void singlezone_jetmain(BhJetClass& bhjet, JetOutput& output) {


    int nel = 100;			    //array size for particle distributions
    int nfreq = 200;            //array size for frequency arrays
    
    double B,R,n;               //plasma quantities of emitting region: bfield, radius, number density
    double gmin,gmax,p;         //details of particle distribution: minimum/maximum gamma factor, powerlaw slope
    double pmin;                //minimum momentum corresponding to gmin
    double delta,gamma,beta;    //plasma speed/beaming factor
    double theta;               //jet viewing angle
    double Pj;                  //total power carried in the jet 
    double Pe, Ue;              //power/energy density in electrons
    double Pb, Ub;              //power/energy density in bfields
    double Pp, Up;              //power/energy density in cold protons
    double equip;               //standard equipartition factor Ue/Ub
    double nus_min,nus_max;     //synchrotron frequency range
    double nuc_min,nuc_max;     //SSC frequency range 

    //splines for electron distribution. These are needed by the radiation codes
    gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
    gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen,nel);

    gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
    gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen,nel);  


    double Mbh = bhjet.Mbh;
    double Eddlum = 1.25e38 * Mbh;
    double Rg = gconst * Mbh * msun / (cee * cee);

	B = 1.5e-3;
	R = 626.*Rg;
	n = 9.5e-3;
	gmin = 4.1e3;
	pmin = pow(pow(gmin,2.)-1.,1./2.)*emgm*cee;
	gmax = 6.4e7;
	p = 3.03;
	delta = 3.3;
	gamma = 3.;
	beta = 0.942809;
	theta = acos((delta*gamma-1.)/(beta*delta*gamma))*180./pi;

    nus_min = 1.e8;
	nus_max = 1.e21;
	nuc_min = 1.e17;
	nuc_max = 1.e28;

    Powerlaw Electrons(nel);
    Electrons.set_p(pmin,gmax); 
    Electrons.set_pspec(p);
    Electrons.set_norm(n);	
    Electrons.set_ndens();

    store_numdens(nel, Electrons.get_p(),Electrons.get_gamma(),Electrons.get_pdens(),Electrons.get_gdens(), output.numdens);

    gsl_spline_init(spline_eldis,Electrons.get_gamma(),Electrons.get_gdens(),nel);
    gsl_spline_init(spline_deriv,Electrons.get_gamma(),Electrons.get_gdens_diff(),nel);

    Cyclosyn Syncro(nfreq);
    Syncro.set_frequency(nus_min,nus_max);    
    Syncro.set_bfield(B);
    Syncro.set_beaming(theta,beta,delta);
    Syncro.set_geometry("sphere",R);
    Syncro.cycsyn_spectrum(gmin,gmax,spline_eldis,acc_eldis,spline_deriv,acc_deriv);
    store_output(nfreq,Syncro.get_energy_obs(),Syncro.get_nphot_obs(), output.cyclosyn_zones, 1., 0.); 

    Compton InvCompton(nfreq,nfreq);
    InvCompton.set_frequency(nuc_min,nuc_max);	    
    InvCompton.set_beaming(theta,beta,delta);
    InvCompton.set_geometry("sphere",R);	
    InvCompton.set_tau(n,Electrons.av_gamma()*511.);			
    InvCompton.cyclosyn_seed(Syncro.get_energy(),Syncro.get_nphot());
    InvCompton.compton_spectrum(gmin,gmax,spline_eldis,acc_eldis);
    store_output(nfreq,InvCompton.get_energy_obs(),InvCompton.get_nphot_obs(), output.compton_zones, 1., 0.);
    
    gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);

}