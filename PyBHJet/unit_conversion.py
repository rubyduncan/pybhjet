CC = 3*10**10
PLANCK_CONSTANT_J_S = 6.626e-34  # Planck's constant in J·s
PLANCK_CONSTANT_ERG_S = 6.626 * 10 **-27 #ergs 
ELECTRONVOLT_TO_JOULE = 1.602e-19  # conversion from eV to Joules
TEV_TO_JOULE = ELECTRONVOLT_TO_JOULE * 1e12  # conversion from TeV to Joules

ONE_KEV_in_ERG = 1.602*10**-9 #1 keV = 10^-9 ergs
MJY_TO_CGS  = 1e-26  # mJy -> erg/cm^2/s/Hz

def hz_to_kev(nu_hz): 
    #e = h (keV/s) * nu (/hz) = kev
    e_keV = nu_hz *(4.135 * 10 **-18) #keV/s
    return e_keV

def kev_to_hz(keV): 
    #freq = E/h
    nu_hz = (2.42*10**17) * keV 
    return nu_hz

def mjy_to_diff_photon_flux(flux_mjy, freq):
    
    '''
    expects flux density in mjy, 
    returns dNp/dA dt dE 
    1) energy flux / unit frequency: conv F_nu [erg/cm2/s/hz] = S_nu * 1e-26 [mJy]
    2) energy flux/ energy: E * N_e(E) (which means):
        a) F_e(E) = E_kev * (one kev in erg) * N_e(E) (relation)
        b) N_e(E) = F_e(E)/(E_kev * one kev in erg) *** this is the unit
    3) energy flux: F_nu dnu = F_e dE; ***> F_e = F_nu dnu/dE
        a) dnu/dE = (nu/E (erg) = 1/h (erg)) *bc it starts in erg  
        b) energy flux F_e = F_nu(E) [erg/cm2/s/hz] * / 1/h [1/erg s] 
            i) unit check: erg/cm2/s/hz * 1/erg s = erg/cm2/ * 1/erg s = cm2/s  
            ii) F_e = [cm2/s] **good (technically erg/cm2/s/erg)
    4) F_e[keV] = F_nu * (1/h) * (1 kev in erg) = erg/cm2/s/keV
    5) N_e(E) = F_e[keV] / (E_kev * one kev in erg)

    '''
    
    kev = hz_to_kev(freq)

    #step 0:
    f_nu = flux_mjy * MJY_TO_CGS #erg/cm2/s/hz

    #step 3: 
    energy_flux = f_nu * 1/PLANCK_CONSTANT_ERG_S #erg/cm2/s/erg 

    #step 4: 
    energy_flux_h = energy_flux * ONE_KEV_in_ERG

    #step 5: 
    number_flux_kev = energy_flux_h / (kev * ONE_KEV_in_ERG)

    return number_flux_kev 

def photon_flux_density_to_mjy(number_flux, kev):
    ''' 
    Expects number flux, and energy in keV
    returns flux in mJy 
    
    '''
    energy_flux_h = number_flux * (kev * ONE_KEV_in_ERG)
    energy_flux = energy_flux_h/ONE_KEV_in_ERG

    f_nu = energy_flux * PLANCK_CONSTANT_ERG_S

    return f_nu / MJY_TO_CGS
