import numpy as np
from astromodels.functions.function import ModelAssertionViolation

import astropy.units as u
from astromodels.functions.function import (
    Function1D,
    FunctionMeta,
    ModelAssertionViolation,
)

from bhjet_plotting import (
    preprocess_component_output,
    preprocess_radiative_zones_output,
    preprocess_numdens_output,
    preprocess_jet_profile,
    preprocess_jet_zone_properties,
    preprocess_jet_base_properties,
    preprocess_spectral_properties,
)

# this needs to be compatible with the location of the library with the model 
import build.pybhjet as pybhjet

# to define a custom model in 3ml, need to include docstring, units setter, evaluate function: 
class BHJetModel(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        BHJet: steady state, multi-zone jet model 
    latex : $ tbd $
    parameters : 
        Mbh :
            desc : Black hole mass, in solar masses 
            initial value : 6.6e8
            min : 1e5
            max : 1e9
            delta : 0.1
        theta : 
            desc : inclination of the jet, line of sight
            initial value : 65
            min : 1
            max : 85
            fix : yes
        dist : 
            desc : distance to the source, in kpc
            initial value : 9000
            fix : yes
        redsh : 
            desc : source redshift 
            initial value : 0.003
            fix : yes 
        jetrat : 
            desc : injected jet power, in Edd. Lum 
            initial value : 1e-2
            min : 1e-7
            max : 1e0
            delta : 0.1
        r_0 : 
            desc : jet base radius, in rg 
            initial value : 10
            min : 2
            max : 50
            delta : 0.1
        z_diss : 
            desc : dissipation region distance along the jet, in rg 
            initial value : 50
            min : 30
            max : 500 
            delta : 0.1
        z_acc : 
            desc : acceleration region distance along the jet, in rg 
            initial value : 50 
            min : 30
            max : 500
            delta : 0.1
        z_max : 
            desc : maximum distance of the jet, in rg
            initial value : 100000 
            fix : yes 
        t_e : 
            desc : temperature of rel. e in the nozzle, keV
            initial value : 400 
            min : 50
            max : 1200 
            delta : 0.1
        f_nth : 
            desc : fraction of particles accelerated into nT tail 
            initial value : 0.1
            min : 0.05
            max : 0.9
            delta : 0.1
        f_pl : 
            desc : dissipation parameter for particles along the jet 
            initial value : 0
            min : 0
            max : 10 
        pspec : 
            desc : slope of nT lepton distribution 
            initial value : 2
            min : 1.5
            max : 3
            delta : 0.1
        f_heat : 
            desc : factor to increase electron temp. 
            initial value : 1 
            delta : 0.1 
        f_beta : 
            desc : adiabatic cooling timescale
            initial value : 0.1
            delta : 0.1
        f_sc : 
            desc : maximum energy of nT particles
            initial value : 2e-9
            min : 1e-9
            max : 1
            delta : 0.1
        p_beta : 
            desc : plasma beta 
            initial value : 0.02
            min : 0.001
            max : 0.04 
            delta : 0.1
        sig_acc : 
            desc : Magnetization at acceleration region
            initial value : 0.01
            min : 0.01
            max : 1
            delta : 0.1
        l_disk : 
            desc : Disk Luminosity, L_edd
            initial value : 1e-5
            min : 1e-7
            max : 1
        r_in : 
            desc : inner disk radius, rg 
            initial value : 20
            min : 10
            max : 3000 
        r_out : 
            desc : outer disk radius, rg
            initial value : 1000
            min : 500
            max : 1500
        compar1 : 
            desc : Compsw=1 Temp of BB [K], Compsw=2 frac of disk photons reprocessed by BLR
            initial value : 1500
            min : 1000
            max : 3000
        compar2 : 
            desc : Compsw=1 Temp of BB [K], Compsw=2 frac of disk photons reprocessed by BLR
            initial value : 5e40
            min : 1e35
            max : 1e43
        compar3 : 
            desc : Compsw=1 BB Energy Density
            initial value : 1e-8
            min : 1e-10
            max : 1e-1
        compsw : 
            desc : Adds Ext. Photon Field, = 1 BB, = 2 BLR & Torus
            initial value : 0
            min : 0
            max : 3
            fix : yes 
        velsw : 
            desc : =0, 1 ->  agnjet(pressure driven), > 1 --> bljet(magnetic driven)
            initial value : 4
            min : 0
            max : 25
            fix : yes 
        infosw : 
            desc : returns information about the code
            initial value : 1
            fix : yes 
        EBLsw : 
            desc : account for ebl 
            initial value : 0
            fix : yes

    """
    cutoff_type = 0

    def _setup(self):
        self.bhjet = pybhjet.PyBHJet()
        self._cached_params   = None
        self._cached_E_keV    = None
        self._cached_flux_ph  = None

        #code switches: 
        self.bhjet.cutoff_type = int(self.cutoff_type)

        # plotting/components cache
        self._last_output              = None
        self._last_components          = None
        self._last_radiative_zones     = None
        self._last_numdens             = None
        self._last_jet_profile         = None
        self._last_jet_zone_properties = None
        self._last_jet_base_properties = None
        self._last_spectral_properties = None

        # toggle to avoid doing heavy diag during big MCMC runs
        self.enable_detailed_output = False
  

    def _set_units(self, x_unit, y_unit):
        # Units for input energy grid (e.g., keV)
        # self.x_unit = x_unit

        # Output unit: Photon flux (ph / cm^2 / s / keV)
        # self.y_unit = y_unit

        #--------- units for parameters? -----
        self.Mbh.unit = u.dimensionless_unscaled
        self.theta.unit = u.dimensionless_unscaled
        self.dist.unit = u.dimensionless_unscaled
        self.redsh.unit = u.dimensionless_unscaled
        self.jetrat.unit = u.dimensionless_unscaled
        self.r_0.unit = u.dimensionless_unscaled
        self.z_diss.unit = u.dimensionless_unscaled
        self.z_acc.unit = u.dimensionless_unscaled
        self.z_max.unit = u.dimensionless_unscaled
        self.t_e.unit = u.dimensionless_unscaled
        self.f_nth.unit = u.dimensionless_unscaled
        self.f_pl.unit = u.dimensionless_unscaled
        self.pspec.unit = u.dimensionless_unscaled
        self.f_heat.unit = u.dimensionless_unscaled
        self.f_beta.unit = u.dimensionless_unscaled
        self.f_sc.unit = u.dimensionless_unscaled
        self.p_beta.unit = u.dimensionless_unscaled
        self.sig_acc.unit = u.dimensionless_unscaled
        self.l_disk.unit = u.dimensionless_unscaled
        self.r_in.unit = u.dimensionless_unscaled
        self.r_out.unit = u.dimensionless_unscaled
        self.compar1.unit = u.dimensionless_unscaled
        self.compar2.unit = u.dimensionless_unscaled
        self.compar3.unit = u.dimensionless_unscaled
        self.compsw.unit = u.dimensionless_unscaled
        self.velsw.unit = u.dimensionless_unscaled
        self.infosw.unit = u.dimensionless_unscaled
        self.EBLsw.unit = u.dimensionless_unscaled


    def evaluate(self, x, Mbh, theta,  dist, redsh, jetrat, r_0, z_diss, z_acc, z_max, t_e, 
                 f_nth, f_pl, pspec, f_heat, f_beta, f_sc, p_beta, sig_acc, l_disk, r_in, r_out, 
                 compar1, compar2, compar3, compsw, velsw,infosw, EBLsw
                 ):
        
        params_vec = (
            float(Mbh), float(theta), float(dist), float(redsh),
            float(jetrat), float(r_0), float(z_diss), float(z_acc), float(z_max),
            float(t_e), float(f_nth), float(f_pl), float(pspec), float(f_heat),
            float(f_beta), float(f_sc), float(p_beta), float(sig_acc),
            float(l_disk), float(r_in), float(r_out),
            float(compar1), float(compar2), float(compar3),
            float(compsw), float(velsw), float(infosw), float(EBLsw),
        )
        
        """
        Map the 3ML parameters to Pybhjet, run the model, and return the interpolated output 
        """

        # only if params changed, call BHJet - this helps a lot with computation time per dataset
        if self._cached_params is None or params_vec != self._cached_params:

            # when jetmain is run (so bhjet.run()), premap parameters to BHJet
            self.bhjet.set_parameter("Mbh", Mbh)
            self.bhjet.set_parameter("theta", theta)
            self.bhjet.set_parameter("dist", dist)
            self.bhjet.set_parameter("redsh", redsh)
            self.bhjet.set_parameter("jetrat", jetrat)
            self.bhjet.set_parameter("r_0", r_0)
            self.bhjet.set_parameter("z_diss", z_diss)
            self.bhjet.set_parameter("z_acc", z_acc)
            self.bhjet.set_parameter("z_max", z_max)
            self.bhjet.set_parameter("t_e", t_e)
            self.bhjet.set_parameter("f_nth", f_nth)
            self.bhjet.set_parameter("f_pl", f_pl)
            self.bhjet.set_parameter("pspec", pspec)
            self.bhjet.set_parameter("f_heat", f_heat)
            self.bhjet.set_parameter("f_beta", f_beta)
            self.bhjet.set_parameter("f_sc", f_sc)
            self.bhjet.set_parameter("p_beta", p_beta)
            self.bhjet.set_parameter("sig_acc", sig_acc)
            self.bhjet.set_parameter("l_disk", l_disk)
            self.bhjet.set_parameter("r_in", r_in)
            self.bhjet.set_parameter("r_out", r_out)
            self.bhjet.set_parameter("compar1", compar1)
            self.bhjet.set_parameter("compar2", compar2)
            self.bhjet.set_parameter("compar3", compar3)
            self.bhjet.set_parameter("compsw", compsw)
            self.bhjet.set_parameter("velsw", velsw)
            self.bhjet.set_parameter("infosw", infosw)
            self.bhjet.set_parameter("EBLsw", EBLsw)
            
            self.bhjet.cutoff_type = int(self.cutoff_type)
            self.bhjet.run()

            out = self.bhjet.get_output() #returning linear arrays for freq, flux 

            if self.enable_detailed_output==True:
                self._last_output              = out
                self._last_components          = preprocess_component_output(out)
                self._last_radiative_zones     = preprocess_radiative_zones_output(out)
                self._last_numdens             = preprocess_numdens_output(out)
                self._last_jet_profile         = preprocess_jet_profile(out)
                self._last_jet_zone_properties = preprocess_jet_zone_properties(out)
                self._last_jet_base_properties = preprocess_jet_base_properties(out)
                self._last_spectral_properties = preprocess_spectral_properties(out)

    # interpolation from BHJet energy grid to 3ml x points : 
            native_energy = np.array([point.energy for point in out.total])  # Hz
            native_flux_mjy = np.array([point.flux for point in out.total])  # mJy 

        # convert energy to expected keV output - tries to catch for different cases of incorrect returns
            if native_energy.size == 0:
                self._cached_E_keV   = None
                self._cached_flux_ph = None
            else:
                E_keV = native_energy * 4.135667696e-18

                h_ergs_s   = 6.62607015e-27
                mjy_to_cgs = 1e-26
                flux_ph = (mjy_to_cgs * native_flux_mjy / h_ergs_s) / E_keV

                mask = ( #trying to not kill threeml 
                    (E_keV > 0)
                    & np.isfinite(E_keV)
                    & np.isfinite(flux_ph)
                )
                E_keV   = E_keV[mask]
                flux_ph = flux_ph[mask]

                if E_keV.size < 2:
                    self._cached_E_keV   = None
                    self._cached_flux_ph = None
                else:
                    order = np.argsort(E_keV)
                    self._cached_E_keV   = E_keV[order]
                    self._cached_flux_ph = flux_ph[order]

            self._cached_params = params_vec

        x = np.asarray(x, dtype=float)
        y = np.zeros_like(x, dtype=float) + 1e-99 #this was an issue before 


        if self._cached_E_keV is None or self._cached_E_keV.size < 2:
            return y

        Egrid = self._cached_E_keV
        Fgrid = self._cached_flux_ph

        inside = (
            np.isfinite(x) & (x > 0) & (x >= Egrid[0]) & (x <= Egrid[-1])
        )
        if np.any(inside):
            y[inside] = np.interp(x[inside], Egrid, Fgrid)

        bad = ~np.isfinite(y) | (y < 0)
        y[bad] = 1e-99

        return y

        