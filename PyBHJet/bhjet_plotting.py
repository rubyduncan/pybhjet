import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math 
from matplotlib import rc, rcParams
from matplotlib.colors import LogNorm, Normalize
from plotting import DEFAULT_STYLE, plot_bhjet_component_data

#for using latex commands for plotting
# rc('text', usetex=True)
# rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif Display']})
# plt.rcParams.update({'font.size': 14})

kev_conv = 2.41*10**17
mjy_conv = 1.e26

def flux_conv(dist):    
    return 4.*math.pi*(dist*3.*10**21)**2


def preprocess_component_output(output, comp=None): 
    """
    Extract and store data from output for specified components.
    For all emission components out of bhjet, they are saved in structs of energy & flux, with units [hz, mjy].
    """
    data = {}
    all_components = ["disk", "presyn", "postsyn", "precom", "postcom", "bb", "total"]
    components_to_extract = [comp] if comp is not None else all_components

    for component in components_to_extract:
        try:
            energy = np.array([point.energy for point in getattr(output, component)])
            flux = np.array([point.flux for point in getattr(output, component)])
            if len(energy) > 0 and len(flux) > 0:
                data[component] = {"energy": energy, "flux": flux}
            else:
                print(f"No data found for component: {component}")
        except AttributeError:
            print(f"Component {component} not found in output.")
    return data


def preprocess_radiative_zones_output(output):

    ''' Returns the flux and frequency for the spectral contribution from 
    cyclo_synchrotron and compton emission from each zone. (#nu [Hz]:  Flux [mJy]:)
    '''
    data = {}
    all_components = ['cyclosyn_zones', "compton_zones"]

    for component in all_components:
        try: 
            energy = np.array([point.energy for point in getattr(output, component)])
            flux = np.array([point.flux for point in getattr(output, component)])
            data[component] = {"energy": energy, "flux": flux}

        except AttributeError:
            print(f"Component {component} not found -> need to turn up infosw to 2.") 

    return data 
            
    


def preprocess_numdens_output(output): 
    """
    Extract number density data from output.numdens (vector of NumDenPoint).
    """
    data = {}
    try:
        numdens = output.numdens

        momentum = np.array([pt.momentum for pt in numdens])
        gamma    = np.array([pt.gamma    for pt in numdens])
        n_p      = np.array([pt.n_p      for pt in numdens])
        n_g      = np.array([pt.n_g      for pt in numdens])

        data = {
            "momentum": momentum,
            "gamma": gamma,
            "n_p": n_p,
            "n_g": n_g
        }
        
    except AttributeError:
        print("[ERROR] Could not access output.numdens or its fields.")
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")

    return data


def _detail_column(detail, name):
    """Return one named column from a detail dictionary or DataFrame."""
    if detail is None:
        return None
    if isinstance(detail, dict):
        value = detail.get(name)
    else:
        try:
            value = detail[name]
        except (KeyError, TypeError):
            return None
    return None if value is None else np.asarray(value, dtype=float)


def plot_lepton_distribution_per_zone(
    numdens,
    jet_profile=None,
    *,
    distribution="momentum",
    weighted=True,
    zones=None,
    every_nth=1,
    points_per_zone=None,
    cmap="magma",
    ax=None,
    colorbar=True,
    linewidth=1.25,
    alpha=0.85,
    title=None,
):
    """Plot BHJet lepton distributions for the individual jet zones.

    Parameters
    ----------
    numdens : dict
        The mapping returned by ``jet.get_detail("numdens", infosw=3)``.
        BHJet stores each zone consecutively in the flattened arrays.
    jet_profile : dict or pandas.DataFrame, optional
        The matching result from ``jet.get_detail("jet_profile", infosw=3)``.
        When supplied, curves are coloured by their distance along the jet,
        ``z_rg``. Otherwise they are coloured by zone index.
    distribution : {"gamma", "momentum"}
        Plot ``n(gamma)`` or ``n(p)``. The default, ``"momentum"``, matches
        the per-zone particle plots previously used in the BHJet notebooks.
    weighted : bool, default True
        Plot ``gamma * n(gamma)`` or ``p * n(p)``. Set False for the
        differential number-density distribution itself.
    zones : iterable of int, optional
        Explicit zone indices to show. By default all zones are eligible.
    every_nth : int, default 1
        Keep every Nth eligible zone, useful for dense detailed output.
    points_per_zone : int, optional
        Override the number of particle-grid samples per zone. This is
        inferred from ``jet_profile`` when possible; otherwise BHJet's
        standard 70-point particle grid is used.

    Returns
    -------
    (matplotlib.figure.Figure, matplotlib.axes.Axes)
        The figure and axis containing the curves.
    """
    if distribution not in {"gamma", "momentum"}:
        raise ValueError("distribution must be 'gamma' or 'momentum'")
    if every_nth < 1:
        raise ValueError("every_nth must be at least 1")

    x_key, density_key = ("gamma", "n_g") if distribution == "gamma" else ("momentum", "n_p")
    x_values = _detail_column(numdens, x_key)
    density = _detail_column(numdens, density_key)
    if x_values is None or density is None:
        raise ValueError(f"numdens must contain '{x_key}' and '{density_key}' arrays")
    if x_values.ndim != 1 or density.ndim != 1 or len(x_values) != len(density):
        raise ValueError("numdens coordinate and density arrays must be one-dimensional and equally sized")

    z_rg = _detail_column(jet_profile, "z_rg")
    n_zones = len(z_rg) if z_rg is not None and len(z_rg) else None
    if points_per_zone is None:
        if n_zones is not None and len(x_values) % n_zones == 0:
            points_per_zone = len(x_values) // n_zones
        elif len(x_values) % 70 == 0:
            points_per_zone = 70  # BHJet's C++ particle grid uses nel = 70.
        else:
            raise ValueError(
                "Cannot infer the number of particle samples per zone. Pass "
                "points_per_zone explicitly or provide the matching jet_profile."
            )
    points_per_zone = int(points_per_zone)
    if points_per_zone < 2 or len(x_values) % points_per_zone:
        raise ValueError("points_per_zone must divide the numdens arrays and be at least 2")

    inferred_zones = len(x_values) // points_per_zone
    if n_zones is not None and n_zones != inferred_zones:
        raise ValueError(
            f"jet_profile has {n_zones} zones but numdens contains {inferred_zones}; "
            "obtain both details from the same BHJet evaluation."
        )
    if z_rg is None or len(z_rg) != inferred_zones:
        z_rg = np.arange(inferred_zones, dtype=float)
        colour_label = "Zone index"
        norm = Normalize(vmin=0, vmax=max(inferred_zones - 1, 1))
    else:
        valid_z = z_rg[np.isfinite(z_rg) & (z_rg > 0)]
        if len(valid_z) < 2 or np.isclose(valid_z.min(), valid_z.max()):
            norm = Normalize(vmin=0, vmax=max(inferred_zones - 1, 1))
            colour_label = "Zone index"
        else:
            norm = LogNorm(vmin=valid_z.min(), vmax=valid_z.max())
            colour_label = r"Zone location $z/r_g$"

    if zones is None:
        selected_zones = list(range(inferred_zones))
    else:
        selected_zones = [int(zone) for zone in zones]
        invalid = [zone for zone in selected_zones if zone < 0 or zone >= inferred_zones]
        if invalid:
            raise IndexError(f"Zone indices out of range: {invalid}")
    selected_zones = selected_zones[::every_nth]

    if ax is None:
        figure, ax = plt.subplots(figsize=(8, 5.5))
    else:
        figure = ax.figure
    cmap_object = plt.colormaps[cmap]
    plotted = []
    for zone in selected_zones:
        start, stop = zone * points_per_zone, (zone + 1) * points_per_zone
        x_zone = x_values[start:stop]
        density_zone = density[start:stop]
        valid = np.isfinite(x_zone) & np.isfinite(density_zone) & (x_zone > 0) & (density_zone > 0)
        if valid.sum() < 2:
            continue
        x_zone, density_zone = x_zone[valid], density_zone[valid]
        order = np.argsort(x_zone)
        x_zone, density_zone = x_zone[order], density_zone[order]
        y_zone = density_zone * x_zone if weighted else density_zone
        colour_value = z_rg[zone] if colour_label.startswith("Zone location") else zone
        ax.plot(x_zone, y_zone, color=cmap_object(norm(colour_value)), linewidth=linewidth, alpha=alpha)
        plotted.append(zone)

    if not plotted:
        raise ValueError("No positive, finite lepton-distribution values were available for the selected zones")

    variable_label = r"$\gamma$" if distribution == "gamma" else r"$p$ (g cm s$^{-1}$)"
    if distribution == "gamma":
        density_label = r"$\gamma\,n(\gamma)$ (cm$^{-3}$)" if weighted else r"$n(\gamma)$ (cm$^{-3}$)"
    else:
        density_label = r"$p\,n(p)$ (cm$^{-3}$)" if weighted else r"$n(p)$ (cm$^{-3}$ (g cm s$^{-1}$)$^{-1}$)"
    ax.set(xscale="log", yscale="log", xlabel=variable_label, ylabel=density_label)
    ax.grid(which="both", alpha=0.2)
    ax.set_title(title or f"Lepton distribution per BHJet zone ({distribution})")

    if colorbar:
        scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap_object)
        scalar_map.set_array([])
        colour_bar = figure.colorbar(scalar_map, ax=ax, pad=0.02)
        colour_bar.set_label(colour_label)
        if len(plotted) <= 8:
            ticks = [z_rg[zone] if colour_label.startswith("Zone location") else zone for zone in plotted]
            colour_bar.set_ticks(ticks)
    return figure, ax


def plot_lepton_distribution_from_jet(jet, *, e_min_keV=1e-9, e_max_keV=1e3,
                                      n_eval=2, infosw=3, **kwargs):
    """Evaluate a BHJet model once, then plot its lepton distributions by zone.

    This is a convenience wrapper around :func:`plot_lepton_distribution_per_zone`.
    The underlying ``numdens`` and jet-profile details remain available on the
    model as ``_last_numdens`` and ``_last_jet_profile`` after the call.
    """
    numdens = jet.get_detail(
        "numdens", e_min_keV=e_min_keV, e_max_keV=e_max_keV,
        n_eval=n_eval, infosw=infosw,
    )
    return plot_lepton_distribution_per_zone(
        numdens, jet_profile=getattr(jet, "_last_jet_profile", None), **kwargs
    )


def preprocess_jet_profile(output, include_descriptions=False):
    """
    Extract JetProfile data into a Pandas DataFrame.
    Each row corresponds to a zone, and columns represent the different information in the jet profile.
    Can set include_descriptions = True for what is included in jet profile information. 
    """
    try:
        # Extract data from output
        z_rg = np.array(output.jetprofile.z_rg)
        zone_rg = np.array(output.jetprofile.zone_rg)
        zone_bfield = np.array(output.jetprofile.zone_bfield)
        zone_lepdens = np.array(output.jetprofile.zone_lepdens)
        zone_gamma = np.array(output.jetprofile.zone_gamma)
        zone_eltemp = np.array(output.jetprofile.zone_eltemp)
        
        # Put data into a DataFrame
        jet_profile_df = pd.DataFrame({
            "z_rg": z_rg,
            "zone_rg": zone_rg, #zone radius
            "zone_bfield": zone_bfield,
            "zone_lepdens": zone_lepdens,
            "zone_gamma": zone_gamma,
            "zone_eltemp": zone_eltemp
        })
        
        # Print descriptions if requested
        if include_descriptions:
            description_text = (
                "Jet Profile Properties:\n"
                "- z_rg: Distance along the jet in r_g.\n"
                "- zone_rg: Zone radius in r_g.\n"
                "- zone_bfield: Magnetic field in the zone.\n"
                "- zone_lepdens: Lepton number density in the zone.\n"
                "- zone_gamma: Lorentz factor of the zone.\n"
                "- zone_eltemp: Electron temperature in the zone.\n"
            )
            print(description_text)
            
        return jet_profile_df

    except AttributeError as e:
        print(f"Error accessing jet profile data: {e}")
        return None



def preprocess_jet_zone_properties(output, include_descriptions=False):
    """
    Extract jet zone properties into a dictionary and optionally include descriptions.
    
    Args:
        output: The output object from the jet model.
        include_descriptions (bool): If True, print a description of all properties.
        
    Returns:
        Dictionary of jet zone properties.
    """
    jet_zone_properties = {
        "jet_bfield": np.array(output.jet_zone_properties.jet_bfield),
        "lepton_ndens": np.array(output.jet_zone_properties.lepton_ndens),
        "speed_gamma": np.array(output.jet_zone_properties.speed_gamma),
        "delta": np.array(output.jet_zone_properties.delta),
        "tshift": np.array(output.jet_zone_properties.tshift),
        "temp_kev": np.array(output.jet_zone_properties.temp_kev),
        "grid_r": np.array(output.jet_zone_properties.grid_r),
        "delz": np.array(output.jet_zone_properties.delz),
        "dist_z": np.array(output.jet_zone_properties.dist_z),
        "z_delz": np.array(output.jet_zone_properties.z_delz),
        "equpar_check": np.array(output.jet_zone_properties.equpar_check),
        "ue_ub": np.array(output.jet_zone_properties.ue_ub),
    }
    
    if include_descriptions:
        description_text = (
            "The jet zone properties include the following parameters:\n"
            "- 'jet_bfield': Magnetic field strength in the jet zones.\n"
            "- 'lepton_ndens': Number density of leptons in the zones.\n"
            "- 'speed_gamma': Lorentz factor of the jet zones.\n"
            "- 'delta': Doppler factor in each zone.\n"
            "- 'tshift': Time shift between zones due to jet motion.\n"
            "- 'temp_kev': Electron temperature in keV.\n"
            "- 'grid_r': Grid position along the jet radius in cm.\n"
            "- 'delz': Incremental distance along the jet axis.\n"
            "- 'dist_z': Distance from the jet base.\n"
            "- 'z_delz': Sum of current zone's position and delz.\n"
            "- 'equpar_check': Equipartition check value.\n"
            "- 'ue_ub': Ratio of internal energy density to magnetic energy density.\n"
            "\n"
        )
        print(description_text)
    
    return jet_zone_properties
    

def preprocess_jet_base_properties(output, include_descriptions=False):
    """
    Extract jet base properties into a dictionary and optionally include descriptions.
    
    Args:
        output: The output object from the jet model.
        include_descriptions (bool): If True, print a description of all properties.
        
    Returns:
        Dictionary of jet base properties.
    """
    jet_base_properties = {
        "pair_content": np.array(output.jet_base_properties.pair_content),
        "init_mag": np.array(output.jet_base_properties.init_mag),
        "particle_avg_lorentz_factor": np.array(output.jet_base_properties.particle_avg_lorentz_factor),
        "jet_nozzle_end": np.array(output.jet_base_properties.jet_nozzle_end),
        "jet_nozzle_optical_depth": np.array(output.jet_base_properties.jet_nozzle_optical_depth),
    }
    
    if include_descriptions:
        description_text = (
            "The jet base properties include the following parameters:\n"
            "- 'pair_content': Pair content (ratio of pairs to protons) in the jet base.\n"
            "- 'init_mag': Initial magnetic field strength in the jet nozzle.\n"
            "- 'particle_avg_lorentz_factor': Average Lorentz factor of particles in the jet base.\n"
            "- 'jet_nozzle_end': Distance where the jet nozzle ends, in R_g.\n"
            "- 'jet_nozzle_optical_depth': Optical depth in the jet nozzle.\n"
        )
        print(description_text)
    
    return jet_base_properties

def preprocess_spectral_properties(output, include_descriptions=False):
    """
    Extract spectral properties into a dictionary and optionally include descriptions.
    
    Args:
        output: The output of spectral properties from jet model.
        include_descriptions (bool): If True, print a description of all properties.
        
    Returns:
        Dictionary of spectral properties.
    """
    spectral_properties = {
        "disk_lum": np.array(output.spectral_properties.disk_lum),
        "IC_lum": np.array(output.spectral_properties.IC_lum),
        "xray_lum": np.array(output.spectral_properties.xray_lum),
        "radio_lum": np.array(output.spectral_properties.radio_lum),
        "xray_index": np.array(output.spectral_properties.xray_index),
        "radio_index": np.array(output.spectral_properties.radio_index),
        "jetbase_compactness": np.array(output.spectral_properties.jetbase_compactness),
    }
    
    if include_descriptions:
        description_text = (
            "The spectral properties include the following parameters:\n"
            "- 'disk_lum': Observed 0.3-5 keV disk luminosity.\n"
            "- 'IC_lum': Observed 0.3-300 keV Inverse Compton luminosity.\n"
            "- 'xray_lum': Observed 1-10 keV total luminosity.\n"
            "- 'radio_lum': Observed 4-6 GHz luminosity.\n"
            "- 'xray_index': X-ray 10-100 keV photon index estimate.\n"
            "- 'radio_index': Radio 10-100 GHz spectral index estimate.\n"
            "- 'jetbase_compactness': Jet base compactness.\n"
        )
        print(description_text)
    
    return spectral_properties



def plot_nufnu_ergshz(fits, fig_output_path=None, title="Emission Components", ax=None,
                      style=DEFAULT_STYLE):
    """Plot preprocessed BHJet components in ``nu F_nu`` space.

    This legacy notebook helper now delegates to the shared component plotter;
    edit ``plotting.DEFAULT_COMPONENT_STYLES`` or pass a ``PlotStyle`` to
    change every component consistently.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(13, 7))
    plot_bhjet_component_data(fits, ax=ax, style=style)
    ax.set_xlabel("Frequency (Hz)", fontsize=18)
    ax.set_ylabel(r"$\nu F_\nu$ (erg/cm2/s)", fontsize=18)
    if fig_output_path:
        ax.figure.savefig(fig_output_path, dpi=300)
    return ax.figure, ax



def plot_flux_mjy(data, output_path=None, title=None, ax=None, style=DEFAULT_STYLE):
    """Plot preprocessed BHJet components in flux-density space.

    The legacy function divided output fluxes by ``mjy_conv``; that scaling is
    retained so existing notebooks reproduce their current figures.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(13, 7))
    plot_bhjet_component_data(
        data, ax=ax, style=style.with_overrides(model_linewidth=1),
        flux_density=True, scale_factor=1 / mjy_conv, legend=True,
    )
    ax.set_xlabel("Frequency (Hz)", fontsize=14)
    ax.set_ylabel(r"$F_\nu$ (mJy)")
    ax.set_title(title, fontsize=16)
    ax.grid(True)
    if output_path:
        ax.figure.savefig(output_path, dpi=300)
    return ax.figure, ax



# def plot_radiative_zones_style_with_sizes(output, colors, size_cyclo_arr, size_com_arr,
#                                           kevconv=1.0, mjy=1e-26, fluxconv=1.0,
#                                           blim_f=1e8, ulim_f=1e26, blim_fl=1e-6, ulim_fl=1e2):
#     """
#     Plot radiative zone spectra using provided size arrays and color gradient.
#     """
#     nzones = 100

#     cyclosyn = output.cyclosyn_zones
#     compton = output.compton_zones

#     fig, ax1 = plt.subplots(1, 1, figsize=(7.5, 6))

#     totindex1 = 0
#     totindex2 = 0

#     for i in range(nzones):
#         n_cyclo = int(size_cyclo_arr[i])
#         n_com = int(size_com_arr[i])

#         if n_cyclo > 0:
#             nu_cyclosyn = np.array([cyclosyn[totindex1 + j].energy / kevconv for j in range(n_cyclo)])
#             lnu_cyclosyn = np.array([cyclosyn[totindex1 + j].flux * cyclosyn[totindex1 + j].energy * mjy * kevconv
#                                      for j in range(n_cyclo)])
#             totindex1 += n_cyclo
#         else:
#             nu_cyclosyn = lnu_cyclosyn = []

#         if n_com > 0:
#             nu_compton = np.array([compton[totindex2 + j].energy / kevconv for j in range(n_com)])
#             lnu_compton = np.array([compton[totindex2 + j].flux * compton[totindex2 + j].energy * mjy * kevconv
#                                     for j in range(n_com)])
#             totindex2 += n_com
#         else:
#             nu_compton = lnu_compton = []

#         # Plot every 3rd zone, or whatever 
#         if i % 3 == 0 and (n_cyclo > 0 or n_com > 0):
#             zorder = nzones - i
#             if len(nu_cyclosyn):
#                 ax1.plot(nu_cyclosyn, lnu_cyclosyn * fluxconv,
#                          linewidth=2.0, color=colors[i], zorder=zorder, linestyle='--')
#             if len(nu_compton):
#                 ax1.plot(nu_compton, lnu_compton * fluxconv,
#                          linewidth=2.0, color=colors[i], zorder=zorder, linestyle='--')

#     # Plot total spectrum if available
#     total = output.total
#     total_nu = np.array([pt.energy / kevconv for pt in total])
#     total_flux = np.array([pt.flux * pt.energy * mjy * fluxconv for pt in total])
#     ax1.plot(total_nu, total_flux, linewidth=2.5, color='black', zorder=nzones + 1)

#     # Axis and scale settings
#     ax1.set_ylim([0.001 * blim_fl * fluxconv, 0.1 * ulim_fl * fluxconv])
#     ax1.set_xlim([blim_f / kevconv, ulim_f / kevconv])
#     ax1.set_xscale('log', base=10)
#     ax1.set_yscale('log', base=10)

#     ax1.set_xlabel('Energy (keV)' if kevconv != 1 else 'Frequency (Hz)', fontsize=18)
#     ax1.set_ylabel('Luminosity (erg/s)' if fluxconv != 1 else 'Flux (erg/s/cm²)', fontsize=18)

#     plt.tight_layout()
