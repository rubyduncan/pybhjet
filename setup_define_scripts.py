import astromodels
import yaml
import numpy as np
import matplotlib.pyplot as plt
import builtins
from threeML import XYLike, OGIPLike
from threeML.utils.OGIP.response import OGIPResponse
from threeML import *
import sys
sys.path.append("PyBHJet/")
from pybhjet_3ml import BHJetModel
from pathlib import Path
from importing_data import * #converts normal data into expected threeml units 
from unit_conversion import *
# from astromodels.xspec import *
import re

#dictionary of all models that might be used for different sources, yaml will grab ones as needed, 
#only uncomment the xspec ones if it is running on helios
TOTAL_MODEL_LIB = {
    "BHJetModel": BHJetModel,
    "TbAbs": TbAbs,
    "ZDust": ZDust,
    # "zpcfabs": XS_zpcfabs, 
    # "pexmon": XS_pexmon, 
    # "apec": XS_apec,
}

def load_data_from_yaml(path_to_data_yaml_file):
    
    with open(path_to_data_yaml_file) as f: 
        yaml_dict = yaml.safe_load(f) 

    #base directory, or the folder that holds the yaml file
    base_dir = Path(path_to_data_yaml_file).parent 
    data_cfg = yaml_dict["data"] #specifically the data section 

    data_dict = {}
 
    rad_data = ir_data = uv_data = None
    plugin = None

    for name, comp in data_cfg.items(): #for data section in yaml file 
        kind = comp['kind'].lower()

        if kind == "combine_dataframes":
            dir = base_dir/comp['directory'] #combine dataframes expects the path to the directory w/ all data files
            extension = comp.get("extension", ".dat") #frequency spec. data file being loaded 
            columns = comp.get("columns", None) #this is to know the components of each file 

            rad, ir, uv = combine_dataframes(str(dir), extension, columns=columns)

            data_dict["rad"] = XYLike.from_dataframe("rad", rad, "x", "y", "yerr", False)
            data_dict["ir"] = XYLike.from_dataframe("ir", ir, "x", "y", "yerr", False)
            data_dict["uv"] = XYLike.from_dataframe("uv", uv, "x", "y", "yerr", False)

        elif kind == 'ogip':
            obs = str(base_dir/comp['observation'])
            arf = str(base_dir/comp['arf_file'])
            rmf = str(base_dir/comp['response'])
            bkg = str(base_dir/comp['background'])

            plugin = OGIPLike(name, observation=obs, arf_file=arf, response=rmf, background=bkg)

            if "energy_range" in comp: 
                plugin.set_active_measurements(comp['energy_range'])
            
            if "rebin_on_source" in comp:
                plugin.rebin_on_source(int(comp['rebin_on_source']))

            data_dict[name] = plugin 
            
    return data_dict

def load_params_priors_from_yaml(func, cfg_block): 

    for name, spec in cfg_block.items():
         #fixed parmeters don't have everything defined for them
        if not hasattr(func, name):
            continue

        p = getattr(func, name)

        # bounds - first so that parameter changes don't just throw errors
        if "bounds" in spec:
            lo, hi = spec["bounds"]
            p.bounds = (lo, hi)

        # value 
        if "value" in spec:
            p.value = spec["value"]

        # free / fixed
        if "free" in spec:
            p.free = bool(spec["free"])

        # priors 
        prior_cfg = spec.get("prior")
        if prior_cfg:
            ptype = prior_cfg["type"].lower()

            if ptype == "uniform":
                p.prior = Uniform_prior(
                    lower_bound=prior_cfg["min"],
                    upper_bound=prior_cfg["max"],
                )
            elif ptype == "log_uniform":
                p.prior = Log_uniform_prior(
                    lower_bound=prior_cfg["min"],
                    upper_bound=prior_cfg["max"],
                )
            elif ptype == "gaussian":
                p.prior = Gaussian(
                    mu=prior_cfg["mu"],
                    sigma=prior_cfg["sigma"],
                )
            elif ptype == 'truncated_gaussian':
                p.prior = Truncated_gaussian(
                mu=prior_cfg["mu"],
                sigma=prior_cfg["sigma"],
                lower_bound=prior_cfg.get("min", p.min_value),
                upper_bound=prior_cfg.get("max", p.max_value),
                )

def build_components_from_yaml(MODEL_SETUP, path_to_model_yaml):

    with open(path_to_model_yaml) as f: 
        yaml_dict = yaml.safe_load(f)

    comp_cfg  = yaml_dict["model"]["components"] 
    param_cfg = yaml_dict.get("parameters", {})

    model_components = {}

    for logical_name, class_name in comp_cfg.items(): #this is how it makes it into threeml exp. format
        if class_name not in MODEL_SETUP:
            raise KeyError(f"forgot to add this model in: {class_name!r} !!")
        cls = MODEL_SETUP[class_name]
        obj = cls()

        if logical_name in param_cfg:
            load_params_priors_from_yaml(obj, param_cfg[logical_name])

        model_components[logical_name] = obj

    return model_components, yaml_dict


_TOKEN = re.compile(r"\s*([A-Za-z_]\w*|\+|\*|\(|\))\s*")

def tokenize(expr):
    tokens = []
    i = 0
    while i < len(expr):
        m = _TOKEN.match(expr, i)
        tokens.append(m.group(1))
        i = m.end()
    return tokens

def to_rpn(tokens):
    precedence = {"+": 1, "*": 2}
    output = []
    operators = []

    for t in tokens:
        if t[0].isalpha() or t[0] == "_":
            output.append(t)

        elif t in ("+", "*"):
            while operators and operators[-1] in precedence:
                if precedence[operators[-1]] >= precedence[t]:
                    output.append(operators.pop())
                else:
                    break
            operators.append(t)

        elif t == "(":
            operators.append(t)

        else:  # for t == ")"
            while operators[-1] != "(":
                output.append(operators.pop())
            operators.pop()

    while operators:
        output.append(operators.pop())

    return output

def eval_rpn(rpn, components):
    stack = []

    for t in rpn:
        if t == "+":
            right = stack.pop()
            left = stack.pop()
            stack.append(left + right)

        elif t == "*":
            right = stack.pop()
            left = stack.pop()
            stack.append(left * right)

        else:
            stack.append(components[t])

    return stack[0]


def build_spectrum(expr, components):
    '''
    this is to build the spectrum for any number of 
    models and addition/multiplicative components
    just expects that the entire line from spectra (in yaml) is fed in, and uses + or * and ()
    these functions are readable by threeml (tested)
    
'''
    tokens = tokenize(expr)
    rpn = to_rpn(tokens)
    return eval_rpn(rpn, components)

# def build_spectrum(expr, components):
#     names = [s.strip() for s in expr.split("*")]
#     f = components[names[0]]
#     for name in names[1:]:
#         f = components[name] * f
#     return f

def build_model_and_data_from_yaml(path_to_data_yaml_file, path_to_model_yaml_file, MODEL_SETUP=None):
    '''
    build_model_and_data_from_yaml --> constructs everything for the model to run & plot
    '''
    #this is how to use the entire dictonary of models and construct the components as needed from there
    # otherwise can still supply the specific models directly 
    if MODEL_SETUP is None:
        MODEL_SETUP = TOTAL_MODEL_LIB

    # 1. load data + yaml dict
    data_dict = load_data_from_yaml(path_to_data_yaml_file)

    # 2. components (jet, gal_ext, intr_ext, dust_ext) with params/prior applied
    model_components, model_yaml_dict = build_components_from_yaml(MODEL_SETUP, path_to_model_yaml_file)

    # 3. spectra (radio, iruv, xray) from yaml model.spectra
    spectra_cfg = model_yaml_dict["model"]["spectra"]
    spectra = {
        spec_name: build_spectrum(expr, model_components)
        for spec_name, expr in spectra_cfg.items()
    }
    # 3.5. this is what the plotting functions are going to use to make individual components 
    sed_components_cfg = model_yaml_dict["model"].get("sed_components", {})

    # 4. build point sources for each dataset from yaml model source 
    sources_cfg = model_yaml_dict["sources"]
    sources = {}

    for src_name, cfg in sources_cfg.items():
        ra = cfg.get("ra", 0.0)
        dec = cfg.get("dec", 0.0)
        spec_name = cfg["spectrum"]         # "radio", "iruv", "xray"
        spectral_shape = spectra[spec_name]

        src = PointSource(
            src_name,
            ra=ra,
            dec=dec,
            spectral_shape=spectral_shape,
        )
        sources[src_name] = src

    # 5. assign datasets to sources
    for src_name, cfg in sources_cfg.items():
        datasets = cfg.get("datasets", [])
        for dname in datasets:
            data_dict[dname].assign_to_source(src_name)

    # 6. threeml Model and DataList
    model_obj = Model(*sources.values())
    data_obj = DataList(*[data_dict[name] for name in data_dict.keys()])

    return model_obj, data_obj, model_components, data_dict, sources, sed_components_cfg


def kev_plot_xylike_with_ind_model(xy, ax=None, label=None, color_data="k", color_model="r"):

    '''this takes the loaded "xylike" data for threeml and plots it 
        with the model that was assigned to it, just over the range for that model
    '''
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    x = np.asarray(xy.x)
    y = np.asarray(xy.y)
    yerr = np.asarray(xy.yerr) if xy.has_errors else None

    # data
    if yerr is not None:
        ax.errorbar(x, y, yerr=yerr, fmt="o", ms=4, lw=1, label=label or xy.name,
                    color=color_data)
    else:
        ax.scatter(x, y, s=10, label=label or xy.name, color=color_data)

    # model evaluated on the same x grid
    m = np.asarray(xy.get_model())
    ax.plot(x, m, lw=1.5, color=color_model, label=f"{label or xy.name} model")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    return fig, ax


def hz_plot_xylike_data(xy, ax=None, label=None, color_data="k", color_model="r"):

    '''this takes the loaded "xylike" data for threeml and plots it, using my conversion.py notebook
    '''
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    x = np.asarray(xy.x)
    y = np.asarray(xy.y)
    yerr = np.asarray(xy.yerr) if xy.has_errors else None

    x_hz = kev_to_hz(x)
    y_mjy = photon_flux_density_to_mjy(y, x)
    y_err_mjy = photon_flux_density_to_mjy(yerr, x)

    # data
    if y_err_mjy is not None:
        ax.errorbar(x_hz, y_mjy*x_hz/1e26, yerr=y_err_mjy*x_hz/1e26, fmt="o", ms=4, lw=1, label=label or xy.name,
                    color=color_data)
        
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    return fig, ax


def plot_ogip_with_model(*ogips, model_obj, fig=None, model_labels=None):

    if len(ogips) == 1 and isinstance(ogips[0], (list, tuple)):
        ogip_list = builtins.list(ogips[0])
    else:
        ogip_list = builtins.list(ogips)
        
    if model_labels is None:
        model_labels = [None] * len(ogip_list)
   
    fig_out = None 

    for i, (ogip, label) in enumerate(zip(ogip_list, model_labels)):
        ogip.set_model(model_obj)
        if fig_out is None:
            if fig is None:
                fig_out = ogip.display_model(model_label=label)
            else:
                fig_out = ogip.display_model(model_subplot=fig.axes, model_label=label)
        else:
            fig_out = ogip.display_model(model_subplot=fig_out.axes, model_label=label)

    return fig_out


def old_hz_plot_model_space_sed(data_dict,model_components,e_min_keV=1e-9,e_max_keV=1e3,n_points=1000,radio_key="rad",ir_key="ir",uv_key="uv"):
    """
    Plot radio/IR/UV data plus BHJet components in hz, nuS_nu 

    Parameters
    ----------
    data_dict : dict
        Mapping dataset name -> plugin (from build_model_and_data_from_yaml).
        Must contain entries for radio_key, ir_key, uv_key.
    components : dict
        Mapping logical component name -> astromodels function
        (from build_components_from_yaml), must contain 'jet', 'gal_ext',
        'intr_ext', 'dust_ext'.
    e_min_keV, e_max_keV : float
        Energy range in keV for computing the model curves.
    n_points : int
        Number of points in log-space energy grid.
    radio_key, ir_key, uv_key : str
        Keys in data_dict for the radio / IR / UV XYLike plugins.
    """

    # pull the XYLike plugins from data_dict
    rad_data = data_dict[radio_key]
    ir_data  = data_dict[ir_key]
    uv_data  = data_dict[uv_key]

    # pull the component functions from the components dict
    jet      = model_components["jet"]
    gal_ext  = model_components["gal_ext"]
    intr_ext = model_components["intr_ext"]
    dust_ext = model_components["dust_ext"]

    fig, ax = plt.subplots(figsize=(10, 6))

    # data + model over data range, to show effects of each
    hz_plot_xylike_data(rad_data, ax=ax, color_data="red",    color_model="red")
    hz_plot_xylike_data(ir_data,  ax=ax, color_data="orange", color_model="orange")
    hz_plot_xylike_data(uv_data,  ax=ax, color_data="yellow", color_model="yellow")

    # BHJet components in model space
    ene = np.logspace(np.log10(e_min_keV), np.log10(e_max_keV), n_points)

    # just the jet
    jet_flux = jet(ene)  # ph / (cm^2 s keV)
    jet_flux_mjy = photon_flux_density_to_mjy(jet_flux, ene)

    # dust * jet
    dust_flux = dust_ext(ene) * jet(ene)
    dust_flux_mjy = photon_flux_density_to_mjy(dust_flux, ene)

    # gal * intr * jet
    abs_flux = gal_ext(ene) * intr_ext(ene) * jet(ene)
    abs_flux_mjy = photon_flux_density_to_mjy(abs_flux, ene)

    ene_hz = kev_to_hz(ene)

    ax.plot(ene_hz, jet_flux_mjy * ene_hz / 1e26,
            ls="-",  lw=1.5, color="darkblue",  label="jet")
    ax.plot(ene_hz, dust_flux_mjy * ene_hz / 1e26,
            ls="--", lw=1.5, color="lightblue", label="ZDust x jet")
    ax.plot(ene_hz, abs_flux_mjy * ene_hz / 1e26,
            ls=":",  lw=1.5, color="black",     label="Nh_Gal x Nh_src x jet")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(ncol=2, fontsize=8)

    return fig, ax

def hz_plot_model_space_sed(
    data_dict,
    model_components,
    sed_components_expr,
    data_keys=("rad", "ir", "uv"),
    e_min_keV=1e-9,
    e_max_keV=1e3,
    n_points=1000,
    data_colors=("red", "orange", "gold"),
    model_lw=1.5,
):
    """
     SED creator which makes individual YAML sed_components per dataset

    Parameters
    ----------
    data_dict : dict
        dataset name -> XYLike plugin (or any plugin hz_plot_xylike_data)
    model_components : dict
        component name -> astromodels function (built from YAML components)
    sed_components_expr : dict
        name -> expression string, e.g. {"jet":"jet", "xray_model":"gal_ext*(zpcfabs*jet+pexmon)"}
    data_keys : tuple[str]
        keys in data_dict to plot as data
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    for k, c in zip(data_keys, data_colors):
        if k in data_dict:
            hz_plot_xylike_data(data_dict[k], ax=ax, color_data=c, color_model=c)

    ene = np.logspace(np.log10(e_min_keV), np.log10(e_max_keV), n_points)
    ene_hz = kev_to_hz(ene)

    for name, expr in sed_components_expr.items():
        model = build_spectrum(expr, model_components)
        flux = model(ene)
        flux_mjy = photon_flux_density_to_mjy(flux, ene)
        ax.plot(ene_hz, flux_mjy * ene_hz / 1e26, lw=model_lw, label=name)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    return fig, ax


def quick_eval(model, data_dict, verbose=True):

    for name, plugin in data_dict.items():
        plugin.set_model(model)

    total_log_like = 0.0
    per_plugin_stats = {}

    for name, plugin in data_dict.items():
        ll = plugin.get_log_like()
        total_log_like += ll
        per_plugin_stats[name] = -2.0 * ll

    stat_total = -2.0 * total_log_like

    if verbose:
        print(f"Global -2 log L: {stat_total:.3f}")
        for name, stat in per_plugin_stats.items():
            print(f"  {name}: -2 log L ≈ {stat:.3f}")

    return stat_total, per_plugin_stats


def hz_eval_and_plot_sed(model,data_dict, components, e_min_keV=1e-9,e_max_keV=1e3,n_points=1000,verbose=True, radio_key='rad', ir_key='ir', uv_key='uv'):

    stat_total, per_plugin_stats = quick_eval(
        model,
        data_dict,
        verbose=verbose,
    )
    fig, ax = hz_plot_model_space_sed(data_dict,components, e_min_keV=e_min_keV,e_max_keV=e_max_keV,n_points=n_points)

    # annotate with fit statistic
    ax.set_title(f"BHJet SED, current params, -2 log L ≈ {stat_total:.1f}", fontsize=10)

    return fig, ax, stat_total, per_plugin_stats

def add_bhjet_radiative_components_to_plot(model_components,ax,e_min_keV=1e-9,e_max_keV=1e3,jet_comp=True,force_rerun=True,n_eval=2,):
    jet = model_components["jet"]

    jet_components_to_plot = ["presyn", "postsyn", "precom", "postcom"]
    all_components_to_plot = ["presyn", "postsyn", "precom", "postcom", "disk", "bb"]

    style_map = {
        "presyn":  dict(color="dodgerblue", ls="-",  lw=1.5, label="Syn, z < z_diss"),
        "postsyn": dict(color="darkblue",   ls="--", lw=1.5, label="Syn, z > z_diss"),
        "precom":  dict(color="lightgreen", ls="-",  lw=1.5, label="IC,  z < z_diss"),
        "postcom": dict(color="green",      ls=":",  lw=1.5, label="IC,  z > z_diss"),
        "disk":    dict(color="red",        ls="-.", lw=1.5, label="Disk"),
        "bb":      dict(color="orange",     ls="-.", lw=1.5, label="BB"),
    }

    plot_comp = jet_components_to_plot if jet_comp else all_components_to_plot

    old_enable = getattr(jet, "enable_detailed_output", False)

    old_infosw = jet.infosw.value

    jet.enable_detailed_output = True
    jet.infosw.value = 2

    if force_rerun:
        jet._cached_params = None

    E_eval = np.logspace(np.log10(e_min_keV), np.log10(e_max_keV), max(int(n_eval), 2))
    _ = jet(E_eval) #this is where it is re-run 

    comps = jet._last_components  # should exist it is populated 

    jet.enable_detailed_output = old_enable
    jet.infosw.value = old_infosw

    for name in plot_comp:
        if name not in comps:
            continue

        nu_hz = np.asarray(comps[name]["energy"], dtype=float)
        Snu_mjy = np.asarray(comps[name]["flux"], dtype=float)

        m = np.isfinite(nu_hz) & np.isfinite(Snu_mjy) & (nu_hz > 0)
        nu_hz = nu_hz[m]
        Snu_mjy = Snu_mjy[m]
        if nu_hz.size < 2:
            continue

        order = np.argsort(nu_hz)
        nu_hz = nu_hz[order]
        Snu_mjy = Snu_mjy[order]

        nuSnu = nu_hz * (Snu_mjy * 1e-26)  # erg / (s cm^2)
        style = style_map.get(name, dict(color="gray", ls="--", lw=1.0, label=name))
        ax.plot(nu_hz, nuSnu, **style)

    ax.legend(ncol=2, fontsize=8)
    return ax



def old_add_bhjet_radiative_components_to_plot(model_components, ax, e_min_keV=1e-9, e_max_keV=1e3, jet_comp=True): 

    jet = model_components['jet']
    jet.enable_detailed_output=True
    jet.infosw.value=2

    jet_components_to_plot = ['presyn', "postsyn", "precom", "postcom"]
    all_components_to_plot = ['presyn', "postsyn", "precom", "postcom", "disk", "bb"]

    style_map = {
        "presyn":  dict(color="dodgerblue",  ls="-",   lw=1.5, label="Syn, z < z_diss"),
        "postsyn": dict(color="darkblue",    ls="--",  lw=1.5, label="Syn, z > z_diss"),
        "precom":  dict(color="lightgreen",  ls="-",   lw=1.5, label="IC,  z < z_diss"),
        "postcom": dict(color="green",       ls=":",   lw=1.5, label="IC,  z > z_diss"),
        "disk":    dict(color="red",         ls="-.",  lw=1.5, label="Disk"),
        "bb":      dict(color="orange",      ls="-.",  lw=1.5, label="BB"),
        }
    
    #need to just re-evaluate the bhjet model again to generate the components @ this current par set 
    E_range = np.logspace(np.log10(e_min_keV), np.log10(e_max_keV), 1000)

    old_flag = getattr(jet, "enable_detailed_output", False) #whatever was earlier 
    _ = jet(E_range) #which will run and give jet.last_components with detailed output enabled 

    comps = getattr(jet, "_last_components", None) #I also don't really know what this line is returning, with getattr and None
    jet.enable_detailed_output = False
    

    if jet_comp==True: 
        plot_comp = jet_components_to_plot
    else: 
        plot_comp = all_components_to_plot

    for name in plot_comp:
        #so they have to be turned into arrays? and assigned a type? 
        energy_hz = np.asarray(comps[name]['energy'], float)
        flux_mjy = np.asarray(comps[name]['flux'], float)

        nu_f_nu = flux_mjy * energy_hz/1e26
        style = style_map.get(name, dict(color="gray", ls="--", lw=1.0, label=name))
        ax.plot(energy_hz, nu_f_nu, **style)

    
    ax.legend(ncol=2, fontsize=8)