import astromodels
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import builtins
from threeML import XYLike, OGIPLike
from threeML.utils.OGIP.response import OGIPResponse
from threeML import *
import sys
import hashlib
from collections import OrderedDict
sys.path.append("PyBHJet/")
from pybhjet_3ml import BHJetModel
from pathlib import Path
from importing_data import * #converts normal data into expected threeml units 
from astromodels.xspec import *
import re

#dictionary of all models that might be used for different sources, yaml will grab ones as needed, 
#only add the xspec ones if you're going to load the astromodels xspec 
TOTAL_MODEL_LIB = {
    "BHJetModel": BHJetModel,
    "TbAbs": TbAbs,
    "ZDust": ZDust,
    "zpcfabs": XS_zpcfabs,
    "zphabs": XS_zphabs,
    "zxipcf": XS_zxipcf,  
    "pexmon": XS_pexmon, 
    "apec": XS_apec,
    "Powerlaw": XS_powerlaw,
    "Constant": Constant, 
}


def load_full_detailed_yaml_as_xylike(base_dir, comp):
    '''
    This will load the detailed data yaml file with citations, convert the units to photon flux, 
    and return the 3 correct dataframes of radio, ir, and uv to use in the load_data_from_yaml function. 

    it will look for sections like "radio" , "submm" and the threeml script just expects that these will be returned
    in three groups (radio, ir, uv) and you can store whichever data inside those three. (e.g. see grouped data below)
    '''

    #opening the second yaml from the first, this just has the 
    yaml_path = base_dir / comp["directory"]

    with open(yaml_path) as f:
        dat = yaml.safe_load(f)
    data_dict = {}
    grouped_data = {
        "rad": dat.get("radio", []) + dat.get("submm", []),
        "ir": dat.get("infrared", []),
        "uv": dat.get("optical", []),
    }

    for name, entries in grouped_data.items():
        df = pd.DataFrame(entries).copy()
        nu_hz = df["frequency_Hz"].astype(float).to_numpy()
        flux_mjy = df["flux_mjy"].astype(float).to_numpy()
        flux_err_mjy = df["flux_err"].astype(float).to_numpy()

        e_keV = hz_to_kev(nu_hz)
        y_ph = mjy_to_diff_photon_flux(flux_mjy, nu_hz)
        yerr_ph = mjy_to_diff_photon_flux(flux_err_mjy, nu_hz)
        new_df = pd.DataFrame({"x": e_keV,"y": y_ph,"yerr": yerr_ph})

        data_dict[name] = XYLike.from_dataframe(name, new_df, "x", "y", "yerr", False)

    return data_dict


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

        elif kind == 'full_detailed_yaml':
            data_dict.update(load_full_detailed_yaml_as_xylike(base_dir, comp))

        elif kind == 'ogip':
            obs = str(base_dir/comp['observation'])
            arf = str(base_dir/comp['arf_file'])
            rmf = str(base_dir/comp['response'])
            
            bkg_val = comp.get("background", None) #retrieves the value for the background key, and otherwise None 
            has_bkg = bkg_val not in (None, "", "none", "null")

            if has_bkg:
                bkg = str(base_dir / bkg_val)
                plugin = OGIPLike(name, observation=obs, arf_file=arf, response=rmf, background=bkg)
            else:
                plugin = OGIPLike(name, observation=obs, arf_file=arf, response=rmf)

            # plugin = OGIPLike(name, observation=obs, arf_file=arf, response=rmf, background=bkg)
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
                lower = prior_cfg.get("min", prior_cfg.get("lower_bound", p.min_value))
                upper = prior_cfg.get("max", prior_cfg.get("upper_bound", p.max_value))

                p.prior = Truncated_gaussian(
                mu=prior_cfg["mu"],
                sigma=prior_cfg["sigma"],
                lower_bound=lower,
                upper_bound=upper,
                )
            elif ptype == 'log_normal':
                mu_lin = prior_cfg['mu']
                if mu_lin <= 0:
                    raise ValueError(f"{name}: log_normal prior requires mu > 0 in linear space, got {mu_lin}")
            
                p.prior = Log_normal(
                F=1.0,
                mu=np.log(prior_cfg["mu"]),
                sigma=prior_cfg["sigma"],
                piv=prior_cfg.get("piv", 1.0),
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


def link_params(model_obj, model_components, driver: str, dependent: str, freeze_dependent: bool = True):
    comp_d, par_d = driver.split(".", 1)
    comp_y, par_y = dependent.split(".", 1)

    d = getattr(model_components[comp_d], par_d)
    y = getattr(model_components[comp_y], par_y)

    model_obj.link(y, d)

    if freeze_dependent:
        y.free = False

def apply_links(model_obj, model_components, model_yaml_dict):
    links = model_yaml_dict.get("links") #will skip for none 
    if not links:
        return

    for entry in links:
        link_params(
            model_obj,
            model_components,
            entry["driver"],
            entry["dependent"],
            bool(entry.get("freeze_dependent", True)),
        )

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


def apply_jet_config(model_components, model_yaml_dict):
    '''This is to set the cutoff type before the model runs '''
    cfg = model_yaml_dict.get("model_switches", {})
    jet_cfg = cfg.get("jet", {})

    jet = model_components["jet"]  # This is meant to map to the BHJet Model that we imported earlier 
    jet.cutoff_type = int(jet_cfg["cutoff_type"])


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
    apply_jet_config(model_components, model_yaml_dict) ## this controls the switches in bhjet 
    apply_component_caching(model_components, model_yaml_dict) #this helps with the xspec models and how they are evaluated



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

    apply_links(model_obj, model_components, model_yaml_dict)

    return model_obj, data_obj, model_components, data_dict, sources, sed_components_cfg


def enable_xspec_cache_evaluate(func):
    '''astromodels.functions.function.Function.fast_call and then goes to the xspec binding
    '''
    if getattr(func, "_cache_enabled", False): #cancel the caching 
        return

    orig_evaluate = func.evaluate
    func._orig_evaluate = orig_evaluate

    func._cache_enabled = True
    func._cache_hits = 0
    func._cache_misses = 0
    func._cache_last_params_key = None
    func._cache_last_x = None
    func._cache_last_y = None

    def cached_evaluate(x, *args, **kwargs): #targets the evaluate function 
        ps = func.parameters
        keys = sorted(ps.keys())
        vals = np.array([float(ps[k].value) for k in keys], dtype=np.float64)
        params_key = (tuple(keys), vals.tobytes())

        x = np.asarray(x, dtype=np.float64)

        if func._cache_last_params_key == params_key and func._cache_last_y is not None:
            func._cache_hits += 1
            return np.interp(x, func._cache_last_x, func._cache_last_y)

        y = orig_evaluate(x, *args, **kwargs)

        func._cache_last_params_key = params_key
        func._cache_last_x = np.array(x, copy=True)
        func._cache_last_y = np.array(y, copy=True)
        func._cache_misses += 1
        return func._cache_last_y

    func.evaluate = cached_evaluate


def enable_xspec_cache_exact(func, max_grids_per_param=4):
    if getattr(func, "_cache_enabled", False):
        return

    if not hasattr(func, "fast_call") or func.fast_call is None:
        raise RuntimeError("Cannot cache: func.fast_call is missing or None")

    func._orig_fast_call = func.fast_call
    func._cache_enabled = True
    func._cache_hits = 0
    func._cache_misses = 0

    # params_key -> OrderedDict[xhash -> y]
    func._cache_store = {}

    def params_key():
        ps = func.parameters
        keys = sorted(ps.keys())
        vals = np.array([float(ps[k].value) for k in keys], dtype=np.float64)
        return (tuple(keys), vals.tobytes())

    def x_key(x):
        x = np.asarray(x, dtype=np.float64)
        h = hashlib.blake2b(digest_size=16)
        h.update(x.tobytes())
        return h.digest()

    def cached_fast_call(x):
        pk = params_key()
        xk = x_key(x)

        store = func._cache_store.get(pk)
        if store is None:
            store = OrderedDict()
            func._cache_store[pk] = store

        if xk in store:
            func._cache_hits += 1
            y = store.pop(xk)
            store[xk] = y
            return y

        y = func._orig_fast_call(x)
        y = np.array(y, copy=True)

        store[xk] = y
        while len(store) > max_grids_per_param:
            store.popitem(last=False)

        func._cache_misses += 1
        return y

    func.fast_call = cached_fast_call

def disable_xspec_cache_exact(func):
    if hasattr(func, "_orig_fast_call") and func._orig_fast_call is not None:
        func.fast_call = func._orig_fast_call
    func._cache_enabled = False

def apply_component_caching(model_components, model_yaml_dict):
    cfg = model_yaml_dict.get("config", {}) or {}
    names = cfg.get("cache_components", []) or []
    if not names:
        return

    for name in names:
        enable_xspec_cache_exact(model_components[name])

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


def lum_plot_xylike_data(model_components, xy, ax=None, label=None, color_data="k"):

    '''this takes the loaded "xylike" data for threeml and plots it, using my conversion.py notebook
    '''

    jet = model_components["jet"]
    dist_kpc = jet.dist.value 

    fluxconv = 4.0 * np.pi * (dist_kpc * 3.085677581e21) ** 2
    mjy_to_cgs = 1e-26
    
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
    Lnu = y_mjy * mjy_to_cgs * fluxconv
    Lnu_err = y_err_mjy * mjy_to_cgs * fluxconv

    if y_err_mjy is not None:
        ax.errorbar(x_hz, x_hz*Lnu, yerr=Lnu_err*x_hz, fmt="o", ms=4, lw=1,
                    color=color_data)
        
    ax.set_xscale("log")
    ax.set_yscale("log")
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


def hz_plot_model_space_sed(
    data_dict,
    model_components,
    sed_components_expr,
    xray_path=None,
    data_keys=("rad", "ir", "uv"),
    e_min_keV=1e-9,
    e_max_keV=1e4,
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


    if xray_path is not None:
        add_xray_to_ax_nu_f_nu(
            ax,
            xray_path,
            model_components,
            lum=False,
            color = 'black'
        )

    ene = np.logspace(np.log10(e_min_keV), np.log10(e_max_keV), n_points)
    ene_hz = kev_to_hz(ene)

    for name, expr in sed_components_expr.items():
        model = build_spectrum(expr, model_components)

        if "zpcfabs" in expr or "gal_ext" in expr or "pexmon" in expr:
            ene_use = np.logspace(np.log10(0.1), np.log10(300.0), n_points)
        else:
            ene_use = ene

        ene_hz_use = kev_to_hz(ene_use)
        ph_flux = model(ene_use)
        fnu_mjy = photon_flux_density_to_mjy(ph_flux, ene_use)

        ax.plot(ene_hz_use, fnu_mjy*ene_hz_use/1e26, lw=model_lw, label=name)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    return fig, ax


def add_xray_to_ax_nu_f_nu(ax,path,model_components,*,color,
    y_col=2,
    yerr_col=3,
    fmt="D",
    markersize=8,
    alpha=0.5,
    zorder=1,
    label=None,
    lum = False,
    **errorbar_kwargs,
):
    jet = model_components["jet"]
    dist_kpc = jet.dist.value 

    fluxconv = 4.0 * np.pi * (dist_kpc * 3.085677581e21) ** 2
    mjy_to_cgs = 1e-26

    path = Path(path)

    def read_data(path):
        if os.path.exists(path):
            return np.genfromtxt(path)
        else:
            print(f"File not found: {path}")
            return None

    data = read_data(str(path))
    if data is None:
        return None

    data = np.asarray(data)
    if data.ndim != 2 or data.shape[1] <= max(1, y_col, yerr_col):
        raise ValueError(
            f"Unexpected shape for {path}: got {data.shape}, need 2D with enough columns."
        )

    freq = 10.0 ** ((np.log10(data[:, 0]) + np.log10(data[:, 1])) / 2.0)

    if lum == True: 
        
        y = data[:, y_col]
        yerr= data[:, yerr_col]

        Lnu = y * fluxconv 
        Lnu_err = yerr * fluxconv

        nuLnu =  Lnu
        nuLnu_err = Lnu_err 

        handle = ax.errorbar(
        freq,
        nuLnu,
        yerr=nuLnu_err,
        fmt=fmt,
        color=color,
        markersize=markersize,
        alpha=alpha,
        zorder=zorder,
        label=label,
        **errorbar_kwargs,
        )
        
    else:
        y = data[:, y_col] 
        y_err = data[:, yerr_col]

        nuFnu = y
        nuFnu_err = y_err
        
        handle = ax.errorbar(
        freq,
        nuFnu,
        yerr=nuFnu_err,
        fmt=fmt,
        color=color,
        markersize=markersize,
        alpha=alpha,
        zorder=zorder,
        label=label,
        **errorbar_kwargs,
        ) 
    return handle

def hz_plot_all_data_w_fake_xray(ax,
    xray_path,
    data_dict,
    model_components,
    sed_components_expr,
    data_keys=("rad", "ir", "uv"),
    e_min_keV=1e-9,
    e_max_keV=1e6,
    n_points=1000,
    data_colors=("red", "orange", "gold"),
    model_lw=1.5,
):


    for k, c in zip(data_keys, data_colors):
        if k in data_dict:
            hz_plot_xylike_data(data_dict[k], ax=ax, color_data=c, color_model=c)

    ene = np.logspace(np.log10(e_min_keV), np.log10(e_max_keV), n_points)
    ene_hz = kev_to_hz(ene)

    add_xray_to_ax_nu_f_nu(ax, xray_path, model_components, lum = False, color = 'black')
    


def lum_plot_all_data_w_fake_xray(
    xray_path,
    data_dict,
    model_components,
    data_keys=("rad", "ir", "uv"),
    e_min_keV=1e-9,
    e_max_keV=1e6,
    n_points=1000,
    data_colors=("red", "orange", "gold"),
    zorder = 10,
    
):

    jet = model_components["jet"]
    dist_kpc = jet.dist.value 

    fluxconv = 4.0 * np.pi * (dist_kpc * 3.085677581e21) ** 2
    mjy_to_cgs = 1e-26

    fig, ax = plt.subplots(figsize=(10, 6))

    for k, c in zip(data_keys, data_colors):
        if k in data_dict:
            lum_plot_xylike_data(model_components, data_dict[k], ax=ax, color_data=c, zorder=zorder)

    ene = np.logspace(np.log10(e_min_keV), np.log10(e_max_keV), n_points)
    ene_hz = kev_to_hz(ene)

    add_xray_to_ax_nu_f_nu(ax, xray_path, model_components, lum=True, color = 'black', zorder=zorder, alpha = 0.7)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    return fig, ax

def quick_eval(model, data_dict, verbose=True, reset_model=False):
    #this resets the plugin model state, so will invalidate the cache if True 

    if reset_model==True:
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


def add_bhjet_radiative_components_to_plot(model_components,ax,e_min_keV=1e-9,e_max_keV=1e3,plot_mode='jet',force_rerun=True,n_eval=2,):
    jet = model_components["jet"]

    total_component = ['total']
    jet_components_to_plot = ["presyn", "postsyn", "precom", "postcom"]
    all_components_to_plot = ["presyn", "postsyn", "precom", "postcom", "disk", "bb"]

    style_map = {
        "total" :  dict(color="black", ls="-",  lw=1.5, label="Total Jet Emission"),
        "presyn":  dict(color="dodgerblue", ls="-",  lw=1.5, label="Syn, z < z_diss"),
        "postsyn": dict(color="darkblue",   ls="--", lw=1.5, label="Syn, z > z_diss"),
        "precom":  dict(color="lightgreen", ls="-",  lw=1.5, label="IC,  z < z_diss"),
        "postcom": dict(color="green",      ls=":",  lw=1.5, label="IC,  z > z_diss"),
        "disk":    dict(color="red",        ls="-.", lw=1.5, label="Disk"),
        "bb":      dict(color="orange",     ls="-.", lw=1.5, label="BB"),
    }

    if plot_mode == "total":
        plot_comp = total_component
    elif plot_mode == "jet":
        plot_comp = jet_components_to_plot
    elif plot_mode == "all":
        plot_comp = all_components_to_plot

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


def bhjet_luminosity_components_to_plot(
    model_components,
    ax,
    total_color=None,
    total_label=None,
    e_min_keV=1e-9,
    e_max_keV=1e3,
    plot_mode="jet",
    force_rerun=True,
    n_eval=2,
    style_map=None,
    show_legend=True,
):
    jet = model_components["jet"]
    dist_kpc = jet.dist.value

    fluxconv = 4.0 * np.pi * (dist_kpc * 3.085677581e21) ** 2
    mjy_to_cgs = 1e-26

    total_component = ["total"]
    jet_components_to_plot = [
        "presyn",
        "postsyn",
        "precom",
        "postcom",
    ]
    all_components_to_plot = [
        "presyn",
        "postsyn",
        "precom",
        "postcom",
        "disk",
        "bb",
    ]

    default_style_map = {
        "total": dict(
            color="black",
            ls="-",
            lw=1.5,
            label="Total Jet Emission",
        ),
        "presyn": dict(
            color="dodgerblue",
            ls="-",
            lw=1.5,
            label=r"Syn, $z < z_{\rm diss}$",
        ),
        "postsyn": dict(
            color="darkblue",
            ls="--",
            lw=1.5,
            label=r"Syn, $z > z_{\rm diss}$",
        ),
        "precom": dict(
            color="lightgreen",
            ls="-",
            lw=1.5,
            label=r"IC, $z < z_{\rm diss}$",
        ),
        "postcom": dict(
            color="green",
            ls=":",
            lw=1.5,
            label=r"IC, $z > z_{\rm diss}$",
        ),
        "disk": dict(
            color="red",
            ls="-.",
            lw=1.5,
            label="Disk",
        ),
        "bb": dict(
            color="orange",
            ls="-.",
            lw=1.5,
            label="BB",
        ),
    }

    if style_map is None:
        style_map = {
            name: style.copy()
            for name, style in default_style_map.items()
        }
    else:
        merged_style_map = {
            name: style.copy()
            for name, style in default_style_map.items()
        }

        for name, style in style_map.items():
            if name in merged_style_map:
                merged_style_map[name].update(style)
            else:
                merged_style_map[name] = style.copy()

        style_map = merged_style_map

    if total_color is not None:
        style_map["total"]["color"] = total_color

    if total_label is not None:
        style_map["total"]["label"] = total_label

    if plot_mode == "total":
        plot_comp = total_component
    elif plot_mode == "jet":
        plot_comp = jet_components_to_plot
    elif plot_mode == "all":
        plot_comp = all_components_to_plot
    else:
        raise ValueError(
            "plot_mode must be 'total', 'jet', or 'all'."
        )

    old_enable = getattr(jet, "enable_detailed_output", False)
    old_infosw = jet.infosw.value

    jet.enable_detailed_output = True
    jet.infosw.value = 2

    try:
        if force_rerun:
            jet._cached_params = None

        E_eval = np.logspace(
            np.log10(e_min_keV),
            np.log10(e_max_keV),
            max(int(n_eval), 2),
        )

        _ = jet(E_eval)
        comps = jet._last_components

    finally:
        jet.enable_detailed_output = old_enable
        jet.infosw.value = old_infosw

    for name in plot_comp:
        if name not in comps:
            continue

        nu_hz = np.asarray(
            comps[name]["energy"],
            dtype=float,
        )
        Snu_mjy = np.asarray(
            comps[name]["flux"],
            dtype=float,
        )

        m = (
            np.isfinite(nu_hz)
            & np.isfinite(Snu_mjy)
            & (nu_hz > 0)
            & (Snu_mjy > 0)
        )

        nu_hz = nu_hz[m]
        Snu_mjy = Snu_mjy[m]

        if nu_hz.size < 2:
            continue

        order = np.argsort(nu_hz)
        nu_hz = nu_hz[order]
        Snu_mjy = Snu_mjy[order]

        Lnu = Snu_mjy * mjy_to_cgs * fluxconv
        nuLnu = nu_hz * Lnu

        style = style_map.get(
            name,
            dict(
                color="gray",
                ls="--",
                lw=1.0,
                label=name,
            ),
        )

        ax.plot(
            nu_hz,
            nuLnu,
            zorder=200,
            **style,
        )

    if show_legend:
        ax.legend(
            ncol=2,
            fontsize=8,
            frameon=False,
        )

    return ax


def old_bhjet_luminosity_components_to_plot(model_components,ax,total_color=None,total_label=None,e_min_keV=1e-9,e_max_keV=1e3,plot_mode='jet',force_rerun=True,n_eval=2, style_map=None,):
    jet = model_components["jet"]
    dist_kpc = jet.dist.value 

    fluxconv = 4.0 * np.pi * (dist_kpc * 3.085677581e21) ** 2
    mjy_to_cgs = 1e-26

    total_component = ['total']
    jet_components_to_plot = ["presyn", "postsyn", "precom", "postcom"]
    all_components_to_plot = ["presyn", "postsyn", "precom", "postcom", "disk", "bb"]

    style_map = {
        # "total" :  dict(color="black", ls="-",  lw=1.5, label="Total Jet Emission"),
        "total" :  dict(color="black", ls="-",  lw=1.5, label="Total Jet Emission"),
        "presyn":  dict(color="dodgerblue", ls="-",  lw=1.5, label="Syn, z < z_diss"),
        "postsyn": dict(color="darkblue",   ls="--", lw=1.5, label="Syn, z > z_diss"),
        "precom":  dict(color="lightgreen", ls="-",  lw=1.5, label="IC,  z < z_diss"),
        "postcom": dict(color="green",      ls=":",  lw=1.5, label="IC,  z > z_diss"),
        "disk":    dict(color="red",        ls="-.", lw=1.5, label="Disk"),
        "bb":      dict(color="orange",     ls="-.", lw=1.5, label="BB"),
    }

    if total_color is not None:
        style_map["total"]["color"] = total_color

    if total_label is not None:
        style_map["total"]["label"] = total_label

    if plot_mode == "total":
        plot_comp = total_component
    elif plot_mode == "jet":
        plot_comp = jet_components_to_plot
    elif plot_mode == "all":
        plot_comp = all_components_to_plot

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

        # nuSnu = nu_hz * (Snu_mjy * 1e-26)  # erg / (s cm^2)
        Lnu = Snu_mjy * mjy_to_cgs * fluxconv
        style = style_map.get(name, dict(color="gray", ls="--", lw=1.0))
        ax.plot(nu_hz, nu_hz*Lnu, **style)

    ax.legend(ncol=2, fontsize=8)
    return ax


#unit conversions ----------------- 

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



def combine_dataframes(directory_path, file_extension, columns=None):
    """
    Will load files that have the columns:
        [nu (Hz), flux (mJy), err (mJy)]
    Then converts each file to a new dataframe with 
        x    = E [keV]
        y    = dN/dE [ph cm^-2 s^-1 keV^-1]
        yerr = same 
    put files into radio / IR / UV based on filename naming
    return (radio_df, ir_df, uv_df) for what 3ML expects
    """

    radio_df_list = []
    ir_df_list    = []
    uv_df_list    = []

    pattern = os.path.join(directory_path, f"*{file_extension}")
    for file_path in glob.glob(pattern):

        fname = os.path.basename(file_path)

        # decide category from filename
        if "radio" in fname.lower():
            category = "radio"
        elif "UV" in fname.upper():
            category = "UV"
        elif "IR" in fname.upper():
            category = "IR"
        else:
            continue

        # read  file
        if columns:
            raw = pd.read_table(file_path, names=columns, sep='\s+', comment="#")
        else:
            raw = pd.read_table(file_path, sep=r"\s+", header=None, comment="#")

        if raw.shape[1] < 3:
            # not enough columns, skip
            continue

        # pull arrays
        nu_hz        = raw.iloc[:, 0].astype(float).to_numpy()  # frequency [Hz]
        flux_mJy     = raw.iloc[:, 1].astype(float).to_numpy()  # S_nu [mJy]
        flux_err_mJy = raw.iloc[:, 2].astype(float).to_numpy()  # err [mJy]

        # convert frequency to photon energy in keV
        e_keV = hz_to_kev(nu_hz)

        # convert mJy to photon differential flux - test 
        # y_ph     = s_nu_mjy_to_diff_photon_flux(flux_mJy,     nu_hz)
        # yerr_ph  = s_nu_mjy_to_diff_photon_flux(flux_err_mJy, nu_hz)

        y_ph     = mjy_to_diff_photon_flux(flux_mJy,     nu_hz)
        yerr_ph  = mjy_to_diff_photon_flux(flux_err_mJy, nu_hz)

        # clean dataframe for file
        clean_df = pd.DataFrame({
            "x":    e_keV,
            "y":    y_ph,
            "yerr": yerr_ph,
        })

        # append to correct data list
        if category == "radio":
            radio_df_list.append(clean_df)
        elif category == "IR":
            ir_df_list.append(clean_df)
        elif category == "UV":
            uv_df_list.append(clean_df)
        else:
            continue

    # concatenate within each category
    if radio_df_list:
        radio_df = pd.concat(radio_df_list, ignore_index=True)
    else:
        radio_df = pd.DataFrame(columns=["x", "y", "yerr"])

    if ir_df_list:
        ir_df = pd.concat(ir_df_list, ignore_index=True)
    else:
        ir_df = pd.DataFrame(columns=["x", "y", "yerr"])

    if uv_df_list:
        uv_df = pd.concat(uv_df_list, ignore_index=True)
    else:
        uv_df = pd.DataFrame(columns=["x", "y", "yerr"])

    return radio_df, ir_df, uv_df

