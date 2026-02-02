import os
import sys
import glob
import pandas as pd
import numpy as np
import yaml 

from unit_conversion import *


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
            raw = pd.read_table(file_path, names=columns, delim_whitespace=True, comment="#")
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



