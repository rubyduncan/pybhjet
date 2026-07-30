"""Composable plotting helpers for BHJet and threeML analyses.

This module deliberately contains no model construction or YAML loading.  Pass it
already-built threeML plugins and model components from ``setup_define_scripts``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np


KEV_TO_HZ = 2.42e17
MJY_TO_CGS = 1e-26
KPC_TO_CM = 3.085677581e21


# This is the one place to change BHJet component colours, dashes, or labels.
# ``style`` is retained as an accepted alias for old notebook dictionaries.
DEFAULT_COMPONENT_STYLES = {
    "total": {"color": "black", "linestyle": "-", "label": "Total"},
    "presyn": {"color": "dodgerblue", "linestyle": "-", "label": r"Syn, $z < z_{\rm diss}$"},
    "postsyn": {"color": "darkblue", "linestyle": (0, (5, 1)), "label": r"Syn, $z > z_{\rm diss}$"},
    "precom": {"color": "lightgreen", "linestyle": "-", "label": r"IC, $z < z_{\rm diss}$"},
    "postcom": {"color": "green", "linestyle": (0, (3, 1, 1, 1)), "label": r"IC, $z > z_{\rm diss}$"},
    "disk": {"color": "red", "linestyle": (0, (3, 2, 1, 2, 1, 2)), "label": "Disk"},
    "bb": {"color": "orange", "linestyle": "-", "label": "Blackbody"},
}


@dataclass(frozen=True)
class PlotStyle:
    """Visual defaults, overridable per figure without changing global rcParams."""

    figsize: tuple[float, float] = (10, 6)
    data_marker: str = "o"
    data_markersize: float = 5
    data_alpha: float = 0.85
    model_linewidth: float = 1.5
    legend_fontsize: float = 9
    confidence_alpha_95: float = 0.16
    confidence_alpha_68: float = 0.32
    data_colors: Mapping[str, str] = field(
        default_factory=lambda: {"rad": "#0072B2", "ir": "#E69F00", "uv": "#009E73"}
    )
    component_styles: Mapping[str, Mapping[str, object]] = field(
        default_factory=lambda: {name: dict(spec) for name, spec in DEFAULT_COMPONENT_STYLES.items()}
    )

    def with_overrides(self, **overrides):
        """Return a copy with figure-wide fields changed, without global state."""
        return replace(self, **overrides)

    def with_component_overrides(self, **component_overrides):
        """Return a copy with only the named component styles changed.

        Example
        -------
        ``paper_style = DEFAULT_STYLE.with_component_overrides(
        presyn={"color": "navy", "label": "Compact-jet synchrotron"})``
        """
        merged = {name: dict(spec) for name, spec in self.component_styles.items()}
        for name, overrides in component_overrides.items():
            overrides = dict(overrides)
            if "style" in overrides and "linestyle" not in overrides:
                overrides["linestyle"] = overrides["style"]
            if "ls" in overrides and "linestyle" not in overrides:
                overrides["linestyle"] = overrides["ls"]
            if "lw" in overrides and "linewidth" not in overrides:
                overrides["linewidth"] = overrides["lw"]
            merged.setdefault(name, {})
            merged[name].update(overrides)
        return replace(self, component_styles=merged)


DEFAULT_STYLE = PlotStyle()


def component_plot_kwargs(style, component, overrides=None):
    """Return Matplotlib-ready style keywords for one BHJet component.

    Existing notebooks sometimes use ``style``, ``ls`` and ``lw``. They are
    normalized here so the editable component mapping can use either familiar
    shorthand or Matplotlib's long names.
    """
    spec = dict(style.component_styles.get(component, {"label": component}))
    if "style" in spec and "linestyle" not in spec:
        spec["linestyle"] = spec["style"]
    if "ls" in spec and "linestyle" not in spec:
        spec["linestyle"] = spec["ls"]
    if "lw" in spec and "linewidth" not in spec:
        spec["linewidth"] = spec["lw"]
    incoming = dict(overrides or {})
    spec.update(incoming)
    if "style" in incoming and "linestyle" not in incoming:
        spec["linestyle"] = spec["style"]
    spec.pop("style", None)
    if "ls" in incoming and "linestyle" not in incoming:
        spec["linestyle"] = spec["ls"]
    spec.pop("ls", None)
    if "lw" in incoming and "linewidth" not in incoming:
        spec["linewidth"] = spec["lw"]
    spec.pop("lw", None)
    spec.setdefault("linewidth", style.model_linewidth)
    return spec


def is_ogip_plugin(plugin):
    """Return whether a threeML plugin is OGIP/count-space data.

    OGIP data are response-folded spectra and are intentionally never
    converted to the flux-space SED representation used in this module.
    """
    return any(cls.__name__ == "OGIPLike" for cls in type(plugin).__mro__)


def is_flux_space_plugin(plugin):
    """Return whether a plugin can be plotted as flux-space XY data."""
    return not is_ogip_plugin(plugin) and hasattr(plugin, "x") and hasattr(plugin, "y")


def split_plot_data(data):
    """Split a plugin mapping into flux-space and OGIP/count-space datasets."""
    flux_data, ogip_data = {}, {}
    for name, plugin in (data or {}).items():
        if is_ogip_plugin(plugin):
            ogip_data[name] = plugin
        elif is_flux_space_plugin(plugin):
            flux_data[name] = plugin
    return flux_data, ogip_data


def kev_to_hz(energy_kev):
    return np.asarray(energy_kev) * KEV_TO_HZ


def photon_flux_to_mjy(number_flux, energy_kev):
    """Convert differential photon flux to mJy at the supplied energy in keV."""
    return np.asarray(number_flux) * np.asarray(energy_kev) * 6.626e-27 / MJY_TO_CGS


def make_sed_axes(ax=None, *, style=DEFAULT_STYLE, luminosity=False, title=None):
    """Return log-log SED axes with consistent labels and no global state changes."""
    if ax is None:
        _, ax = plt.subplots(figsize=style.figsize)
    ax.set(xscale="log", yscale="log", xlabel=r"Frequency $\nu$ (Hz)")
    ax.set_ylabel(r"$\nu L_\nu$ (erg s$^{-1}$)" if luminosity else r"$\nu F_\nu$ (erg cm$^{-2}$ s$^{-1}$)")
    if title:
        ax.set_title(title)
    return ax


def _luminosity_factor(model_components):
    distance_kpc = model_components["jet"].dist.value
    return 4 * np.pi * (distance_kpc * KPC_TO_CM) ** 2


def plot_flux_points(points, *, ax=None, label=None, color=None, style=DEFAULT_STYLE,
                     luminosity=False, distance_kpc=None, scale_factor=1, annotate_dates=False,
                     **kwargs):
    """Plot a raw observation mapping with ``frequency_Hz``, ``flux_mjy``, and optional errors.

    This handles the compact SED dictionaries returned by the analysis notebooks,
    without requiring conversion to a threeML plugin first.
    """
    ax = make_sed_axes(ax, style=style, luminosity=luminosity)
    frequency = np.asarray(points["frequency_Hz"], dtype=float)
    flux = np.asarray(points["flux_mjy"], dtype=float)
    error = points.get("flux_err")
    valid = np.isfinite(frequency) & np.isfinite(flux) & (frequency > 0)
    frequency, flux = frequency[valid], flux[valid]
    factor = 1
    if luminosity:
        if distance_kpc is None:
            raise ValueError("distance_kpc is required for raw-data luminosity plots")
        factor = 4 * np.pi * (distance_kpc * KPC_TO_CM) ** 2
    y = frequency * flux * MJY_TO_CGS * factor * scale_factor
    yerr = None
    if error is not None:
        error = np.asarray(error, dtype=float)[valid]
        yerr = frequency * error * MJY_TO_CGS * factor * scale_factor
    artist = ax.errorbar(frequency, y, yerr=yerr, fmt=style.data_marker,
                         ms=style.data_markersize, alpha=style.data_alpha,
                         color=color, label=label, **kwargs)
    if annotate_dates and "date" in points:
        for x, yy, date in zip(frequency, y, np.asarray(points["date"])[valid]):
            if date is not None:
                ax.annotate(str(date), (x, yy), fontsize=8)
    return artist


def plot_xylike(plugin, *, ax=None, label=None, color=None, style=DEFAULT_STYLE, luminosity=False, model_components=None, scale_factor=1, **kwargs):
    """Plot one threeML ``XYLike`` dataset in frequency SED space."""
    if not is_flux_space_plugin(plugin):
        raise TypeError(
            "Flux-space plotting accepts XYLike-style plugins only. "
            "Plot OGIPLike data separately with plot_ogip_with_model()."
        )
    ax = make_sed_axes(ax, style=style, luminosity=luminosity)
    energy = np.asarray(plugin.x)
    flux = photon_flux_to_mjy(plugin.y, energy)
    y = kev_to_hz(energy) * flux * MJY_TO_CGS * scale_factor
    if luminosity:
        if model_components is None:
            raise ValueError("model_components is required for luminosity plots")
        y *= _luminosity_factor(model_components)
    yerr = None
    if getattr(plugin, "has_errors", False):
        yerr = kev_to_hz(energy) * photon_flux_to_mjy(plugin.yerr, energy) * MJY_TO_CGS * scale_factor
        if luminosity:
            yerr *= _luminosity_factor(model_components)
    return ax.errorbar(kev_to_hz(energy), y, yerr=yerr, fmt=style.data_marker,
                       ms=style.data_markersize, alpha=style.data_alpha,
                       color=color, label=label or plugin.name, **kwargs)


def plot_xray_file(path, *, ax=None, model_components=None, luminosity=False,
                   y_col=2, yerr_col=3, value_kind="nuFnu", style=DEFAULT_STYLE,
                   scale_factor=1, **kwargs):
    """Plot a binned X-ray SED table.

    The legacy LLAGN tables use lower/upper frequency edges followed by
    ``nu F_nu`` and its uncertainty, so ``value_kind="nuFnu"`` is the default.
    Set ``value_kind="flux_density_mjy"`` for a table whose value columns are
    instead flux densities in mJy.
    """
    ax = make_sed_axes(ax, style=style, luminosity=luminosity)
    data = np.genfromtxt(Path(path))
    if data.ndim != 2 or data.shape[1] <= max(y_col, yerr_col):
        raise ValueError(f"Expected a 2-D table with columns through {max(y_col, yerr_col)}: {path}")
    frequency = np.sqrt(data[:, 0] * data[:, 1])
    if value_kind == "nuFnu":
        scale = 1.0
    elif value_kind == "flux_density_mjy":
        scale = frequency * MJY_TO_CGS
    else:
        raise ValueError(
            "value_kind must be 'nuFnu' or 'flux_density_mjy', "
            f"not {value_kind!r}"
        )
    if luminosity:
        if model_components is None:
            raise ValueError("model_components is required for luminosity plots")
        scale *= _luminosity_factor(model_components)
    return ax.errorbar(frequency, data[:, y_col] * scale * scale_factor,
                       yerr=data[:, yerr_col] * scale * scale_factor, **kwargs)


def evaluate_bhjet_components(model_components, *, energy_range=(1e-9, 1e3), n_eval=2, rerun=True):
    """Evaluate BHJet detailed output and return its component cache, restoring state."""
    jet = model_components["jet"]
    previous_detailed = getattr(jet, "enable_detailed_output", False)
    previous_infosw = jet.infosw.value
    try:
        jet.enable_detailed_output = True
        jet.infosw.value = 2
        if rerun:
            jet._cached_params = None
        grid = np.logspace(*np.log10(energy_range), max(2, int(n_eval)))
        jet(grid)
        return jet._last_components
    finally:
        jet.enable_detailed_output = previous_detailed
        jet.infosw.value = previous_infosw


def plot_bhjet_component_data(component_data, *, ax=None, components=None,
                              style=DEFAULT_STYLE, luminosity=False,
                              distance_kpc=None, scale_factor=1,
                              flux_density=False, legend=False, **kwargs):
    """Plot BHJet detailed-output dictionaries without evaluating a model.

    ``component_data`` is the mapping returned by
    :func:`bhjet_plotting.preprocess_component_output` or stored in a fitted
    model's ``_last_components`` cache. Values are expected in Hz and mJy.
    Set ``flux_density=True`` for an ``F_nu`` plot; the default is ``nu F_nu``
    (or ``nu L_nu`` when ``luminosity=True``).
    """
    ax = make_sed_axes(ax, style=style, luminosity=luminosity)
    if luminosity and distance_kpc is None:
        raise ValueError("distance_kpc is required for raw-component luminosity plots")
    factor = 4 * np.pi * (distance_kpc * KPC_TO_CM) ** 2 if luminosity else 1
    selected = tuple(component_data) if components is None else components
    custom_label = kwargs.pop("label", None)
    for index, name in enumerate(selected):
        if name not in component_data:
            continue
        frequency = np.asarray(component_data[name]["energy"], dtype=float)
        flux = np.asarray(component_data[name]["flux"], dtype=float)
        valid = np.isfinite(frequency) & np.isfinite(flux) & (frequency > 0)
        if valid.sum() < 2:
            continue
        order = np.argsort(frequency[valid])
        spec = component_plot_kwargs(style, name, kwargs)
        if custom_label is not None:
            spec["label"] = custom_label if index == 0 else "_nolegend_"
        y = flux[valid] * MJY_TO_CGS * factor * scale_factor
        if not flux_density:
            y = frequency[valid] * y
        ax.plot(frequency[valid][order], y[order], **spec)
    if legend:
        ax.legend(fontsize=style.legend_fontsize)
    return ax


def plot_bhjet_components(model_components, *, ax=None, components=("presyn", "postsyn", "precom", "postcom"), style=DEFAULT_STYLE, luminosity=False, energy_range=(1e-9, 1e3), n_eval=2, rerun=True, scale_factor=1, **kwargs):
    """Evaluate BHJet then add selected radiative components to an SED axis."""
    detailed = evaluate_bhjet_components(
        model_components, energy_range=energy_range, n_eval=n_eval, rerun=rerun,
    )
    distance_kpc = model_components["jet"].dist.value if luminosity else None
    return plot_bhjet_component_data(
        detailed, ax=ax, components=components, style=style, luminosity=luminosity,
        distance_kpc=distance_kpc, scale_factor=scale_factor, **kwargs,
    )


def plot_sed(*, data=None, model_components=None, model_expressions=None, xray_path=None, ax=None, style=DEFAULT_STYLE, luminosity=False, title=None, energy_range=(1e-9, 1e3), n_points=1000, legend=True):
    """Compose a publication-ready SED from data, analytic expressions, and X-rays.

    ``model_expressions`` maps display labels to already-built callable models.
    Use ``setup_define_scripts.build_spectrum`` outside this module to construct
    expressions from your YAML component dictionary.
    """
    ax = make_sed_axes(ax, style=style, luminosity=luminosity, title=title)
    flux_data, _ = split_plot_data(data)
    for name, plugin in flux_data.items():
        plot_xylike(plugin, ax=ax, color=style.data_colors.get(name), style=style,
                    luminosity=luminosity, model_components=model_components)
    if xray_path:
        plot_xray_file(xray_path, ax=ax, model_components=model_components, luminosity=luminosity,
                       fmt="D", ms=style.data_markersize, color="C3", label="X-ray")
    energy = np.logspace(*np.log10(energy_range), n_points)
    factor = _luminosity_factor(model_components) if luminosity else 1
    for label, model in (model_expressions or {}).items():
        flux = photon_flux_to_mjy(model(energy), energy)
        ax.plot(kev_to_hz(energy), kev_to_hz(energy) * flux * MJY_TO_CGS * factor,
                linewidth=style.model_linewidth, label=label)
    if legend:
        ax.legend(fontsize=style.legend_fontsize)
    return ax


def plot_source_sed(*, data, model_components, ax=None, xray_path=None, color=None,
                    label=None, style=DEFAULT_STYLE, luminosity=True, scale_factor=1,
                    components=("total",), energy_range=(1e-9, 1e3)):
    """Overlay one source's data, optional X-ray table, and BHJet components.

    Designed for comparison figures: call it once per source on the same axis.
    ``scale_factor`` offsets a source vertically while preserving all uncertainties.
    """
    ax = make_sed_axes(ax, style=style, luminosity=luminosity)
    for dataset in data.values():
        if is_flux_space_plugin(dataset):
            plot_xylike(dataset, ax=ax, color=color, style=style, luminosity=luminosity,
                        model_components=model_components, scale_factor=scale_factor)
    if xray_path:
        plot_xray_file(xray_path, ax=ax, model_components=model_components,
                       luminosity=luminosity, scale_factor=scale_factor, fmt="D",
                       ms=style.data_markersize, color=color)
    plot_bhjet_components(model_components, ax=ax, components=components, style=style,
                          luminosity=luminosity, energy_range=energy_range,
                          scale_factor=scale_factor, color=color, label=label)
    return ax


def plot_density_comparison(samples_by_label, *, parameters=None, log_parameters=(),
                            bins=80, style=DEFAULT_STYLE, alpha=0.7):
    """Return one consistent posterior-density figure per shared parameter.

    ``samples_by_label`` maps a run label to a pandas DataFrame (or mapping of
    one-dimensional arrays).  This keeps loading sampler outputs separate from
    the visualization, and makes run comparisons reusable across notebooks.
    """
    if not samples_by_label:
        return {}
    log_parameters = set(log_parameters)
    columns = [set(values.columns if hasattr(values, "columns") else values)
               for values in samples_by_label.values()]
    shared = set.intersection(*columns)
    selected = sorted(shared) if parameters is None else [name for name in parameters if name in shared]
    figures = {}
    for name in selected:
        series = {}
        for label, values in samples_by_label.items():
            column = values[name].to_numpy() if hasattr(values[name], "to_numpy") else np.asarray(values[name])
            column = column[np.isfinite(column)]
            if name in log_parameters:
                column = np.log10(column[column > 0])
            if column.size:
                series[label] = column
        if not series:
            continue
        all_values = np.concatenate(list(series.values()))
        low, high = np.nanmin(all_values), np.nanmax(all_values)
        if low == high:
            low, high = low - 0.5, high + 0.5
        edges = np.linspace(low, high, bins + 1)
        centers = (edges[:-1] + edges[1:]) / 2
        _, ax = plt.subplots(figsize=style.figsize)
        for label, values in series.items():
            density, _ = np.histogram(values, bins=edges, density=True)
            ax.plot(centers, density, label=label, alpha=alpha, linewidth=style.model_linewidth)
        ax.set(xlabel=f"log10({name})" if name in log_parameters else name, ylabel="Posterior density")
        ax.legend(fontsize=style.legend_fontsize)
        ax.figure.tight_layout()
        figures[name] = ax.figure
    return figures


def plot_confidence_band(ax, x, bands, *, label=None, color=None, show_95=True,
                         show_68=True, alpha_95=None, alpha_68=None,
                         style=DEFAULT_STYLE, transform=None, **kwargs):
    """Draw percentile confidence bands from a ``{percentile: values}`` mapping.

    Supply ``transform`` when spectra need conversion before display, e.g.
    ``lambda y: frequency * y * 1e-26`` for a BHJet mJy spectrum.
    """
    x = np.asarray(x)
    convert = transform or (lambda values: values)
    alpha_95 = style.confidence_alpha_95 if alpha_95 is None else alpha_95
    alpha_68 = style.confidence_alpha_68 if alpha_68 is None else alpha_68
    handles = []
    if show_95 and 2.5 in bands and 97.5 in bands:
        handles.append(ax.fill_between(x, convert(np.asarray(bands[2.5])),
                                      convert(np.asarray(bands[97.5])), color=color,
                                      alpha=alpha_95, label=label, **kwargs))
    if show_68 and 16.0 in bands and 84.0 in bands:
        handles.append(ax.fill_between(x, convert(np.asarray(bands[16.0])),
                                      convert(np.asarray(bands[84.0])), color=color,
                                      alpha=alpha_68,
                                      label="_nolegend_" if label else None, **kwargs))
    return handles
