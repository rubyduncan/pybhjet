# BHJet plotting interface

`PyBHJet/plotting.py` is the canonical plotting layer. It deliberately keeps
model/YAML construction in `setup_define_scripts.py`, so figures can be reused
from notebooks, scripts, and fitted-model output.

## Change one component

The default colours, labels, and dash patterns live in
`plotting.DEFAULT_COMPONENT_STYLES`. For a one-off notebook figure, copy the
default style and override only the component you want:

```python
from plotting import DEFAULT_STYLE, plot_bhjet_components

style = DEFAULT_STYLE.with_component_overrides(
    presyn={"color": "navy", "label": "Compact-jet synchrotron"},
    disk={"style": "--"},  # ``style`` and ``linestyle`` both work
)
plot_bhjet_components(components, ax=ax, style=style)
```

## Reuse a paper style

Create a named style in the notebook or figure script, then pass it to every
plotting call. No global Matplotlib settings are changed.

```python
PAPER_STYLE = DEFAULT_STYLE.with_overrides(
    figsize=(8, 5), model_linewidth=2.0, legend_fontsize=8
).with_component_overrides(
    total={"linewidth": 2.5},
    bb={"color": "darkorange"},
)

ax = plot_sed_from_components(data, components, shown_models, style=PAPER_STYLE)
plot_bhjet_components(components, ax=ax, style=PAPER_STYLE, components=("presyn", "postsyn", "precom", "postcom"))
```

## Flux-space versus OGIP plots

`plot_sed` and `plot_sed_from_components` plot only flux-space (`XYLike`)
data. OGIP spectra are response-folded count-space data and stay separate:

```python
flux_data, ogip_data = split_plot_data(data_dict)
ax = plot_sed(data=flux_data, model_components=components)
ogip_figure = plot_ogip_with_model(*ogip_data.values(), model_obj=model)
```

The legacy notebook helpers (`hz_plot_model_space_sed`,
`add_bhjet_radiative_components_to_plot`, `plot_nufnu_ergshz`, and
`plot_flux_mjy`) remain available and delegate to the same implementation.
