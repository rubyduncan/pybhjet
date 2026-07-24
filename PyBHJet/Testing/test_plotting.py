"""Lightweight tests for plotting helpers; no ThreeML or HEASoft required."""

import os
import sys
import unittest

os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import matplotlib.pyplot as plt
import numpy as np

from plotting import (
    DEFAULT_STYLE,
    component_plot_kwargs,
    plot_bhjet_component_data,
    plot_bhjet_components,
    plot_confidence_band,
    plot_sed,
    plot_xylike,
    split_plot_data,
)
from bhjet_plotting import plot_flux_mjy, plot_nufnu_ergshz


class FakeXYLike:
    name = "radio"
    has_errors = False
    x = np.array([1.0, 2.0, 4.0])
    y = np.array([1e-4, 5e-5, 2e-5])


class OGIPLike:
    pass


class Value:
    def __init__(self, value):
        self.value = value


class FakeJet:
    def __init__(self):
        self.enable_detailed_output = False
        self.infosw = Value(1)
        self.dist = Value(10.0)
        self._cached_params = "cached"
        self._last_components = {"total": {"energy": np.array([1e9, 1e10]), "flux": np.array([1.0, 0.5])}}

    def __call__(self, grid):
        self.grid = grid
        return grid


class PlottingTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_component_override_preserves_other_defaults(self):
        style = DEFAULT_STYLE.with_component_overrides(presyn={"color": "navy"})
        self.assertEqual(component_plot_kwargs(style, "presyn")["color"], "navy")
        self.assertEqual(component_plot_kwargs(style, "postcom")["color"], "green")
        self.assertEqual(component_plot_kwargs(style, "postsyn")["linestyle"], (0, (5, 1)))
        legacy_style = DEFAULT_STYLE.with_component_overrides(disk={"style": "--"})
        self.assertEqual(component_plot_kwargs(legacy_style, "disk")["linestyle"], "--")

    def test_raw_component_plot_and_confidence_band_render(self):
        components = {
            "presyn": {"energy": np.array([1e9, 1e10]), "flux": np.array([2.0, 1.0])},
            "disk": {"energy": np.array([1e10, 1e11]), "flux": np.array([1.0, 0.5])},
        }
        ax = plot_bhjet_component_data(components, legend=True)
        self.assertEqual(len(ax.lines), 2)
        self.assertEqual(ax.lines[0].get_color(), "dodgerblue")
        handles = plot_confidence_band(
            ax, np.array([1e9, 1e10]),
            {2.5: np.array([1.0, 1.0]), 16.0: np.array([1.5, 1.5]),
             84.0: np.array([2.0, 2.0]), 97.5: np.array([2.5, 2.5])},
        )
        self.assertEqual(len(handles), 2)
        ax.figure.canvas.draw()

    def test_legacy_component_helpers_return_figure_and_axes(self):
        components = {"total": {"energy": np.array([1e9, 1e10]), "flux": np.array([2.0, 1.0])}}
        figure, ax = plot_nufnu_ergshz(components)
        self.assertIs(figure, ax.figure)
        self.assertEqual(len(ax.lines), 1)
        figure, ax = plot_flux_mjy(components, title="Test")
        self.assertIs(figure, ax.figure)
        self.assertEqual(len(ax.lines), 1)

    def test_model_component_plotter_evaluates_and_restores_jet_state(self):
        jet = FakeJet()
        ax = plot_bhjet_components({"jet": jet}, components=("total",), n_eval=3)
        self.assertEqual(len(jet.grid), 3)
        self.assertFalse(jet.enable_detailed_output)
        self.assertEqual(jet.infosw.value, 1)
        self.assertEqual(len(ax.lines), 1)

    def test_ogip_is_never_sent_to_flux_sed(self):
        flux, ogip = split_plot_data({"radio": FakeXYLike(), "xray": OGIPLike()})
        self.assertEqual(set(flux), {"radio"})
        self.assertEqual(set(ogip), {"xray"})
        ax = plot_sed(data={"radio": FakeXYLike(), "xray": OGIPLike()})
        self.assertEqual(len(ax.lines), 1)
        with self.assertRaises(TypeError):
            plot_xylike(OGIPLike())


if __name__ == "__main__":
    unittest.main()
