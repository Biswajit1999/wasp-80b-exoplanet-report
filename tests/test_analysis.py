"""Regression tests against the committed real MAST light curve."""

import csv

import numpy as np
import analyze_transit as analysis


def read_summary():
    with (analysis.FIG_DIR / "summary_statistics.csv").open(encoding="utf-8") as handle:
        return {row["quantity"]: float(row["value"]) for row in csv.DictReader(handle)}


def test_real_fits_loads_and_contains_transit_windows():
    time, flux, error, clipped = analysis.load_light_curve()
    assert len(time) > 100
    assert np.all(np.isfinite(flux))
    assert np.all(error > 0)
    result = analysis.compare_models(time, flux, error)
    assert result["n_in_transit"] >= 5
    assert result["n_out_of_transit"] >= 20
    assert result["dof_flat"] == result["n_points"] - 1
    assert result["dof_box"] == result["n_points"] - 2


def test_pipeline_reproduces_committed_statistics():
    analysis.main()
    values = read_summary()
    assert np.isclose(values["depth_ppm"], 25866.426518, rtol=0, atol=0.1)
    assert np.isclose(values["chi_square_flat"], 133043.495216, rtol=1e-9)
    assert np.isclose(values["chi_square_box"], 27276.5416366, rtol=1e-9)
    assert np.isclose(values["delta_chi_square"], 105766.95358, rtol=1e-9)
