"""Regression tests for the timing-adjusted real-data transit fit."""

import json
import numpy as np
import analyze_transit as analysis


def test_real_fits_loads_and_physical_fit_converges():
    time, flux, error, clipped = analysis.load_light_curve()
    assert len(time) > 100 and np.all(np.isfinite(flux)) and np.all(error > 0)
    result = analysis.compare_models(time, flux, error)
    assert result["fit_success"]
    assert result["n_in_transit"] >= 3 and result["n_out_of_transit"] >= 20
    assert result["dof_null"] == result["n_points"] - 2
    assert result["dof_transit"] == result["n_points"] - 5
    assert 0.002 <= result["radius_ratio"] <= .48
    assert 0 <= result["impact_parameter"] <= .98
    assert abs(result["timing_shift_hours"]) <= 1.5 * analysis.DURATION_HOURS + 1e-6


def test_pipeline_reproduces_frozen_fit():
    result = analysis.main()
    expected = json.loads((analysis.ROOT / "tests" / "expected_fit.json").read_text(encoding="utf-8"))
    for key, tolerance in (("depth_ppm", .2), ("timing_shift_hours", 2e-4),
                           ("chi_square_transit", 1e-3), ("delta_bic", 1e-3)):
        assert np.isclose(result[key], expected[key], rtol=0, atol=tolerance)
    assert analysis.FIGURE_FILE.stat().st_size > 10_000
