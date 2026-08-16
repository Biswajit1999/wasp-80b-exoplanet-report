import csv
import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analyze_methane_consistency.py"


def load_module():
    spec = importlib.util.spec_from_file_location("methane", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_stats():
    with (ROOT / "figures" / "methane_consistency_statistics.csv").open(newline="", encoding="utf-8") as handle:
        return {row["metric"]: float(row["value"]) for row in csv.DictReader(handle)}


def test_offset_fit_recovers_known_shift():
    module = load_module()
    model = np.asarray([1.0, 2.0, 3.0])
    chi2, offset = module.fitted_offset_chi2(model + 0.25, np.ones(3), model)
    assert np.isclose(offset, 0.25)
    assert np.isclose(chi2, 0.0)


def test_spectrum_is_public_data_not_placeholder():
    path = ROOT / "data" / "spectra" / "wasp80b_nircam_source_data_rebinned.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 23
    wavelength = np.asarray([float(row["wavelength_micron"]) for row in rows])
    assert np.all(np.diff(wavelength) > 0)
    assert wavelength[0] < 2.5 and wavelength[-1] > 4.0


def test_methane_shape_diagnostic_survives_stress_tests():
    stats = read_stats()
    assert stats["transmission_delta_chi2_no_ch4_minus_full"] > 20
    assert stats["emission_delta_chi2_no_ch4_minus_full"] > 20
    assert stats["transmission_jackknife_min_delta_chi2"] > 10
    assert stats["emission_jackknife_min_delta_chi2"] > 10
    assert stats["transmission_rebin_4_delta_chi2"] > 10
    assert stats["emission_rebin_4_delta_chi2"] > 10


def test_published_abundances_are_consistent_between_geometries():
    stats = read_stats()
    assert abs(stats["published_log_ch4_difference_emission_minus_transmission"] - 0.3) < 1e-6
    assert abs(stats["published_log_ch4_tension"]) < 1.0
    assert stats["one_scale_height"] > 0
