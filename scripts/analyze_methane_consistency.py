"""Quantify whether the published WASP-80 b methane morphology is robust.

This is a transparent model-shape diagnostic using Bell et al. (2023) source
data. It is not an atmospheric retrieval, and delta chi-square values are not
converted into detection significances because the supplied curves are
posterior summaries rather than nested maximum-likelihood fits.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPECTRUM = ROOT / "data" / "spectra" / "wasp80b_nircam_source_data_rebinned.csv"
ABUNDANCES = ROOT / "data" / "published_methane_abundances.csv"
SYSTEM = ROOT / "data" / "system_parameters.csv"
FIGURE = ROOT / "figures" / "wasp80b_methane_consistency.png"
STATS = ROOT / "figures" / "methane_consistency_statistics.csv"
CONTRIBUTIONS = ROOT / "figures" / "methane_bin_contributions.csv"

MODEL_COLUMNS = {
    "transmission": {
        "full": ("trans_full_lower", "trans_full_upper"),
        "no_ch4": ("trans_no_ch4_lower", "trans_no_ch4_upper"),
        "no_h2o": ("trans_no_h2o_lower", "trans_no_h2o_upper"),
    },
    "emission": {
        "full": ("emission_full_lower", "emission_full_upper"),
        "no_ch4": ("emission_no_ch4_lower", "emission_no_ch4_upper"),
        "no_h2o": ("emission_no_h2o_lower", "emission_no_h2o_upper"),
    },
}


def read_numeric_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {key: np.asarray([float(row[key]) for row in rows]) for key in rows[0]}


def fitted_offset_chi2(y: np.ndarray, error: np.ndarray, model: np.ndarray) -> tuple[float, float]:
    """Return chi-square and the fitted additive data-minus-model offset."""
    weight = 1.0 / error**2
    offset = float(np.sum(weight * (y - model)) / np.sum(weight))
    chi2 = float(np.sum(((y - (model + offset)) / error) ** 2))
    return chi2, offset


def combine_bins(values: np.ndarray, errors: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray]:
    combined, uncertainty = [], []
    for start in range(0, len(values), factor):
        part = slice(start, min(start + factor, len(values)))
        weight = 1.0 / errors[part] ** 2
        combined.append(np.sum(values[part] * weight) / np.sum(weight))
        uncertainty.append(1.0 / np.sqrt(np.sum(weight)))
    return np.asarray(combined), np.asarray(uncertainty)


def rebinned_delta_chi2(y: np.ndarray, error: np.ndarray, full: np.ndarray, no_ch4: np.ndarray, factor: int) -> float:
    yb, eb = combine_bins(y, error, factor)
    fb, _ = combine_bins(full, error, factor)
    nb, _ = combine_bins(no_ch4, error, factor)
    return fitted_offset_chi2(yb, eb, nb)[0] - fitted_offset_chi2(yb, eb, fb)[0]


def geometry_diagnostics(y: np.ndarray, error: np.ndarray, models: dict[str, np.ndarray]) -> dict[str, object]:
    fits = {name: fitted_offset_chi2(y, error, model) for name, model in models.items()}
    full_chi2, full_offset = fits["full"]
    no_ch4_chi2, no_ch4_offset = fits["no_ch4"]
    full_resid = ((y - (models["full"] + full_offset)) / error) ** 2
    no_ch4_resid = ((y - (models["no_ch4"] + no_ch4_offset)) / error) ** 2
    per_bin = no_ch4_resid - full_resid
    jackknife = []
    for omitted in range(len(y)):
        keep = np.arange(len(y)) != omitted
        jackknife.append(
            fitted_offset_chi2(y[keep], error[keep], models["no_ch4"][keep])[0]
            - fitted_offset_chi2(y[keep], error[keep], models["full"][keep])[0]
        )
    return {
        "fits": fits,
        "delta": no_ch4_chi2 - full_chi2,
        "per_bin": per_bin,
        "jackknife": np.asarray(jackknife),
        "rebin": {factor: rebinned_delta_chi2(y, error, models["full"], models["no_ch4"], factor) for factor in (1, 2, 4)},
    }


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    data = read_numeric_csv(SPECTRUM)
    wavelength = data["wavelength_micron"]
    models = {
        geometry: {name: 0.5 * (data[lower] + data[upper]) for name, (lower, upper) in columns.items()}
        for geometry, columns in MODEL_COLUMNS.items()
    }
    observations = {
        "transmission": (data["transit_depth"], data["transit_error"]),
        "emission": (data["eclipse_depth"], data["eclipse_error"]),
    }
    diagnostics = {
        geometry: geometry_diagnostics(*observations[geometry], models[geometry])
        for geometry in observations
    }

    with ABUNDANCES.open(newline="", encoding="utf-8") as handle:
        abundance = list(csv.DictReader(handle))
    transmission_abundance, emission_abundance = abundance
    abundance_difference = float(emission_abundance["log10_ch4_mixing_ratio"]) - float(transmission_abundance["log10_ch4_mixing_ratio"])
    transmission_sigma = np.mean([float(transmission_abundance["lower_uncertainty_dex"]), float(transmission_abundance["upper_uncertainty_dex"])])
    emission_sigma = np.mean([float(emission_abundance["lower_uncertainty_dex"]), float(emission_abundance["upper_uncertainty_dex"])])
    abundance_tension = abundance_difference / np.hypot(transmission_sigma, emission_sigma)

    with SYSTEM.open(newline="", encoding="utf-8") as handle:
        system = next(csv.DictReader(handle))
    planet_mass = float(system["pl_bmasse"]) * 5.9722e24
    planet_radius = float(system["pl_rade"]) * 6.371e6
    gravity = 6.67430e-11 * planet_mass / planet_radius**2
    scale_height = 1.380649e-23 * float(system["pl_eqt"]) / (2.3 * 1.66053906660e-27 * gravity)
    stellar_radius = float(system["st_rad"]) * 6.957e8
    radius_ratio = planet_radius / stellar_radius
    one_scale_height_ppm = 2.0 * radius_ratio * scale_height / stellar_radius * 1e6
    methane_window = (wavelength >= 3.05) & (wavelength <= 3.65)
    transmission_shape_ppm = np.ptp(models["transmission"]["full"][methane_window] - models["transmission"]["no_ch4"][methane_window]) * 1e6

    stats_rows: list[dict[str, object]] = []
    for geometry in observations:
        result = diagnostics[geometry]
        for model_name, (chi2, offset) in result["fits"].items():
            stats_rows.extend([
                {"metric": f"{geometry}_{model_name}_chi2", "value": f"{chi2:.6f}", "unit": "chi2", "interpretation": "one fitted additive offset"},
                {"metric": f"{geometry}_{model_name}_offset_ppm", "value": f"{offset * 1e6:.6f}", "unit": "ppm", "interpretation": "data minus supplied model midpoint"},
            ])
        stats_rows.extend([
            {"metric": f"{geometry}_delta_chi2_no_ch4_minus_full", "value": f"{result['delta']:.6f}", "unit": "chi2", "interpretation": "shape diagnostic; not detection sigma"},
            {"metric": f"{geometry}_jackknife_min_delta_chi2", "value": f"{np.min(result['jackknife']):.6f}", "unit": "chi2", "interpretation": "minimum after omitting one rebinned channel"},
            {"metric": f"{geometry}_jackknife_max_delta_chi2", "value": f"{np.max(result['jackknife']):.6f}", "unit": "chi2", "interpretation": "maximum after omitting one rebinned channel"},
        ])
        for factor, value in result["rebin"].items():
            stats_rows.append({"metric": f"{geometry}_rebin_{factor}_delta_chi2", "value": f"{value:.6f}", "unit": "chi2", "interpretation": "additional inverse-variance binning"})
    stats_rows.extend([
        {"metric": "published_log_ch4_difference_emission_minus_transmission", "value": f"{abundance_difference:.6f}", "unit": "dex", "interpretation": "published posterior medians"},
        {"metric": "published_log_ch4_tension", "value": f"{abundance_tension:.6f}", "unit": "sigma", "interpretation": "symmetrized independent-error comparison"},
        {"metric": "surface_gravity", "value": f"{gravity:.6f}", "unit": "m s-2", "interpretation": "from saved composite mass and radius"},
        {"metric": "scale_height_mu_2p3", "value": f"{scale_height / 1000:.6f}", "unit": "km", "interpretation": "isothermal at saved equilibrium temperature"},
        {"metric": "one_scale_height", "value": f"{one_scale_height_ppm:.6f}", "unit": "ppm", "interpretation": "transmission signal approximation"},
        {"metric": "methane_window_model_shape", "value": f"{transmission_shape_ppm:.6f}", "unit": "ppm", "interpretation": "peak-to-peak full minus no-CH4 model difference, 3.05-3.65 micron"},
    ])
    write_rows(STATS, ["metric", "value", "unit", "interpretation"], stats_rows)

    contribution_rows = []
    for index, wave in enumerate(wavelength):
        contribution_rows.append({
            "wavelength_micron": f"{wave:.9f}",
            "transmission_delta_chi2_contribution": f"{diagnostics['transmission']['per_bin'][index]:.9f}",
            "emission_delta_chi2_contribution": f"{diagnostics['emission']['per_bin'][index]:.9f}",
            "in_methane_window": str(bool(methane_window[index])).lower(),
        })
    write_rows(CONTRIBUTIONS, list(contribution_rows[0]), contribution_rows)

    colors = {"full": "#8e44ad", "no_ch4": "#f39c12", "no_h2o": "#2d9cdb"}
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6), constrained_layout=True)
    for ax, geometry, scale, ylabel in [
        (axes[0, 0], "transmission", 100, "Transit depth (%)"),
        (axes[0, 1], "emission", 1e6, "Eclipse depth (ppm)"),
    ]:
        y, error = observations[geometry]
        ax.errorbar(wavelength, y * scale, error * scale, fmt="o", ms=4, color="#17202a", ecolor="#85929e", label="Bell et al. source data")
        for name, model in models[geometry].items():
            offset = diagnostics[geometry]["fits"][name][1]
            ax.plot(wavelength, (model + offset) * scale, lw=2, color=colors[name], label=name.replace("_", " "))
        ax.axvspan(3.05, 3.65, color="#8e44ad", alpha=0.09)
        ax.set(xlabel="Wavelength (micron)", ylabel=ylabel, title=f"{geometry.title()} spectrum")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)

    ax = axes[1, 0]
    ax.axhline(0, color="#566573", lw=1)
    ax.bar(wavelength - 0.009, diagnostics["transmission"]["per_bin"], width=0.016, color="#8e44ad", label="transmission")
    ax.bar(wavelength + 0.009, diagnostics["emission"]["per_bin"], width=0.016, color="#2d9cdb", label="emission")
    ax.axvspan(3.05, 3.65, color="#8e44ad", alpha=0.09)
    ax.set(xlabel="Wavelength (micron)", ylabel="Per-bin delta chi-square", title="Where the no-CH4 mismatch accumulates")
    ax.legend()
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1, 1]
    factors = np.asarray([1, 2, 4])
    width = 0.34
    for shift, geometry, color in [(-width / 2, "transmission", "#8e44ad"), (width / 2, "emission", "#2d9cdb")]:
        values = [diagnostics[geometry]["rebin"][int(factor)] for factor in factors]
        ax.bar(factors.astype(float) + shift, values, width=width, color=color, label=geometry)
        jackknife = diagnostics[geometry]["jackknife"]
        ax.hlines(np.min(jackknife), 0.55, 4.45, color=color, ls=":", lw=1.5)
    ax.set_xticks(factors, ["native", "2x", "4x"])
    ax.set(xlabel="Additional binning", ylabel="delta chi-square (no CH4 - full)", title="Rebinning and leave-one-bin robustness")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle("WASP-80 b: methane morphology appears in both viewing geometries", fontsize=15, weight="bold")
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(f"Wrote {FIGURE}")
    print(f"Transmission delta chi2: {diagnostics['transmission']['delta']:.2f}")
    print(f"Emission delta chi2: {diagnostics['emission']['delta']:.2f}")
    print(f"Published abundance tension: {abundance_tension:.2f} sigma")


if __name__ == "__main__":
    main()
