"""Timing-adjusted, limb-darkened transit fit to a real TESS light curve.

The NASA Exoplanet Archive period and historical epoch define the predicted
phase.  The TESS data then fit a bounded timing correction, radius ratio,
impact parameter, and local linear baseline.  A circular orbit and fixed,
representative quadratic TESS-band limb-darkening coefficients are adopted;
the result is a reproducible diagnostic fit, not a precision global retrieval.
"""

from __future__ import annotations

import csv
from pathlib import Path

from astropy.io import fits
import batman
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
DATA_FILE = DATA_DIR / "tess2022190063128-s0054-0000000243921117-0227-s_lc.fits"
PLANET = "WASP-80 b"
PERIOD_DAYS = 3.06785234
EPOCH_BJD = 2456487.425006
DURATION_HOURS = 2.131
A_OVER_RS = 12.6230389078
LIMB_DARKENING = (0.450, 0.250)
SECTOR = 54
FIGURE_FILE = FIG_DIR / "wasp80b_tess_transit.png"
COLOR = "#7048e8"


def load_light_curve(path: Path = DATA_FILE) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return quality-filtered, normalized BJD time, PDCSAP flux, and error."""
    with fits.open(path, memmap=False) as hdul:
        data = hdul[1].data; header = hdul[1].header
        time = np.asarray(data["TIME"], dtype=float)
        flux = np.asarray(data["PDCSAP_FLUX"], dtype=float)
        error = np.asarray(data["PDCSAP_FLUX_ERR"], dtype=float)
        quality = np.asarray(data["QUALITY"], dtype=int)
        bjdref = float(header.get("BJDREFI", 2457000.0)) + float(header.get("BJDREFF", 0.0))
    good = (quality == 0) & np.isfinite(time) & np.isfinite(flux) & np.isfinite(error) & (error > 0)
    time, flux, error = time[good] + bjdref, flux[good], error[good]
    scale = np.nanmedian(flux); flux, error = flux / scale, error / scale
    median = np.nanmedian(flux); mad_sigma = 1.4826 * np.nanmedian(np.abs(flux - median))
    threshold = max(8.0 * mad_sigma, 0.05); keep = np.abs(flux - median) <= threshold
    return time[keep], flux[keep], error[keep], int((~keep).sum())


def phase_offset_days(time_bjd: np.ndarray) -> np.ndarray:
    """Signed time from the nearest published mid-transit, in days."""
    return ((time_bjd - EPOCH_BJD + PERIOD_DAYS / 2.0) % PERIOD_DAYS) - PERIOD_DAYS / 2.0


def transit_profile(offset: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    """Evaluate a circular quadratic-limb-darkened transit and linear baseline."""
    timing, radius_ratio, impact, baseline, slope = parameters
    model = batman.TransitParams()
    model.t0 = timing; model.per = PERIOD_DAYS; model.rp = radius_ratio; model.a = A_OVER_RS
    model.inc = np.degrees(np.arccos(np.clip(impact / A_OVER_RS, 0, .999999)))
    model.ecc = 0.; model.w = 90.; model.u = list(LIMB_DARKENING); model.limb_dark = "quadratic"
    transit = batman.TransitModel(model, np.asarray(offset), supersample_factor=5,
                                 exp_time=2 / 1440).light_curve(model)
    return baseline + slope * (np.asarray(offset) - timing) + transit - 1


def weighted_linear_null(offset: np.ndarray, flux: np.ndarray, error: np.ndarray):
    design = np.column_stack([np.ones(len(offset)), offset])
    weights = 1 / error**2
    covariance = np.linalg.inv(design.T @ (weights[:, None] * design))
    coefficients = covariance @ (design.T @ (weights * flux))
    return design @ coefficients


def model_depth(parameters: np.ndarray) -> float:
    timing, _, _, baseline, _ = parameters
    return float((baseline - transit_profile(np.asarray([timing]), parameters)[0]) * 1e6)


def transit_duration_hours(parameters: np.ndarray) -> float:
    _, radius_ratio, impact, _, _ = parameters
    inclination = np.arccos(np.clip(impact / A_OVER_RS, 0, .999999))
    numerator = np.sqrt(max((1 + radius_ratio)**2 - impact**2, 0))
    argument = np.clip(numerator / (A_OVER_RS * np.sin(inclination)), 0, 1)
    return float(PERIOD_DAYS / np.pi * np.arcsin(argument) * 24)


def compare_models(time: np.ndarray, flux: np.ndarray, error: np.ndarray) -> dict[str, object]:
    """Compare a local linear null with a timing-adjusted transit profile."""
    offset = phase_offset_days(time); duration_days = DURATION_HOURS / 24
    near = np.abs(offset) <= 3.0 * duration_days
    x_time, y, sigma = offset[near], flux[near], error[near]
    if len(y) < 40:
        raise ValueError("The archived sector does not contain enough transit-window samples")
    grid = np.linspace(-1.25 * duration_days, 1.25 * duration_days, 151)
    scores = []
    for center in grid:
        core = np.abs(x_time - center) < .30 * duration_days
        scores.append(np.median(y[core]) if core.sum() >= 3 else np.inf)
    center0 = float(grid[int(np.argmin(scores))])
    core = np.abs(x_time - center0) < .30 * duration_days
    wings = np.abs(x_time - center0) > .70 * duration_days
    if core.sum() < 3 or wings.sum() < 20:
        raise ValueError("The archived sector does not sample a usable fitted transit window")
    depth0 = np.clip(np.median(y[wings]) - np.median(y[core]), 1e-5, .18)
    initial = np.asarray([center0, np.sqrt(depth0), .4, np.median(y[wings]), 0.])
    lower = np.asarray([-1.5 * duration_days, .002, 0., .94, -.08])
    upper = np.asarray([1.5 * duration_days, .48, .98, 1.06, .08])
    fit = least_squares(lambda pars: (y - transit_profile(x_time, pars)) / sigma,
                        initial, bounds=(lower, upper), x_scale="jac", max_nfev=700)
    parameters = fit.x; model = transit_profile(x_time, parameters)
    null_model = weighted_linear_null(x_time, y, sigma)
    chi2_null = float(np.sum(((y - null_model) / sigma)**2))
    chi2_transit = float(np.sum(((y - model) / sigma)**2))
    dof_null, dof_transit = len(y) - 2, len(y) - len(parameters)
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
    parameter_errors = np.sqrt(np.maximum(np.diag(covariance), 0))
    gradient = np.empty(len(parameters))
    for index, value in enumerate(parameters):
        step = max(abs(value) * 1e-5, 1e-7)
        plus, minus = parameters.copy(), parameters.copy()
        plus[index] += step; minus[index] -= step
        gradient[index] = (model_depth(plus) - model_depth(minus)) / (2 * step)
    depth = model_depth(parameters)
    formal_depth_error = float(np.sqrt(max(gradient @ covariance @ gradient, 0)))
    reduced_chi_square = chi2_transit / dof_transit
    scatter_scale = np.sqrt(max(reduced_chi_square, 1.0))
    depth_error = formal_depth_error * scatter_scale
    reported_parameter_errors = parameter_errors * scatter_scale
    fitted_duration = transit_duration_hours(parameters)
    in_transit = np.abs(x_time - parameters[0]) <= fitted_duration / 48
    delta = chi2_null - chi2_transit; delta_dof = len(parameters) - 2
    bic_null = chi2_null + 2 * np.log(len(y))
    bic_transit = chi2_transit + len(parameters) * np.log(len(y))
    return {
        "n_points": len(y), "n_in_transit": int(in_transit.sum()),
        "n_out_of_transit": int((~in_transit).sum()), "depth_ppm": depth,
        "depth_error_ppm": depth_error, "depth_snr": depth / depth_error if depth_error > 0 else np.nan,
        "formal_depth_error_ppm": formal_depth_error, "reduced_chi_square": reduced_chi_square,
        "radius_ratio": parameters[1], "radius_ratio_error": reported_parameter_errors[1],
        "impact_parameter": parameters[2], "impact_parameter_error": reported_parameter_errors[2],
        "timing_shift_hours": parameters[0] * 24, "timing_error_minutes": reported_parameter_errors[0] * 1440,
        "fitted_duration_hours": fitted_duration, "baseline": parameters[3], "baseline_slope_per_day": parameters[4],
        "chi_square_null": chi2_null, "dof_null": dof_null,
        "p_value_null": chi2.sf(chi2_null, dof_null),
        "chi_square_transit": chi2_transit, "dof_transit": dof_transit,
        "p_value_transit": chi2.sf(chi2_transit, dof_transit),
        "delta_chi_square": delta, "delta_dof": delta_dof,
        "p_value_improvement": chi2.sf(max(delta, 0), delta_dof),
        "bic_null": bic_null, "bic_transit": bic_transit, "delta_bic": bic_null - bic_transit,
        "transit_supported": bool(bic_null - bic_transit >= 10),
        # Compatibility aliases used by the multi-sector robustness module.
        "chi_square_flat": chi2_null, "dof_flat": dof_null, "p_value_flat": chi2.sf(chi2_null, dof_null),
        "chi_square_box": chi2_transit, "dof_box": dof_transit, "p_value_box": chi2.sf(chi2_transit, dof_transit),
        "offset": x_time, "flux": y, "error": sigma, "in_transit": in_transit,
        "model": model, "null_model": null_model, "parameters": parameters,
        "parameter_errors": parameter_errors, "fit_success": bool(fit.success),
    }


def binned_curve(offset: np.ndarray, flux: np.ndarray, error: np.ndarray, bins: int = 70):
    edges = np.linspace(offset.min(), offset.max(), bins + 1); centers, means, uncertainties = [], [], []
    for left, right in zip(edges[:-1], edges[1:]):
        chosen = (offset >= left) & (offset < right)
        if not chosen.any(): continue
        weights = 1 / error[chosen]**2
        centers.append((left + right) / 2); means.append(np.sum(weights * flux[chosen]) / np.sum(weights))
        uncertainties.append(np.sqrt(1 / np.sum(weights)))
    return np.asarray(centers), np.asarray(means), np.asarray(uncertainties)


def main() -> dict[str, object]:
    FIG_DIR.mkdir(exist_ok=True); time, flux, error, clipped = load_light_curve()
    result = compare_models(time, flux, error)
    rows = [("sector", SECTOR, "TESS sector"), ("period_days", PERIOD_DAYS, "days; NASA Exoplanet Archive"),
            ("published_duration_hours", DURATION_HOURS, "hours; NASA Exoplanet Archive"),
            ("a_over_rstar", A_OVER_RS, "fixed from saved a and stellar radius"),
            ("limb_darkening_u1", LIMB_DARKENING[0], "fixed quadratic coefficient"),
            ("limb_darkening_u2", LIMB_DARKENING[1], "fixed quadratic coefficient"),
            ("quality_filtered_cadences", len(time), "count"), ("symmetric_outliers_clipped", clipped, "count")]
    keys = ("n_points", "n_in_transit", "n_out_of_transit", "depth_ppm", "depth_error_ppm", "depth_snr",
            "formal_depth_error_ppm", "reduced_chi_square", "radius_ratio", "radius_ratio_error", "impact_parameter", "impact_parameter_error",
            "timing_shift_hours", "timing_error_minutes", "fitted_duration_hours",
            "chi_square_null", "dof_null", "p_value_null", "chi_square_transit", "dof_transit",
            "p_value_transit", "delta_chi_square", "delta_dof", "p_value_improvement",
            "bic_null", "bic_transit", "delta_bic", "transit_supported")
    rows.extend((key, result[key], "") for key in keys)
    with (FIG_DIR / "summary_statistics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle); writer.writerow(["quantity", "value", "unit"])
        for quantity, value, unit in rows:
            writer.writerow([quantity, f"{value:.12g}" if isinstance(value, float) else value, unit])
    centers, means, uncertainties = binned_curve(result["offset"], result["flux"], result["error"])
    dense = np.linspace(result["offset"].min(), result["offset"].max(), 1400)
    dense_model = transit_profile(dense, result["parameters"])
    fig, ax = plt.subplots(figsize=(9.6, 5.6))
    ax.scatter(result["offset"] * 24, (result["flux"] - 1) * 1e6, s=5, alpha=.12,
               color=COLOR, label="quality-filtered cadences")
    ax.errorbar(centers * 24, (means - 1) * 1e6, yerr=uncertainties * 1e6, fmt="o", ms=4,
                color="#17212b", ecolor="#52606d", linewidth=1, label="inverse-variance binned")
    supported = result["transit_supported"]
    ax.plot(dense * 24, (dense_model - 1) * 1e6, color=COLOR, linewidth=2.4,
            linestyle="-" if supported else "--",
            label="timing-adjusted limb-darkened fit" if supported else "best transit profile (not BIC-supported)")
    if not supported:
        order = np.argsort(result["offset"])
        ax.plot(result["offset"][order] * 24, (result["null_model"][order] - 1) * 1e6,
                color="#52606d", linewidth=1.8, label="preferred linear null")
    ax.axvline(0, color="#52606d", ls=":", lw=1.3, label="published prediction")
    ax.axvline(result["timing_shift_hours"], color=COLOR, ls="--", lw=1.3, alpha=.7,
               label=f"fitted midpoint ({result['timing_shift_hours']:+.2f} h)")
    ax.set(xlabel="Hours from published mid-transit", ylabel="Normalized flux - 1 [ppm]",
           title=f"{PLANET}: timing-adjusted TESS Sector {SECTOR} transit fit")
    ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=8); fig.tight_layout()
    fig.savefig(FIGURE_FILE, dpi=200); plt.close(fig)
    return result


if __name__ == "__main__":
    stats = main()
    print(f"{PLANET}: depth={stats['depth_ppm']:.1f} +/- {stats['depth_error_ppm']:.1f} ppm; "
          f"timing shift={stats['timing_shift_hours']:+.3f} h; supported={stats['transit_supported']}")
    print(f"transit chi2/dof={stats['chi_square_transit']:.1f}/{stats['dof_transit']}; "
          f"Delta BIC={stats['delta_bic']:.1f}")
