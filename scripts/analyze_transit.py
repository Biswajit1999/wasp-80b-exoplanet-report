"""Transit-window analysis of a real TESS SPOC light curve for WASP-80 b.

The archived PDCSAP_FLUX series is read directly from the unmodified MAST FITS
file. Period, epoch, and duration are fixed to the saved NASA Exoplanet Archive
composite values. A constant model is compared with a two-level box model; the
box depth is fitted, but the transit timing and width are not searched.
"""

from __future__ import annotations

import csv
from pathlib import Path

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
DATA_FILE = DATA_DIR / "tess2022190063128-s0054-0000000243921117-0227-s_lc.fits"
PLANET = "WASP-80 b"
PERIOD_DAYS = 3.06785234
EPOCH_BJD = 2456487.425006
DURATION_HOURS = 2.131
SECTOR = 54
FIGURE_FILE = FIG_DIR / "wasp80b_tess_transit.png"


def load_light_curve(path: Path = DATA_FILE) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return quality-filtered, normalized BJD time, PDCSAP flux, and error."""
    with fits.open(path, memmap=False) as hdul:
        data = hdul[1].data
        header = hdul[1].header
        time = np.asarray(data["TIME"], dtype=float)
        flux = np.asarray(data["PDCSAP_FLUX"], dtype=float)
        error = np.asarray(data["PDCSAP_FLUX_ERR"], dtype=float)
        quality = np.asarray(data["QUALITY"], dtype=int)
        bjdref = float(header.get("BJDREFI", 2457000.0)) + float(header.get("BJDREFF", 0.0))
    good = (quality == 0) & np.isfinite(time) & np.isfinite(flux) & np.isfinite(error) & (error > 0)
    time, flux, error = time[good] + bjdref, flux[good], error[good]
    scale = np.nanmedian(flux)
    flux, error = flux / scale, error / scale
    median = np.nanmedian(flux)
    mad_sigma = 1.4826 * np.nanmedian(np.abs(flux - median))
    threshold = max(8.0 * mad_sigma, 0.05)
    keep = np.abs(flux - median) <= threshold
    return time[keep], flux[keep], error[keep], int((~keep).sum())


def phase_offset_days(time_bjd: np.ndarray) -> np.ndarray:
    """Signed time from the nearest published mid-transit, in days."""
    return ((time_bjd - EPOCH_BJD + PERIOD_DAYS / 2.0) % PERIOD_DAYS) - PERIOD_DAYS / 2.0


def compare_models(time: np.ndarray, flux: np.ndarray, error: np.ndarray) -> dict[str, float]:
    """Compare a weighted flat model to a fixed-window, fitted-depth box."""
    offset = phase_offset_days(time)
    duration_days = DURATION_HOURS / 24.0
    near = np.abs(offset) <= 2.5 * duration_days
    in_transit = np.abs(offset[near]) <= duration_days / 2.0
    x_time, y, sigma = offset[near], flux[near], error[near]
    if in_transit.sum() < 5 or (~in_transit).sum() < 20:
        raise ValueError("The archived sector does not contain enough transit-window samples")
    weights = 1.0 / sigma**2
    flat_level = np.sum(weights * y) / np.sum(weights)
    chi2_flat = np.sum(((y - flat_level) / sigma) ** 2)

    design = np.column_stack([np.ones(len(y)), in_transit.astype(float)])
    normal = design.T @ (weights[:, None] * design)
    covariance = np.linalg.inv(normal)
    coefficients = covariance @ (design.T @ (weights * y))
    model = design @ coefficients
    chi2_box = np.sum(((y - model) / sigma) ** 2)
    depth = -coefficients[1]
    depth_error = np.sqrt(covariance[1, 1])
    delta = chi2_flat - chi2_box
    return {
        "n_points": len(y),
        "n_in_transit": int(in_transit.sum()),
        "n_out_of_transit": int((~in_transit).sum()),
        "flat_level": flat_level,
        "box_baseline": coefficients[0],
        "depth_ppm": depth * 1e6,
        "depth_error_ppm": depth_error * 1e6,
        "depth_snr": depth / depth_error,
        "chi_square_flat": chi2_flat,
        "dof_flat": len(y) - 1,
        "p_value_flat": chi2.sf(chi2_flat, len(y) - 1),
        "chi_square_box": chi2_box,
        "dof_box": len(y) - 2,
        "p_value_box": chi2.sf(chi2_box, len(y) - 2),
        "delta_chi_square": delta,
        "delta_dof": 1,
        "p_value_improvement": chi2.sf(max(delta, 0.0), 1),
        "offset": x_time,
        "flux": y,
        "error": sigma,
        "in_transit": in_transit,
        "model": model,
    }


def binned_curve(offset: np.ndarray, flux: np.ndarray, error: np.ndarray, bins: int = 70):
    edges = np.linspace(offset.min(), offset.max(), bins + 1)
    centers, means, uncertainties = [], [], []
    for left, right in zip(edges[:-1], edges[1:]):
        chosen = (offset >= left) & (offset < right)
        if not chosen.any():
            continue
        weights = 1.0 / error[chosen] ** 2
        centers.append((left + right) / 2.0)
        means.append(np.sum(weights * flux[chosen]) / np.sum(weights))
        uncertainties.append(np.sqrt(1.0 / np.sum(weights)))
    return np.asarray(centers), np.asarray(means), np.asarray(uncertainties)


def main() -> dict[str, float]:
    FIG_DIR.mkdir(exist_ok=True)
    time, flux, error, clipped = load_light_curve()
    result = compare_models(time, flux, error)
    rows = [
        ("sector", SECTOR, "TESS sector"),
        ("period_days", PERIOD_DAYS, "days; NASA Exoplanet Archive"),
        ("duration_hours", DURATION_HOURS, "hours; NASA Exoplanet Archive"),
        ("quality_filtered_cadences", len(time), "count"),
        ("symmetric_outliers_clipped", clipped, "count"),
    ]
    for key in (
        "n_points", "n_in_transit", "n_out_of_transit", "depth_ppm",
        "depth_error_ppm", "depth_snr", "chi_square_flat", "dof_flat",
        "p_value_flat", "chi_square_box", "dof_box", "p_value_box",
        "delta_chi_square", "delta_dof", "p_value_improvement",
    ):
        rows.append((key, result[key], ""))
    with (FIG_DIR / "summary_statistics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        for quantity, value, unit in rows:
            writer.writerow([quantity, f"{value:.12g}" if isinstance(value, float) else value, unit])

    centers, means, uncertainties = binned_curve(result["offset"], result["flux"], result["error"])
    hours = result["offset"] * 24.0
    fig, ax = plt.subplots(figsize=(9.6, 5.6))
    ax.scatter(hours, (result["flux"] - 1) * 1e6, s=5, alpha=0.12, color="#7048e8", label="quality-filtered cadences")
    ax.errorbar(centers * 24.0, (means - 1) * 1e6, yerr=uncertainties * 1e6,
                fmt="o", ms=4, color="#17212b", ecolor="#52606d", linewidth=1,
                label="inverse-variance binned")
    order = np.argsort(hours)
    ax.plot(hours[order], (result["model"][order] - 1) * 1e6,
            color="#7048e8", linewidth=2.2, label="fitted box model")
    ax.axvspan(-DURATION_HOURS / 2.0, DURATION_HOURS / 2.0, color="#7048e8", alpha=0.08)
    ax.set(xlabel="Hours from published mid-transit", ylabel="Normalized flux − 1 [ppm]",
           title=f"{PLANET}: real TESS Sector {SECTOR} PDCSAP light curve")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_FILE, dpi=200)
    plt.close(fig)
    return result


if __name__ == "__main__":
    stats = main()
    print(f"{PLANET}: depth={stats['depth_ppm']:.1f} +/- {stats['depth_error_ppm']:.1f} ppm")
    print(f"flat chi2/dof={stats['chi_square_flat']:.1f}/{stats['dof_flat']}; p={stats['p_value_flat']:.4g}")
    print(f"box chi2/dof={stats['chi_square_box']:.1f}/{stats['dof_box']}; p={stats['p_value_box']:.4g}")
    print(f"Delta chi2={stats['delta_chi_square']:.1f} for 1 dof; p={stats['p_value_improvement']:.4g}")
