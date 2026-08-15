"""Multi-sector transit-depth consistency and correlated-noise diagnostics.

The fixed NASA ephemeris used by ``analyze_transit`` is applied independently
to every committed standard-cadence SPOC light curve.  Formal depth errors are
inflated by both sqrt(reduced chi-square) and a time-averaging beta factor; no
sector is forced to agree with another.
"""

from __future__ import annotations

import csv
from pathlib import Path

from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

import analyze_transit as base


STATS_FILE = base.FIG_DIR / "multisector_statistics.csv"
STEM = base.FIGURE_FILE.stem.replace("_tess_transit", "")
TRANSITS_FIGURE = base.FIG_DIR / f"{STEM}_multisector_transits.png"
CONSISTENCY_FIGURE = base.FIG_DIR / f"{STEM}_depth_consistency.png"
NOISE_FIGURE = base.FIG_DIR / f"{STEM}_noise_diagnostics.png"


def sector_number(path: Path) -> int:
    with fits.open(path, memmap=False) as hdul:
        return int(hdul[0].header.get("SECTOR", hdul[1].header.get("SECTOR", -1)))


def rms_binned(residuals: np.ndarray, size: int) -> float:
    count = len(residuals) // size
    if count < 2:
        return float("nan")
    means = residuals[:count * size].reshape(count, size).mean(axis=1)
    return float(np.std(means, ddof=1))


def noise_curve(residuals: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    sizes = np.asarray([1, 2, 5, 10, 20, 30], dtype=int)
    measured = np.asarray([rms_binned(residuals, int(size)) for size in sizes])
    unbinned = measured[0]
    expected = []
    for size in sizes:
        groups = len(residuals) // size
        correction = np.sqrt(groups / (groups - 1)) if groups > 1 else np.nan
        expected.append(unbinned / np.sqrt(size) * correction)
    expected = np.asarray(expected)
    valid = (sizes >= 5) & np.isfinite(measured) & np.isfinite(expected) & (expected > 0)
    beta = float(max(1.0, np.max(measured[valid] / expected[valid]))) if valid.any() else 1.0
    return sizes, measured, expected, beta


def analyze_file(path: Path) -> dict[str, object]:
    time, flux, error, clipped = base.load_light_curve(path)
    result = base.compare_models(time, flux, error)
    residuals = np.asarray(result["flux"] - result["model"])
    sizes, measured, expected, beta = noise_curve(residuals)
    reduced_chi2 = float(result["chi_square_box"] / result["dof_box"])
    scale = np.sqrt(max(reduced_chi2, 1.0)) * beta
    robust_error = float(result["depth_error_ppm"] * scale)
    near = np.abs(base.phase_offset_days(time)) <= 2.5 * base.DURATION_HOURS / 24.0
    events = len(np.unique(np.rint((time[near] - base.EPOCH_BJD) / base.PERIOD_DAYS).astype(int)))
    return {"path": path, "sector": sector_number(path), "clipped": clipped,
            "events": events, "result": result, "residuals": residuals,
            "sizes": sizes, "measured": measured, "expected": expected,
            "beta": beta, "reduced_chi2": reduced_chi2,
            "robust_error_ppm": robust_error}


def main() -> dict[str, object]:
    base.FIG_DIR.mkdir(exist_ok=True)
    files = sorted(base.DATA_DIR.glob("tess*_lc.fits"), key=sector_number)
    if not files:
        raise FileNotFoundError("No committed TESS SPOC light curves were found")
    sectors, skipped = [], []
    for path in files:
        try:
            sectors.append(analyze_file(path))
        except ValueError as error:
            skipped.append({"path": path, "sector": sector_number(path), "reason": str(error)})
    if not sectors:
        raise ValueError("None of the committed sectors samples both sides of the fixed transit window")
    depths = np.asarray([item["result"]["depth_ppm"] for item in sectors])
    errors = np.asarray([item["robust_error_ppm"] for item in sectors])
    weights = 1.0 / errors**2
    combined = float(np.sum(weights * depths) / np.sum(weights))
    combined_error = float(np.sqrt(1.0 / np.sum(weights)))
    q = float(np.sum(weights * (depths - combined) ** 2))
    q_dof = len(sectors) - 1
    q_p = float(chi2.sf(q, q_dof)) if q_dof else float("nan")

    with STATS_FILE.open("w", newline="", encoding="utf-8") as handle:
        fields = ["sector", "filename", "transit_events", "n_points", "depth_ppm",
                  "formal_error_ppm", "reduced_chi_square", "beta", "robust_error_ppm"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in sectors:
            result = item["result"]
            writer.writerow({"sector": item["sector"], "filename": item["path"].name,
                             "transit_events": item["events"], "n_points": result["n_points"],
                             "depth_ppm": f"{result['depth_ppm']:.10g}",
                             "formal_error_ppm": f"{result['depth_error_ppm']:.10g}",
                             "reduced_chi_square": f"{item['reduced_chi2']:.10g}",
                             "beta": f"{item['beta']:.10g}",
                             "robust_error_ppm": f"{item['robust_error_ppm']:.10g}"})

    fig, axes = plt.subplots(len(sectors), 1, figsize=(9.6, 3.5 * len(sectors)), squeeze=False)
    for ax, item in zip(axes[:, 0], sectors):
        result = item["result"]
        centers, means, uncertainties = base.binned_curve(result["offset"], result["flux"], result["error"], bins=55)
        ax.errorbar(centers * 24, (means - 1) * 1e6, yerr=uncertainties * 1e6,
                    fmt="o", ms=3.5, color="#17212b", ecolor="#78909c", alpha=.9)
        order = np.argsort(result["offset"])
        ax.plot(result["offset"][order] * 24, (result["model"][order] - 1) * 1e6,
                color="#0b7285", lw=2)
        ax.set_title(f"Sector {item['sector']}: {result['depth_ppm']:.0f} +/- {item['robust_error_ppm']:.0f} ppm (scaled)", loc="left")
        ax.set_ylabel("Flux - 1 [ppm]")
        ax.grid(alpha=.2)
    axes[-1, 0].set_xlabel("Hours from fixed published mid-transit")
    fig.suptitle(f"{base.PLANET}: independent fixed-ephemeris TESS sector fits", y=.998, weight="bold")
    fig.tight_layout()
    fig.savefig(TRANSITS_FIGURE, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.9))
    labels = [f"S{item['sector']}" for item in sectors]
    positions = np.arange(len(sectors))
    ax.errorbar(positions, depths, yerr=errors, fmt="o", ms=7, capsize=4,
                color="#0b7285", ecolor="#52606d", label="sector depth (scaled error)")
    ax.axhspan(combined - combined_error, combined + combined_error, color="#0b7285", alpha=.14)
    ax.axhline(combined, color="#0b7285", lw=1.8, label="inverse-variance combination")
    ax.set(xticks=positions, xticklabels=labels, ylabel="Box depth [ppm]",
           title=f"{base.PLANET}: depth consistency across {len(sectors)} sector(s)")
    note = f"combined = {combined:.0f} +/- {combined_error:.0f} ppm"
    note += f"\nCochran Q = {q:.2f}, p = {q_p:.3g}" if q_dof else "\nQ test requires two or more sectors"
    ax.text(.02, .03, note, transform=ax.transAxes, va="bottom", fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": .8, "edgecolor": "#dce3e8"})
    ax.grid(axis="y", alpha=.2)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(CONSISTENCY_FIGURE, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    for item in sectors:
        ax.plot(item["sizes"], item["measured"] * 1e6, "o-", label=f"S{item['sector']} measured; beta={item['beta']:.2f}")
        ax.plot(item["sizes"], item["expected"] * 1e6, "--", alpha=.55)
    ax.set(xscale="log", yscale="log", xlabel="Cadences per bin", ylabel="Residual RMS [ppm]",
           title=f"{base.PLANET}: time-averaging residual-noise diagnostic")
    ax.grid(which="both", alpha=.2)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(NOISE_FIGURE, dpi=180)
    plt.close(fig)

    return {"sectors": sectors, "skipped": skipped, "combined_depth_ppm": combined,
            "combined_error_ppm": combined_error, "q": q, "q_dof": q_dof, "q_p": q_p}


if __name__ == "__main__":
    summary = main()
    print(f"{base.PLANET}: {len(summary['sectors'])} sector(s); robust combined depth "
          f"{summary['combined_depth_ppm']:.1f} +/- {summary['combined_error_ppm']:.1f} ppm")
