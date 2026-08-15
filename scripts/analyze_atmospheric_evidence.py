"""Render the source-graded molecular-evidence audit.

This script visualizes reported literature evidence. It does not perform an
atmospheric retrieval or convert repository diagnostics into detection sigma.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "atmospheric_evidence.csv"
OUTPUT = ROOT / "figures" / "molecular_evidence.png"

COLORS = {
    "reported detection": "#19b58a",
    "reported evidence": "#79c267",
    "repository diagnostic": "#a56eff",
    "reported non-detection": "#4aa8d8",
    "model space excluded": "#4aa8d8",
    "not detected": "#78889c",
    "not established here": "#78889c",
    "unconstrained": "#78889c",
    "model-dependent hint": "#e2a84a",
    "no evidence": "#78889c",
}


def main() -> None:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    height = max(3.2, 0.72 * len(rows) + 1.55)
    fig, ax = plt.subplots(figsize=(11.2, height), constrained_layout=True)
    fig.patch.set_facecolor("#071018")
    ax.set_facecolor("#071018")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.axis("off")
    for index, row in enumerate(reversed(rows)):
        y = index
        color = COLORS.get(row["status"], "#78889c")
        ax.scatter(0.035, y, s=150, color=color, edgecolor="#dce8f2", linewidth=0.7)
        ax.text(0.075, y + 0.14, row["species"], color="#f4f8fb", fontsize=11.5, weight="bold", va="center")
        ax.text(0.075, y - 0.18, row["status"], color=color, fontsize=9.5, va="center")
        ax.text(0.44, y + 0.08, row["evidence"], color="#f4f8fb", fontsize=10.2, va="center")
        ax.text(0.44, y - 0.20, row["basis"], color="#9eb0bf", fontsize=8.7, va="center")
        if y:
            ax.plot([0.02, 0.98], [y - 0.48, y - 0.48], color="#213241", lw=0.7)
    ax.text(0.02, len(rows) - 0.12, "SOURCE-GRADED ATMOSPHERIC EVIDENCE", color="#72d4f7", fontsize=10, weight="bold")
    ax.text(0.02, -0.68, "Reported literature evidence — not an independent retrieval", color="#9eb0bf", fontsize=8.8)
    fig.savefig(OUTPUT, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
