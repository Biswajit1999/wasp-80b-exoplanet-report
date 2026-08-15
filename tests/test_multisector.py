"""Scientific checks for the committed multi-sector robustness analysis."""

from pathlib import Path

import numpy as np
import analyze_multisector as multi


def test_all_committed_spoc_files_are_accounted_for():
    expected = list(multi.base.DATA_DIR.glob("tess*_lc.fits"))
    summary = multi.main()
    assert len(summary["sectors"]) + len(summary["skipped"]) == len(expected) >= 1
    assert len(summary["sectors"]) >= 1
    if summary["supported"]:
        assert np.isfinite(summary["combined_depth_ppm"])
        assert summary["combined_error_ppm"] > 0
    else:
        assert np.isnan(summary["combined_depth_ppm"])
    for item in summary["sectors"]:
        assert item["result"]["n_in_transit"] >= 5
        assert item["result"]["n_out_of_transit"] >= 20
        assert item["beta"] >= 1
        assert item["robust_error_ppm"] >= item["result"]["depth_error_ppm"] - 1e-9


def test_reproducible_outputs_are_nonempty():
    multi.main()
    for path in (multi.STATS_FILE, multi.TRANSITS_FIGURE,
                 multi.CONSISTENCY_FIGURE, multi.NOISE_FIGURE):
        assert Path(path).is_file()
        assert Path(path).stat().st_size > 100
