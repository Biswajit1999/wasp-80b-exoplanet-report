# Data sources

## TESS light curve

- File: `tess2022190063128-s0054-0000000243921117-0227-s_lc.fits`
- Archive: Mikulski Archive for Space Telescopes (MAST), TESS SPOC light-curve product
- TESS sector: 54
- TIC target ID: 243921117
- MAST observation ID: 91555059
- MAST data URI: `mast:TESS/product/tess2022190063128-s0054-0000000243921117-0227-s_lc.fits`
- Exact download URL: <https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS%2Fproduct%2Ftess2022190063128-s0054-0000000243921117-0227-s_lc.fits>
- Collection DOI: [10.17909/t9-nmc8-f686](https://doi.org/10.17909/t9-nmc8-f686) (TESS 2-minute light curves, all sectors; sector 54 used here)
- Retrieved: 2026-08-15
- SHA-256: `cc6a01fbc5090a11aa09d7c9c785152405dea520d6008460242ad2c2448964c8`

The FITS file is stored unmodified. The analysis reads `TIME`, `PDCSAP_FLUX`,
`PDCSAP_FLUX_ERR`, and `QUALITY`. PDCSAP flux is the SPOC light curve with common
instrumental trends removed and aperture/crowding corrections applied; this does
not make it free of residual stellar or instrumental systematics.

## System parameters

- File: `system_parameters.csv`
- Service: NASA Exoplanet Archive TAP, `pscomppars` table
- Exact query: <https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name%2Chostname%2Cra%2Cdec%2Cpl_orbper%2Cpl_tranmid%2Cpl_trandur%2Cpl_rade%2Cpl_bmasse%2Cpl_eqt%2Cpl_orbsmax%2Csy_dist%2Csy_tmag%2Cst_teff%2Cst_rad%2Cst_mass%2Cdisc_year%2Cdiscoverymethod%2Cdisc_refname%2Cdisc_pubdate%2Cdisc_facility+from+pscomppars+where+pl_name%3D%27WASP-80+b%27&format=csv>
- Retrieved: 2026-08-15

The saved row is the input actually used by `scripts/analyze_transit.py`; the
analysis does not query a changing live service at run time.

## JWST/NIRCam transmission and emission spectra

- Committed derivative: `spectra/wasp80b_nircam_source_data_rebinned.csv`
- Primary publication: [Bell et al. (2023), Nature](https://doi.org/10.1038/s41586-023-06687-0)
- Original workbook: [Source Data Fig. 3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-023-06687-0/MediaObjects/41586_2023_6687_MOESM3_ESM.xlsx)
- Retrieved: 2026-08-16
- Original XLSX SHA-256: `d595847ff4bd6e4cc7878bef1e666cfea190ef780e0846135809ab5a9930a2a2`

The workbook contains the published fiducial transmission and emission spectra
and posterior contours for the full, no-CH4, and no-H2O models. To keep the
repository compact, five adjacent data channels were inverse-variance combined;
the model contours were linearly interpolated onto the original channel centers
and combined with the matching weights. The final partial group is retained. All
depths and uncertainties in the committed CSV are dimensionless. Model values in
the workbook were converted from percent to dimensionless depth before binning.
This deterministic compression is intended for robustness diagnostics, not for
reproducing the paper's retrieval or its reported detection significances.

`published_methane_abundances.csv` transcribes the two CH4 abundance summaries
and reported detection significances from Bell et al. It is used only for an
explicit consistency calculation; it is not a new abundance retrieval.


## Additional TESS sectors for robustness analysis

All are unmodified standard-cadence SPOC light curves from the same [MAST TESS collection](https://doi.org/10.17909/t9-nmc8-f686).

- Sector 54: `tess2022190063128-s0054-0000000243921117-0227-s_lc.fits` (1,918,080 bytes)
  - MAST URI: `mast:TESS/product/tess2022190063128-s0054-0000000243921117-0227-s_lc.fits`
  - SHA-256: `cc6a01fbc5090a11aa09d7c9c785152405dea520d6008460242ad2c2448964c8`
