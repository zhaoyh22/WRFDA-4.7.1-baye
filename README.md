# WRFV4.7.1 with Bayesian Radar RQV Assimilation

This repository is based on the official WRF Model `v4.7.1` release and adds a
Bayesian radar humidity (`rqv`) assimilation option to WRFDA.  The new option is
designed for weak-echo and clear/missing-echo radar situations where the original
WRFDA radar `rqv` retrieval can be too limited or can leave false background
convection insufficiently suppressed.

For the detailed algorithm notes, parameter definitions, and code map, see
[doc/README.radar_rqv_bayes.md](doc/README.radar_rqv_bayes.md).

## What changed relative to original WRFDA

The original WRFDA radar humidity option uses reflectivity thresholds to create
estimated water-vapor observations.  This version keeps that original behavior
as the fallback and adds the following features:

| New feature | Purpose | Main controls |
| --- | --- | --- |
| Weak-echo Bayesian `rqv` adjustment | Estimate pseudo relative humidity from nearby background points with similar reflectivity. | `radar_rqv_bayes_opt`, `radar_rqv_bayes_radius`, `radar_rqv_bayes_sigma`, `radar_rqv_bayes_obs_thresh` |
| Clear/missing-echo false-convection suppression | Dry the background where model reflectivity is strong but radar indicates no, weak, or missing echo. | `radar_rqv_bayes_bg_thresh`, `radar_rqv_bayes_clear_radius`, `radar_rqv_bayes_noise`, `radar_rqv_bayes_min_pts`, `radar_rqv_bayes_rh_min` |
| Internal `rqv`-only pseudo radar observations | Add pseudo observations at strong-background-reflectivity grid cells inside the radar diagnosis area but without a radar record. | `radar_rqv_bayes_pseudo_opt`, `radar_rqv_bayes_pseudo_h_stride`, `radar_rqv_bayes_pseudo_v_stride`, `radar_rqv_bayes_pseudo_max` |
| SCR diagnostic NetCDF output | Write masks showing the diagnosis area, candidate cells, and cells accepted by Bayesian `rqv` QC. | `radar_rqv_bayes_scr_diag` |
| Configurable Bayesian `rqv` observation error | Use either the original WRFDA-style error or a qsat-based configurable error. | `radar_rqv_bayes_err_opt`, `radar_rqv_bayes_err_min`, `radar_rqv_bayes_err_factor` |

After this repository is pushed to GitHub, the full code difference from the
official base can be inspected with:

```bash
git diff v4.7.1..main
```

On GitHub, use the compare page:

```text
https://github.com/<owner>/<repo>/compare/v4.7.1...main
```

## Main modified files

- `Registry/registry.var`: adds the new `&radar_da` namelist variables.
- `var/README.namelist`: documents the new namelist variables.
- `var/da/da_radar/da_radar.f90`: connects the new radar Bayesian routines.
- `var/da/da_radar/da_radar_rqv_bayes.inc`: implements the Bayesian `rqv`
  estimator.
- `var/da/da_radar/da_get_innov_vector_radar.inc`: applies Bayesian `rqv`
  innovations during radar innovation generation.
- `var/da/da_radar/da_write_radar_rqv_scr_diag.inc`: writes `scr_diag.nc`.
- `var/da/da_setup_structures/da_setup_radar_rqv_bayes_pseudo.inc`: creates
  internal `rqv` pseudo radar observations.
- `var/da/da_setup_structures/da_setup_obs_structures_radar.inc`: calls the
  pseudo-observation setup after radar observations are read.

## Quick use

Build WRFDA as usual after configuring the WRFDA option:

```bash
./clean -a
./configure wrfda
./compile all_wrfvar
```

Enable radar `rqv` assimilation and the Bayesian extension in `namelist.input`:

```namelist
&wrfvar4
 use_radarobs = .true.,
 use_radar_rqv = .true.,
 use_radar_rf  = .false.,
 use_radar_rhv = .false.,
/

&radar_da
 radar_rqv_bayes_opt = 1,
 radar_rqv_bayes_radius = 20,
 radar_rqv_bayes_clear_radius = 33,
 radar_rqv_bayes_sigma = 5.0,
 radar_rqv_bayes_bg_thresh = 30.0,
 radar_rqv_bayes_obs_thresh = 25.0,
 radar_rqv_bayes_noise = 0.0,
 radar_rqv_bayes_min_pts = 3,
 radar_rqv_bayes_rh_min = 10.0,
 radar_rqv_bayes_err_opt = 1,
 radar_rqv_bayes_pseudo_opt = 0,
 radar_rqv_bayes_scr_diag = 1,
/
```

When `radar_rqv_bayes_scr_diag = 1`, WRFDA writes `scr_diag.nc` during the first
outer iteration before minimization.  It contains `scr_valid`, `scr_candidate`,
and `scr_assimilated` masks.

## Recommended workflow for experiments

1. Run a control experiment with `radar_rqv_bayes_opt = 0`.
2. Run the Bayesian experiment with `radar_rqv_bayes_opt = 1` and
   `radar_rqv_bayes_scr_diag = 1`.
3. Inspect `scr_diag.nc` to confirm the clear/missing-echo suppression area.
4. If too many cells are selected, increase `radar_rqv_bayes_bg_thresh`, increase
   `radar_rqv_bayes_noise`, reduce `radar_rqv_bayes_clear_radius`, or add
   height bounds with `radar_rqv_h_lbound` and `radar_rqv_h_ubound`.
5. Enable `radar_rqv_bayes_pseudo_opt = 1` only after the direct radar-observation
   path is verified, because pseudo observations can increase the observation
   count substantially.