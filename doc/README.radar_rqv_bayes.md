# Bayesian Radar RQV Assimilation in WRFDA

This document describes the Bayesian radar humidity (`rqv`) assimilation changes
added on top of WRF `v4.7.1`.

## Goal

The original WRFDA radar `rqv` option estimates water vapor from radar
reflectivity using fixed reflectivity and relative-humidity thresholds.  That
logic is useful for saturated echoes, but it is less flexible in two common
cases:

- weak observed echo, where a fixed saturated/unsaturated threshold can be too
  crude;
- no, weak, or missing observed echo with strong background reflectivity, where
  false model convection should be dried more explicitly.

The new option keeps the original WRFDA `rqv` retrieval as a fallback and adds
a local Bayesian estimator that derives pseudo relative humidity from nearby
background reflectivity and RH samples.

## New namelist variables

The variables are defined in `Registry/registry.var` under `namelist,radar_da`
and documented in `var/README.namelist`.

| Variable | Default | Meaning |
| --- | ---: | --- |
| `radar_rqv_bayes_opt` | `0` | `0`: off; `1`: enable Bayesian `rqv` adjustment. |
| `radar_rqv_bayes_radius` | `20` | Horizontal radius, in grid points, for weak-echo Bayesian `rqv`. |
| `radar_rqv_bayes_clear_radius` | `33` | Horizontal radius, in grid points, for clear/missing-echo suppression. |
| `radar_rqv_bayes_sigma` | `5.0` | Reflectivity Gaussian weighting scale, in dBZ, for weak echoes. |
| `radar_rqv_bayes_bg_thresh` | `30.0` | Background reflectivity threshold, in dBZ, used to identify false convection candidates. |
| `radar_rqv_bayes_obs_thresh` | `25.0` | Observed reflectivity upper threshold, in dBZ, for weak-echo adjustment. |
| `radar_rqv_bayes_noise` | `0.0` | Reflectivity value, in dBZ, treated as no or weak echo. |
| `radar_rqv_bayes_min_pts` | `3` | Minimum number of environmental samples required for clear/missing-echo suppression. |
| `radar_rqv_bayes_rh_min` | `10.0` | Minimum pseudo RH, in percent, for clear/missing-echo suppression. |
| `radar_rqv_bayes_err_opt` | `1` | `1`: original-style error, `max(0.001, 0.20*rqvo)`; `2`: qsat-based configurable error. |
| `radar_rqv_bayes_err_min` | `0.0005` | Minimum Bayesian `rqv` error, in kg/kg, for `err_opt=2`. |
| `radar_rqv_bayes_err_factor` | `0.10` | Error factor multiplied by model saturation mixing ratio for `err_opt=2`. |
| `radar_rqv_bayes_pseudo_opt` | `0` | `0`: off; `1`: append internal `rqv`-only pseudo radar observations. |
| `radar_rqv_bayes_pseudo_h_stride` | `1` | Horizontal stride, in grid points, for internal pseudo observations. |
| `radar_rqv_bayes_pseudo_v_stride` | `1` | Vertical stride, in model levels, for internal pseudo observations. |
| `radar_rqv_bayes_pseudo_max` | `0` | Maximum internal pseudo observations per MPI task; `0` means no limit. |
| `radar_rqv_bayes_scr_diag` | `0` | `0`: off; `1`: write SCR diagnostic NetCDF file. |

## Algorithm

The implementation lives mainly in
`var/da/da_radar/da_get_innov_vector_radar.inc` and
`var/da/da_radar/da_radar_rqv_bayes.inc`.

### Preprocessing

When `use_radar_rqv = .true.` and `radar_rqv_bayes_opt > 0`, WRFDA first:

1. computes background reflectivity from background hydrometeors at model grid
   points;
2. collects background reflectivity and background RH into global arrays with
   `da_patch_to_global`;
3. broadcasts those arrays for MPI-parallel use;
4. marks model grid cells where valid observed reflectivity echo exists, so the
   clear/missing-echo suppression path does not dry cells that already have
   direct echo evidence.

### Weak-echo adjustment, `bayes_type = 1`

This path is used when:

- observed reflectivity is valid;
- observed reflectivity is greater than `radar_rqv_bayes_noise`;
- observed reflectivity is smaller than `radar_rqv_bayes_obs_thresh`;
- background reflectivity is greater than `radar_rqv_bayes_noise`;
- the observation height passes the normal `rqv` height checks.

For nearby grid points within `radar_rqv_bayes_radius`, the scheme weights each
background RH sample by:

```text
weight = exp(-((bg_rf - obs_rf)^2) / (2 * radar_rqv_bayes_sigma^2))
```

The pseudo RH is the weighted mean.  The adjustment is accepted only if the
reflectivity innovation and RH adjustment have consistent signs:

```text
(obs_rf - bg_rf) * (pseudo_rh - bg_rh) > 0
```

### Clear/missing-echo suppression, `bayes_type = 2`

This path is used when:

- observed reflectivity is missing, no echo, weak echo, or the radar file marks
  a non-precipitating observation;
- background reflectivity is at least `radar_rqv_bayes_bg_thresh`;
- the observation height passes the normal `rqv` height checks;
- the grid cell is not already marked as containing a valid observed echo.

For nearby grid points within `radar_rqv_bayes_clear_radius`, the scheme uses
environmental samples satisfying:

```text
bg_rf <= radar_rqv_bayes_noise
bg_rh < target_bg_rh
```

If at least `radar_rqv_bayes_min_pts` samples are found, the pseudo RH is their
mean, clipped to at least `radar_rqv_bayes_rh_min`.  The adjustment is accepted
only if:

```text
pseudo_rh < target_bg_rh
```

The resulting `rqv` innovation is clipped to be non-positive in this
clear/missing-echo path, so it acts as a drying constraint.

## Internal pseudo observations

If `radar_rqv_bayes_pseudo_opt > 0`, WRFDA appends internal pseudo radar
profiles after radar observations are read.  These pseudo observations are
created by `var/da/da_setup_structures/da_setup_radar_rqv_bayes_pseudo.inc`.

A grid point can become an internal pseudo observation when:

- it is inside the expanded radar diagnosis area;
- it has no existing radar record at that grid level;
- it passes optional `radar_rqv_h_lbound` and `radar_rqv_h_ubound` height bounds;
- background reflectivity is at least `radar_rqv_bayes_bg_thresh`.

Pseudo observations are labeled as:

```text
platform = RQV_BAYES
name     = RQV_BAYES_PSEUDO
id       = RQV_BAYES_PSEUDO
```

For FGAT runs, pseudo observations are currently appended after all radar files
are read and are assigned to the final FGAT time slot.

## SCR diagnostic output

Set:

```namelist
radar_rqv_bayes_scr_diag = 1
```

WRFDA writes diagnostics during outer iteration 1 before minimization, on the
root process after gathering the masks. With `num_fgat_time = 1`, the legacy
filename `scr_diag.nc` is retained. With `num_fgat_time > 1`, each diagnosed slot
writes `scr_diag_t<time_index>.nc` (for example, `scr_diag_t1.nc`,
`scr_diag_t2.nc`, and `scr_diag_t3.nc`). Slot 1 also has a suffix; no unsuffixed
file is produced by a multi-slot run.

`time_index` is the one-based FGAT slot (`iv%time`), not outer iteration or a
calendar timestamp. The driver in `da_get_innov_vector.inc` visits slots from
`num_fgat_time` down to 1 and selects observations using `plocal(time_index-1)`
and `plocal(time_index)`. Every file retains the text global attributes
`time_index` and `analysis_date`; `analysis_date` is the analysis date, not the
individual slot's valid time. Slots without usable non-pseudo radar observations
can be skipped. Internal pseudo observations still belong to the final slot.

Each file contains three byte masks:

| Variable | Meaning |
| --- | --- |
| `scr_valid` | Grid cells inside the SCR diagnosis area. |
| `scr_candidate` | Background reflectivity exceeds `radar_rqv_bayes_bg_thresh` and there is no/weak/missing observed echo. |
| `scr_assimilated` | Candidate cells accepted by Bayesian `rqv` QC and assimilated. |

## Minimal namelist example

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

## Suggested tuning

- Increase `radar_rqv_bayes_bg_thresh` if suppression is too broad.
- Increase `radar_rqv_bayes_clear_radius` if too few environmental samples are
  found in clear/missing-echo areas.
- Increase `radar_rqv_bayes_min_pts` to make clear/missing-echo suppression more
  conservative.
- Increase `radar_rqv_bayes_sigma` to allow weaker reflectivity similarity in
  the weak-echo estimator; decrease it to make the estimator more selective.
- Use `radar_rqv_h_lbound` and `radar_rqv_h_ubound` to restrict the scheme to a
  vertical layer.
- Keep `radar_rqv_bayes_pseudo_opt = 0` for first tests.  Turn it on only after
  the direct radar observation path and the SCR diagnostic masks look reasonable.

## Build notes

The WRFDA dependency list was updated so the new `.inc` files are dependencies
of `da_radar.o` and `da_setup_structures.o`.

`var/build/da.make` also adds `da_radar.o` to the same special compile rule used
for `da_radiance1.o`, `module_radiance.o`, and `da_setup_structures.o`.

`var/obsproc/src/Makefile` adds a no-optimization compile rule for
`module_err_afwa.o`.  This is a build-stability change and is not part of the
Bayesian radar algorithm.

## SCR regression checks

Run `python3 -m unittest discover -s var/test/radar_rqv_bayes -v` with
GNU Fortran, NetCDF Fortran (`nf-config`), NumPy and SciPy installed. The test
compiles the production NetCDF writer in a minimal host module and reads the
resulting files. It checks single-slot compatibility, reverse-order multi-slot
writes, distinct contents for all three masks, the time-index attribute, and
indices above fixed-width naming limits. It also checks the source wiring of
the slot loop and module import. CI runs this isolated regression; it does not
build or run the full WRFDA model.

For operational validation, run a real FGAT case with radar data in at least two
slots and SCR diagnostics enabled. Check the logs for each expected output,
inspect each file's `time_index`, dimensions and three masks, and compare masks
with the corresponding background and observations. Repeat in MPI mode to
validate collective gathering and root-only output. This requires a configured
WRFDA build and case data and is separate from the isolated writer test.
Repeated runs in the same directory (or another write of the same slot) still
replace the same filename via `nf_clobber`; use separate run directories.
