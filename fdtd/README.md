# Finite Si mark: FDTD source data

`run_fdtd.py` launches a new hidden Ansys Lumerical FDTD v252 session and runs the native non-MPI CPU engine with eight threads. It does not attach to an existing CAD window. The finite geometry has PML on every boundary, an ordinary plane-wave source at normal incidence, and no periodic boundary conditions.

The physical 2D coordinates are x (across the grooves) and y (surface normal). The electric field is Ez, the incident wave travels toward negative y, and reflected waves travel toward positive y. The Si ridge top is y=0, the groove floor and substrate top are y=-100 nm. The period is 2 um and each ridge width is 1 um. The finite N-cell mark is centered on x=0; for every even N used here, x=0 is a groove center and ridge centers are at odd multiples of 1 um. Outside the mark, the flat Si surface stays at the groove-floor height.

The source spans the simulation width N*2 um+8 um. The extra 4 um on each side separates the measured mark aperture from the source and PML edges. This finite open-domain injection still has a measurable nonuniform incident field near the edges; therefore each geometry size has an air-reference simulation with exactly the same source, mesh, boundaries, and monitor coordinates. The source spans 500--565 nm; DFT data are sampled only at 532 nm.

The additional `N010_sourceThroughPML_diagnostic` air run increases the source span by 2 um on each side, while keeping the FDTD region fixed. The Ez arrays at I, R, and R2 are bitwise identical to the original N10 air run (`comparison.json`). Extending the requested source rectangle outside the FDTD user span therefore does not change the effective injected fields in this installation. The production geometry and source remain the original ones.

The material is the bundled `Si (Silicon) - Palik`. `material.json` records both the database value (`getindex`) and the dispersive fit actually used by FDTD over the source band (`getfdtdindex`). The RCWA comparison uses the latter so both solvers compare the same complex permittivity.

Reproduction on the installed workstation:

```powershell
& 'D:\Ansys Inc\v252\Lumerical\python\python.exe' 'D:\Si532_spectrum_native_20261031\fdtd\run_fdtd.py' --N 10 20 50 100 200 400
& 'D:\Ansys Inc\v252\Lumerical\python\python.exe' 'D:\Si532_spectrum_native_20261031\fdtd\run_fdtd.py' --N 20 --dx 10 --dy 5
```

The ASCII path is a junction to the same project directory, required by this local Lumerical installation. On another computer, set `LUMERICAL_ROOT` to the Lumerical installation directory and `SI532_NATIVE_ROOT` to the ASCII path of this project root (or an ASCII junction to it). The default paths reproduce this workstation's run. The script verifies that the selected ASCII path points to the same project before saving a model. An installed solver license is required for reproduction.

Each `cases/<settings>/<air|mark>/` directory contains the following evidence:

- `model.fsp`: the genuine Lumerical project including the solved frequency-domain monitor data.
- `model.lsf`: the equivalent geometry, source, boundary, and monitor commands.
- `input.json`: the exact numerical and physical settings.
- `material.json`: database and fitted optical constants.
- `engine_stdout.log`, `engine_stderr.log`: unchanged local solver output and convergence evidence, retained locally and excluded from the public repository.
- `engine_stdout_public.log`: the public numerical excerpt, retaining solver version, CPU threads, mesh size, progress, autoshutoff, and runtime/memory statistics. Machine, account, licensing, and absolute-path lines are omitted. If the solver creates a separate `model_p0.log`, its public excerpt uses the same `_public.log` suffix.
- `fields_R.npz`, `fields_R2.npz`: six complex E/H field components at y=1.0 um and y=1.2 um, both above the source plane at y=0.7 um.
- `air/fields_I.npz`: incident-wave reference measured at y=0 in air. No incident monitor is placed inside the mark material.
- `DONE.json`: solver return code, measured runtime, full field-array schema, mesh counts, simulated duration, autoshutoff status, and material values.

NPZ files preserve the original native `getdata` arrays, including x/y/z/f, Ex/Ey/Ez/Hx/Hy/Hz, and any available Yee-grid offsets. Aliases `x_m`, `y_m`, and `frequency_hz` make the units explicit. Frequency-domain field arrays have shape `(Nx,1,1,1)` in this single-frequency 2D experiment. The monitor uses Lumerical's `specified position` interpolation, which co-locates E/H components for subsequent wave decomposition.

The reflected field used in the comparison is `R_mark - R_air`, which removes source back leakage. It preserves the reflected zero order and all grating orders. Divide by the complex mean of `I_air` over the same N-period crop to fix incident amplitude and phase at the ridge-top plane y=0. Remove the outgoing propagation factor from the measurement plane to that reference plane. The final report and the root analysis script derive and test the required finite-aperture normalization.

`quick_extract.py` provides a small independent diagnostic extraction over x in `[-NP/2,NP/2)` and writes `quick_summary.json`. The formal analysis, plots, and uncertainty discussion are generated by the project-level scripts.

`run_checks.py` runs the explicit numerical controls: N200 at 1200 fs with otherwise identical coarse settings, N50 with a 10 nm by 5 nm mesh, and a three-level mesh comparison at N20. The coarse and fine meshes use 24 and 48 PML layers, respectively, so their physical PML thickness is fixed. The finest 5 nm by 2.5 nm mesh requests 96 layers, but the saved actual solver setting is 64 because the installed GUI applies its default maximum. `observed_settings.json` records the effective value directly; the report and manifest use 64 and retain 96 separately as the requested value. The finest control therefore does not fully hold the physical PML thickness fixed. The earlier N20 10 nm by 5 nm run with 24 PML layers is retained as a separate PML-thickness control; doubling those layers to 48 changes the extracted coefficients by at most 1.204e-5. `summarize_fdtd.py` collects runtime, convergence, requested/effective PML layers, and independent R-to-R2 propagation diagnostics into `solver_manifest.json` and `.csv`.

The ordinary N200 run with a 1200 fs maximum stops the mark simulation at 621.30 fs when its autoshutoff criterion is met; the corresponding air run completes all 1200 fs. Therefore a separate `N200_dx20_dy10_t1200_strictmark` run reduces the mark's autoshutoff threshold to 1e-12. `run_strict_time.py` performs this additional mark solve and is also called by `run_checks.py`. It reuses the existing air reference that truly completed all 1200 fs, after verifying identical source, material band, mesh, boundaries, and geometry size. `AIR_REFERENCE.json` links to that genuine reference and records the differing stopping thresholds and actual durations. No new air run, new air DONE file, or new air runtime is claimed for this strict case.

The strict mark run actually completes 1200.010 fs, without early termination; its final autoshutoff value is 4.247e-10 and its additional measured runtime is 435.46 s. Compared with the actual 600.0051 fs baseline, the maximum change among the seven complex diffraction coefficients is 1.375e-5; their relative complex L2 change is 0.003195%. This temporal change is far smaller than the measured coarse-to-fine mesh change of about 0.00516 in absolute complex amplitude. All 24 genuine FDTD engine jobs, including the air-source diagnostic and this strict mark control, return exit code zero. The complete numerical controls and actual durations are in `numerical_controls.json` and `solver_manifest.json`.

`publish_diagnostics.py` generates those public log excerpts without changing the originals. Its `audit_saved_projects` function loads solved projects and queries the actual boundary conditions, 2D dimension, source properties, mesh steps, and PML layers. The resulting `saved_project_audit.json` directly verifies the all-PML configuration; the source's UI name `Bloch/periodic` does not imply periodic simulation boundaries.
