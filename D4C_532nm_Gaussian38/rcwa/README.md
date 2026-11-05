# Matched-material TE RCWA and the 2026-11-05 D4C target

This directory performs actual Fourier modal calculations for the eight structures in `../experiment.json`. It does not substitute a local binary reflection mask for Maxwell scattering, and does not use finite-FDTD output to fit the RCWA fields.

## Physical and phase definitions

- Two-dimensional TE: electric field Ez, invariant z, exp(-i omega t).
- y points upwards. Si substrate occupies y <= 0; Si relief is 0 < y < 300 nm.
- For cap000 there is no SiO2, including no groove fill. For the other cases SiO2 fills the grooves and extends to y = 300 nm + cap.
- Incident illumination travels in -y. Unit incident electric amplitude is defined at y = 0, even when that plane intersects the structure; this is the incident plane wave extrapolated in air.
- `r_total_y0` is the reflected propagating electric coefficient extrapolated to y = 0 using the same incident reference. `r_common_y085` propagates its outgoing field to y = .85 um, retaining incident amplitude at y0 as denominator.
- The transverse origin is the center of a symmetric mother cell. No half-cell phase alignment is applied after solving.
- Material n+ik comes from `../fdtd/materials.json`, specifically the actual FDTD dispersive-fit value at 532 nm with a 500–600 nm fitting band. Database interpolation is not used as a different RCWA material.

`te_planarized_rcwa.py` expands the exact rectangular Si/air or Si/SiO2 permittivity coefficients and solves the coupled TE layer modes. The uniform cap is another modal layer joined by input admittances. Outgoing square roots satisfy Im(ky)>=0. Both propagating and evanescent modes participate in layer matching; only the exported far-field `r_total_y0` is zero outside the air light cone.

`provenance/original_si_rcwa.py` is a byte copy of the earlier independently validated single patterned-layer implementation, retained for a direct complex regression. It is not modified. `solver_validation.json` records 32 checks, including uniform thin-film Fresnel amplitudes at several incidence angles and cap thicknesses, lossless-film flux balance, original solver regression, mirror symmetry and passivity.

The semi-infinite absorbing Si substrate has no transmitted far-depth beam. `T_interface` is inward power at y=0 and is ultimately absorbed in that substrate; `finite_layer_absorption` is absorption in the finite relief+cap layers. Thus R + T_interface + finite_layer_absorption = 1. This algebraic sum alone is not an independent accuracy proof; harmonic comparisons and analytic checks are also supplied.

## Numerical outputs

Run `python validate_rcwa.py` and `python run_reference.py` using NumPy 2 and SciPy. The latter caches each actual calculation in `cases/<variant>/N<order>_q0.npz` with solver/config hashes and wall time. Max order N means 2N+1 retained harmonics.

At least N160, N240 and N320 are computed. If the N240->N320 difference is not below relative complex L2 0.0005 for collected m=-18..18 and 0.001 for its nonzero subset, N480/N640 controls follow. These numerical cutoffs are a stated working precision, not proof of an arbitrarily exact solution.

- `all_orders.csv`: signed diffraction orders, complex electric amplitudes, |r|^2 and power efficiencies.
- `convergence.csv`: all truncations and changes on the common collected orders.
- `structure_summary.csv`, `reference_summary.json`: I0, zero terms, selected truncation and provenance.
- Per structure `reference.npz`: selected full complex coefficients and flux outputs.
- Per structure `d4c_targets.csv/.npz`: complete primary `S_SUM`, `S_DIFF_same_kernel` and `AC` curves, plus symbolic K2/K4 sensitivity coefficients. The user confirmed T_m=|r_m|² and authorized the interpolation during this task. Older provisional-named arrays remain identical compatibility aliases.

## Supplied D4C formula, including unresolved definitions

The authority is `../inputs/D4C_FDTD_user_formula_20261105.md`. This calculation uses the NEW formula, not the previous K0_port and outside-1/2 formula:

\[
S_\chi(s)=2\sum_{m=1}^{18}K_mT_m[1-\cos(2mGs+\phi_m)]
+[T_0+L_\chi(1-T_0)]I_0.
\]

The user confirmed T_m=|r_(+m)|^2 after supplying the original note. We use I0=|r0|^2 and phi_m=arg(r_(+m) conj(r_(-m))). Both tested geometries are mirror symmetric, so the computed pair phases are numerically zero. More general asymmetric marks require separate derivation, not this symmetric intensity shortcut. The original supplied note is preserved without silently rewriting its then-unresolved definitions.

L_SUM=.00138, L_DIFF=.000431 and

\[
T_0=\operatorname{erfc}\left(\frac{5\times10^{-4}\sqrt2}{.778\times3.76\times10^{-4}}\right)
=0.000629728295091705.
\]

The resulting zero coefficients per I0 are .002008859270044478 (SUM) and .0010604568821965205 (DIFF). The DIFF data retain the SAME supplied 1-cos kernel only as a leakage-parameter comparison. No claim is made that the physical DIFF interferometer has this sign; the original note leaves it unresolved.

The original note did not list K2 and K4, and the user subsequently authorized linear interpolation. The complete primary target uses K2=.4947 and K4=.49135. For traceable sensitivity, the decomposition remains

\[
S_\chi=S_{\chi,\rm known}+K_2C_2(s)+K_4C_4(s),\quad
C_m(s)=2|r_m|^2[1-\cos(2mGs+\phi_m)].
\]

Columns `coefficient_K2` and `coefficient_K4` allow future changes without rerunning RCWA; they are a sensitivity appendix, not a missing part of the primary target. The original new K13=.4478 and K14=.4380 are used, with the corrected/prior-reading/interpolated provenance preserved in `K_table.csv`. Zero-order conclusions do not depend on K2/K4. Legacy provisional-named NPZ arrays are retained only as compatibility aliases of the now-authorized complete target.

## Gaussian and zero-stop controls

`gaussian_angle_controls.py` computes small-q responses at N160, sharing the same patterned-layer eigensystem across the three uniform-cap thicknesses for each geometry and fill state. Its ±q reflection symmetry is exact for these centered symmetric structures. The 25k rad/m input-q sampling and ±200k range are diagnostics of the frozen-q0 approximation; they are not a convergence-certified finite-mark or complete Gaussian library. The full coefficient arrays remain available for inspection.

For an ideal isolated Gaussian zero-order field envelope exp[-(x/w0)^2], the power outside an opaque symmetric angular stop is erfc(h*w0/sqrt(2)). With w0=19 um, the supplied T0 corresponds to h=.17991979895041987 rad/um. The old h=.2 rad/um corresponds instead to .00014469608785023995. The angular diagnostics stop at ±.2 rad/um and omit that latter input tail fraction; it is small compared with total power but about 23% of the desired zero leakage. They must not be used directly to predict precise blocked-zero leakage without extending angular support. A frozen-q0 Gaussian zero template can retain the full tail analytically.

An angular stop acts on every field occupying its angular interval. A modal removal of the entire m0 component is a different operation. Neither operation alone implements the supplied additive blocked-light leakage term. In particular, assigning the stop an internal power transmission tau and separately adding T0*P0 can double count transmitted zero light when they describe the same channel.
