"""D4C/B-method angular-spectrum candidates, with no electromagnetic solver.

The supplied 2026-11-05 formula is the contract. SUM/DIFF here identify the
two supplied leakage constants; they do not assert complementary hardware
ports. All functions use the entire image plane through angular quadrature.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence
import hashlib
import json
import math

import numpy as np
from scipy.optimize import least_squares

Z0 = 376.730313668
T0 = math.erfc(5e-4 * math.sqrt(2) / (0.778 * 3.76e-4))
LEAKAGE = {"SUM": 0.00138, "DIFF": 0.000431}
K_VALUES = {1: .4957, 2: .4947, 3: .4937, 4: .49135, 5: .4890, 6: .4860, 7: .4875,
            8: .4783, 9: .4704, 10: .4708, 11: .4629, 12: .4540,
            13: .4478, 14: .4380, 15: .4293, 16: .4200, 17: .4056,
            18: .2480}
K_CONFIRMED = (3, 6, 7, 9, 10, 11, 12, 14, 15)
K_INTERPOLATED = (2, 4)
K_PROVISIONAL = tuple(sorted(set(K_VALUES) - set(K_CONFIRMED) - set(K_INTERPOLATED)))
KNOWN_ORDERS = tuple(sorted(K_VALUES))


def zero_coefficient(channel="SUM", transmission=T0, leakage=None):
    """Dimensionless scalar c_chi = T0 + L_chi (1-T0)."""
    channel = channel.upper()
    leak = LEAKAGE[channel] if leakage is None else float(leakage)
    if not 0 <= transmission <= 1 or not 0 <= leak <= 1:
        raise ValueError("Power transmission and leakage must lie in [0, 1]")
    return float(transmission + leak * (1 - transmission))


def nearest_order(kx_rad_m, period_um=16.):
    """Signed nearest diffraction band, with half-integer ties away from 0."""
    q = np.asarray(kx_rad_m, dtype=float) * period_um * 1e-6 / (2*np.pi)
    return (np.sign(q) * np.floor(np.abs(q) + .5 + 1e-12)).astype(int)


def quadrature(kx_rad_m):
    k = np.asarray(kx_rad_m, dtype=float)
    if k.ndim != 1 or len(k) < 3 or not np.all(np.isfinite(k)) or np.any(np.diff(k) <= 0):
        raise ValueError("k must be finite, one-dimensional and increasing")
    w = np.empty_like(k)
    w[1:-1] = (k[2:] - k[:-2])/2
    w[0], w[-1] = (k[1]-k[0])/2, (k[-1]-k[-2])/2
    return w


def angular_weights(kx_rad_m, Pin_W_m, wavelength_nm=532., NA=.6):
    """TE E-spectrum weights, including dk and 1/Pin exactly once.

    A(k)=integral E(x) exp(-ikx) dx. A has V units for raw E in V/m.
    If A is divided by complex Ein, Pin must be divided by |Ein|^2.
    This gamma is an electromagnetic flux factor, never a substitute for K_m.
    """
    if not np.isfinite(Pin_W_m) or Pin_W_m <= 0 or not 0 < NA <= 1:
        raise ValueError("Need actual positive incident power and 0 < NA <= 1")
    k = np.asarray(kx_rad_m, dtype=float)
    k0 = 2*np.pi/(wavelength_nm*1e-9)
    gamma = np.sqrt(np.maximum(0, 1-(k/k0)**2))
    return quadrature(k)*gamma/(4*np.pi*Z0*Pin_W_m)*(np.abs(k) <= NA*k0)


def stop_field(kx_rad_m, halfwidth_rad_um=0., center_rad_um=0., power_transmission=0.):
    """Actual angular mask; sqrt(power transmission) acts on the field.

    With halfwidth=0 no stop is applied, including at the single k=0 sample.
    A finite angular stop acts on every overlapping physical order component.
    """
    if halfwidth_rad_um < 0 or not 0 <= power_transmission <= 1:
        raise ValueError("Invalid stop width or power transmission")
    k = np.asarray(kx_rad_m, dtype=float)
    if halfwidth_rad_um == 0:
        return np.ones_like(k)
    inside = np.abs(k*1e-6-center_rad_um) <= halfwidth_rad_um
    return np.where(inside, math.sqrt(power_transmission), 1.)


def stop_cell_field(kx_rad_m, halfwidth_rad_um=0., power_transmission=0.):
    """Quadrature representation of a centered hard stop, not a soft optic.

    Each trapezoidal cell receives the exact overlap fraction with [-h,h].
    The field multiplier sqrt(cell-averaged power transmission) is valid for
    centered symmetric +/-k pairing, where dc and cross have the SAME weight.
    Density within a boundary cell is approximated by its node value. This
    removes the O(dk) binary-edge error but retains boundary-density error.
    """
    k=np.asarray(kx_rad_m,float)
    if halfwidth_rad_um < 0 or not 0 <= power_transmission <= 1:
        raise ValueError("Invalid stop width or power transmission")
    if halfwidth_rad_um==0:return np.ones_like(k)
    bounds=np.r_[k[0],(k[:-1]+k[1:])/2,k[-1]]
    half=halfwidth_rad_um*1e6
    overlap=np.maximum(0.,np.minimum(bounds[1:],half)-np.maximum(bounds[:-1],-half))
    fractions=np.clip(overlap/np.diff(bounds),0.,1.)
    return np.sqrt(1-(1-power_transmission)*fractions)


def instrument_weights(orders, period_um=16., wavelength_nm=532.,
                       convention="none", K=None):
    """Return bandwise instrument weights, requiring an explicit T_m mapping.

    none: pure electromagnetic power, no empirical K.
    K: T_m is interpreted as one-sided power efficiency.
    K_over_gamma: T_m is |r_m|^2 for a scalar E-field coefficient, so the
      continuous electromagnetic gamma(k) remains in quadrature and a separate
      K_m/gamma_m compensates its order-center conversion. This is NOT a fitted
      angular equalizer and finite-width gamma(k)/gamma_m is retained.
    K_2/K_4 use the user's subsequently authorized linear interpolation, recorded
    separately from directly read/calibrated coefficients; never fitted to marks.
    """
    orders = np.asarray(orders, dtype=int)
    if convention == "none":
        return np.ones(len(orders))
    if convention not in ("K", "K_over_gamma"):
        raise ValueError("convention must be none, K or K_over_gamma")
    table = K_VALUES if K is None else K
    missing = [int(m) for m in orders if int(m) not in table]
    if missing:
        raise ValueError(f"No declared K for orders {missing}; explicitly select a common subset")
    weights = np.array([table[int(m)] for m in orders], dtype=float)
    if convention == "K_over_gamma":
        gamma_m = np.sqrt(1-(orders*wavelength_nm*1e-3/period_um)**2)
        if np.any(~np.isfinite(gamma_m)) or np.any(gamma_m <= 0):
            raise ValueError("Selected order is not propagating")
        weights /= gamma_m
    return weights


@dataclass
class BFeatures:
    """All powers normalized to the supplied actual Pin, never peak-normalized."""
    orders: np.ndarray
    dc: np.ndarray
    cross: np.ndarray
    zero_pair_dark: float
    zero_band_power: float
    zero_band_after_stop_power: float
    zero_raw_power: float
    zero_pass_power: float
    zero_blocked_power: float
    nonzero_raw_power: float
    nonzero_pass_power: float
    zero_nonzero_cross_raw_power: float
    zero_nonzero_cross_pass_power: float
    collected_power: float
    k0_sample_power: float
    period_um: float
    wavelength_nm: float
    metadata: dict


def b_features(kx_rad_m, A_ref_V, Pin_W_m, *, wavelength_nm=532.,
               period_um=16., NA=.6, selected_orders=None,
               stop_halfwidth_rad_um=0., stop_center_rad_um=0.,
               stop_power_transmission=0., A_zero_ref_V=None,
               stop_boundary_cells=True):
    """Integrate entire positive/negative angular-band pairs once.

    Returns dc_m=P(+m)+P(-m), cross_m=integral A(+k) A*(-k) w dk.
    scan_m=dc_m-2 Re[cross_m exp(-2 i m G s)]. The zero quadrature sample
    is retained in powers; its subtraction-port interference is exactly zero.
    Reflection reference-plane propagation must be applied before this function.
    If A_zero_ref_V is supplied, zero powers integrate this physical-mode
    estimate over the FULL pupil and nonzero interference uses A_total-A_zero.
    Otherwise zero powers use the nearest-order B0 angular band. These are
    distinct definitions and are recorded. The optional decomposition must be
    inferred from fields/source templates, never from the target scan curve.
    """
    k = np.asarray(kx_rad_m, float)
    a = np.asarray(A_ref_V, complex)
    if a.shape != k.shape or not np.all(np.isfinite(a)):
        raise ValueError("A must be finite and match k")
    w = angular_weights(k, Pin_W_m, wavelength_nm, NA)
    step = float(np.median(np.diff(k)))
    if np.max(np.abs(k+k[::-1])) > abs(step)*1e-6:
        raise ValueError("Need a symmetric k grid for unambiguous +/-k pairing")
    if np.max(np.abs(np.diff(k)-step)) > abs(step)*1e-6:
        raise ValueError("Need uniform k grid; resample before mirror pairing")
    band = nearest_order(k, period_um)
    available = np.arange(1, int(np.max(np.abs(band[w>0])))+1)
    orders = available if selected_orders is None else np.asarray(selected_orders, int)
    if len(set(orders)) != len(orders) or np.any(orders <= 0):
        raise ValueError("Select distinct positive orders")
    centered_cells=bool(stop_boundary_cells and abs(stop_center_rad_um)<1e-14)
    mask = (stop_cell_field(k,stop_halfwidth_rad_um,stop_power_transmission) if centered_cells
        else stop_field(k, stop_halfwidth_rad_um, stop_center_rad_um, stop_power_transmission))
    zero = band == 0
    if A_zero_ref_V is None:
        a_zero = a*zero
        a_nonzero = a
        zero_definition = "nearest-order B0 angular band, not an isolated physical m=0 mode"
    else:
        a_zero = np.asarray(A_zero_ref_V, complex)
        if a_zero.shape != a.shape or not np.all(np.isfinite(a_zero)):
            raise ValueError("A_zero_ref_V must be finite and match A_total")
        a_nonzero = a-a_zero
        zero_definition = "supplied physical-mode field estimate integrated over the full NA pupil"
    b = a_nonzero*mask
    mirror = b[::-1]
    pos = k > abs(step)*1e-8
    pair_density = w*(np.abs(b)**2+np.abs(mirror)**2)
    cross_density = w*b*np.conj(mirror)
    dc = np.array([np.sum(pair_density[pos & (band == m)]) for m in orders])
    cross = np.array([np.sum(cross_density[pos & (band == m)]) for m in orders])
    # Evaluate by a norm to avoid cancellation when a symmetric zero band is dark.
    zero_after = a_zero*mask
    zero_dark = float(.5*np.sum(w*np.abs(zero_after-zero_after[::-1])**2))
    zero_raw = float(np.sum(w*np.abs(a_zero)**2))
    zero_pass = float(np.sum(w*np.abs(zero_after)**2))
    # Only an explicit physical-mode decomposition is an additive field split.
    # Its component POWERS do not add unless this cross term vanishes.
    an = a-a_zero
    atzero = np.abs(k) <= abs(step)*1e-8
    return BFeatures(orders=orders, dc=dc, cross=cross,
        zero_pair_dark=zero_dark,
        zero_band_power=float(np.sum(w[zero]*np.abs(a[zero])**2)),
        zero_band_after_stop_power=float(np.sum(w[zero]*np.abs((a*mask)[zero])**2)),
        zero_raw_power=zero_raw, zero_pass_power=zero_pass,
        zero_blocked_power=max(0.,zero_raw-zero_pass),
        nonzero_raw_power=float(np.sum(w*np.abs(an)**2)),
        nonzero_pass_power=float(np.sum(w*np.abs(an*mask)**2)),
        zero_nonzero_cross_raw_power=float(2*np.real(np.sum(w*a_zero*np.conj(an)))),
        zero_nonzero_cross_pass_power=float(2*np.real(np.sum(w*zero_after*np.conj(an*mask)))),
        collected_power=float(np.sum(w*np.abs(a*mask)**2)),
        k0_sample_power=float(np.sum(w[atzero]*np.abs(a[atzero])**2)),
        period_um=float(period_um), wavelength_nm=float(wavelength_nm),
        metadata=dict(Pin_W_m=float(Pin_W_m), NA=float(NA),
            selected_orders=orders.tolist(),
            stop_halfwidth_rad_um=float(stop_halfwidth_rad_um),
            stop_center_rad_um=float(stop_center_rad_um),
            stop_power_transmission=float(stop_power_transmission),
            stop_quadrature="Centered hard-stop overlap of trapezoid cells; piecewise-constant node density, not physical soft edges" if centered_cells else "Point-sampled hard stop; edge discretization must be checked",
            zero_definition=zero_definition,
            nonzero_definition="A_total minus supplied A_zero before stop" if A_zero_ref_V is not None else "m>0 nearest bands of total angular spectrum",
            pair_prefactor="positive-k norm |A_s(+k)-A_s(-k)|^2; no hardware splitter prefactor",
            k0_sample_retained=True, image_plane="infinite; angular quadrature only"))


def zero_accounting(features, channel="SUM", transmission=T0, leakage=None):
    """Return transparent scalar zero terms without changing nonzero harmonics."""
    leak = LEAKAGE[channel.upper()] if leakage is None else float(leakage)
    coeff = zero_coefficient(channel,transmission,leak)
    raw,passed,blocked=features.zero_raw_power,features.zero_pass_power,features.zero_blocked_power
    return dict(raw_power=raw,passed_power=passed,blocked_power=blocked,
        measured_effective_transmission=passed/raw if raw else None,
        supplied_effective_transmission=float(transmission),leakage=leak,
        physical_zero_once=passed+leak*blocked,
        scalar_zero_once=coeff*raw,double_counted_zero=passed+coeff*raw,
        zero_definition=features.metadata["zero_definition"],
        ignored_zero_nonzero_cross_power=features.zero_nonzero_cross_pass_power,
        explanation="Scalar zero addition follows the supplied D4C readout contract; a coherent modal decomposition has a separate cross term")


def scan_from_features(features, s_um, *, convention="none", K=None, include_zero_pair=False):
    s = np.asarray(s_um, float)
    phase = np.exp(-4j*np.pi*np.outer(features.orders, s)/features.period_um)
    weights = instrument_weights(features.orders, features.period_um, features.wavelength_nm, convention, K)
    y = np.sum(weights[:,None]*(features.dc[:,None]-2*np.real(features.cross[:,None]*phase)), axis=0)
    if include_zero_pair:
        y += features.zero_pair_dark
    return y


def candidate_scan(features, s_um, *, model="scalar_zero_once", channel="SUM",
                   convention="none", K=None, gain=1., zero_transmission=T0,
                   leakage=None, zero_power_override=None):
    """Explicit alternative zero-accounting models; no mark-specific offset.

    All three principal models use ONLY m>0 interference, as the supplied note:
    physical_zero_once: P0_pass + L_chi*(P0_raw-P0_pass).
    scalar_zero_once: [T0+L_chi*(1-T0)]*P0_raw.
    double_counted_zero: P0_pass + [T0+L_chi*(1-T0)]*P0_raw; deliberately
      exposes repeated accounting if the scalar T0 already describes that stop.

    physical_stop_only: full stopped B dark-pair field, no scalar zero addition.
    stop_plus_d4c: full stopped B dark-pair field plus supplied D4C zero term;
      overlap/double accounting must be audited, not assumed physically correct.
    band_remove_plus_d4c: omit the B0 coherent pair and add supplied scalar zero
      once. Use a no-stop feature object for pure modal/band-operator comparison,
      or a stopped one to separately test loss of nonzero tails.
    """
    if model not in ("physical_zero_once", "scalar_zero_once", "double_counted_zero",
                     "physical_stop_only", "stop_plus_d4c", "band_remove_plus_d4c"):
        raise ValueError("Unknown zero accounting model")
    y = scan_from_features(features, s_um, convention=convention, K=K,
                           include_zero_pair=(model in ("physical_stop_only", "stop_plus_d4c")))
    pzero = features.zero_raw_power if zero_power_override is None else float(zero_power_override)
    if pzero < 0:
        raise ValueError("Zero power must be nonnegative")
    leak = LEAKAGE[channel.upper()] if leakage is None else float(leakage)
    if model == "physical_zero_once":
        if zero_power_override is not None:
            raise ValueError("physical_zero_once needs measured raw/pass powers, not an override")
        y += features.zero_pass_power + leak*features.zero_blocked_power
    elif model == "double_counted_zero":
        y += features.zero_pass_power + zero_coefficient(channel, zero_transmission, leakage)*pzero
    elif model != "physical_stop_only":
        y += zero_coefficient(channel, zero_transmission, leakage)*pzero
    return float(gain)*y


def d4c_target(s_um, orders, T_m, I0, *, phi_m=None, period_um=16.,
                channel="SUM", K=None, transmission=T0, leakage=None,
                T_m_definition):
    """Supplied D4C formula; T_m definition must be explicitly supplied.

    The user subsequently confirmed T_m=|r_m|^2. The caller still records the
    mapping explicitly; this function does not infer it from a complex array.
    """
    if not isinstance(T_m_definition, str) or not T_m_definition.strip():
        raise ValueError("An explicit T_m definition is mandatory")
    orders = np.asarray(orders, int)
    tm = np.asarray(T_m, float)
    phi = np.zeros(len(orders)) if phi_m is None else np.asarray(phi_m, float)
    if tm.shape != orders.shape or phi.shape != orders.shape or np.any(tm < 0):
        raise ValueError("T_m/phi must match orders and T_m must be nonnegative")
    kval = instrument_weights(orders, convention="K", K=K)
    angle = 4*np.pi*np.outer(orders, np.asarray(s_um, float))/period_um+phi[:,None]
    return 2*np.sum((kval*tm)[:,None]*(1-np.cos(angle)), axis=0)+zero_coefficient(channel, transmission, leakage)*I0


def curve_harmonics(s_um, y, orders, period_um=16.):
    """Least-squares DC/cos/sin decomposition over any supplied scan interval."""
    s, y, orders = np.asarray(s_um,float), np.asarray(y,float), np.asarray(orders,int)
    angle = 4*np.pi*np.outer(s, orders)/period_um
    design = np.column_stack([np.ones(len(s)),np.cos(angle),np.sin(angle)])
    coeff, _, rank, singular = np.linalg.lstsq(design,y,rcond=None)
    n=len(orders)
    return dict(dc=float(coeff[0]), orders=orders.tolist(), cosine=coeff[1:n+1].tolist(),
        sine=coeff[n+1:].tolist(), rank=int(rank), columns=int(design.shape[1]),
        condition=float(singular[0]/singular[-1]),
        max_reconstruction_error=float(np.max(np.abs(design@coeff-y))))


def fit_global_parameters(records: Sequence[Mapping], *, free=("zero_transmission",), fixed=None,
                          bounds=None, allow_rank_deficient=False):
    """One parameter set across ALL supplied marks/channels; no per-mark offset.

    record: id, channel, target array, nonzero array, zero_power scalar. The
    nonzero curves are precomputed for one declared stop/K convention.
    Model: gain*[nonzero + (tau+L_chi*(1-tau))*zero_power]. Both gain and zero
    settings are shared. Scales balance records, not independently rescale them.
    tau,L_SUM,L_DIFF cannot all be identified from two effective coefficients;
    rank deficiency raises by default and must not be hidden by bounds.
    """
    defaults=dict(gain=1., zero_transmission=T0, leak_sum=LEAKAGE["SUM"], leak_diff=LEAKAGE["DIFF"])
    if fixed:
        defaults.update(fixed)
    allowed=set(defaults)
    if not set(free) <= allowed or len(set(free)) != len(free):
        raise ValueError("Unsupported or repeated free parameter")
    if not records:
        raise ValueError("Need at least one training record")
    default_bounds=dict(gain=(.5,1.5),zero_transmission=(0.,1.),leak_sum=(0.,.1),leak_diff=(0.,.1))
    if bounds:
        default_bounds.update(bounds)
    def unpack(x):
        return {**defaults,**dict(zip(free,map(float,x)))}
    def prediction(record,p):
        leak=p["leak_"+record["channel"].lower()]
        c=zero_coefficient(record["channel"],p["zero_transmission"],leak)
        return p["gain"]*(np.asarray(record["nonzero"],float)+c*float(record["zero_power"]))
    scales=[max(float(np.sqrt(np.mean(np.asarray(r["target"],float)**2))),1e-12) for r in records]
    def residual(x):
        p=unpack(x)
        return np.concatenate([(prediction(r,p)-np.asarray(r["target"],float))/scale/math.sqrt(len(r["target"]))
            for r,scale in zip(records,scales)])
    def jacobian(x):
        p=unpack(x); blocks=[]
        for r,scale in zip(records,scales):
            channel=r["channel"].lower();leak=p["leak_"+channel]
            pz=float(r["zero_power"]);nz=np.asarray(r["nonzero"],float)
            c=zero_coefficient(channel,p["zero_transmission"],leak)
            deriv={"gain":nz+c*pz,
                "zero_transmission":np.full_like(nz,p["gain"]*(1-leak)*pz),
                "leak_sum":np.full_like(nz,p["gain"]*(1-p["zero_transmission"])*pz if channel=="sum" else 0.),
                "leak_diff":np.full_like(nz,p["gain"]*(1-p["zero_transmission"])*pz if channel=="diff" else 0.)}
            blocks.append(np.column_stack([deriv[name] for name in free])/scale/math.sqrt(len(nz)))
        return np.vstack(blocks)
    x0=np.array([defaults[name] for name in free])
    if len(free):
        lo=np.array([default_bounds[name][0] for name in free]); hi=np.array([default_bounds[name][1] for name in free])
        result=least_squares(residual,x0,jac=jacobian,bounds=(lo,hi),xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=500)
        params=unpack(result.x); singular=np.linalg.svd(result.jac,compute_uv=False)
        rank=int(np.sum(singular > (singular[0]*1e-7 if len(singular) else 0)))
        if rank < len(free) and not allow_rank_deficient:
            raise ValueError(f"Global parameters unidentifiable: Jacobian rank {rank}/{len(free)}; fix T0 or leakage")
        success=bool(result.success)
    else:
        params=unpack([]); singular=np.array([]);rank=0;success=True
    diagnostics=[]
    for record in records:
        pred=prediction(record,params);truth=np.asarray(record["target"],float);err=pred-truth
        diagnostics.append(dict(id=str(record["id"]),channel=str(record["channel"]),
            max_abs_error=float(np.max(np.abs(err))),rmse=float(np.sqrt(np.mean(err**2))),
            dc_error=float(np.mean(err)), relative_l2=float(np.linalg.norm(err)/max(np.linalg.norm(truth),1e-30))))
    return dict(parameters=params,free=list(free),success=success,jacobian_rank=rank,
        jacobian_singular_values=singular.tolist(),records=diagnostics,
        shared_across_marks=True,per_mark_offsets=False,
        formula="gain * [nonzero + {tau + L_channel*(1-tau)} * zero_power]",
        identifiability="Only c_SUM and c_DIFF are identifiable if tau and both L are all free")


def evaluate_global_parameters(records, parameters):
    """Evaluate held-out marks with frozen training parameters, without fitting."""
    return fit_global_parameters(records,free=(),fixed=parameters)


def write_math_audit(output_directory=None):
    """Reproducible algebra checks; synthetic checks never claim new FDTD data."""
    out=Path(output_directory or Path(__file__).parent)
    out.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1]
    source=root/"inputs"/"D4C_FDTD_user_formula_20261105.md"
    k=np.linspace(-6e6,6e6,24001);G=2*np.pi/16e-6;w0=19e-6
    # Actual finite Gaussian shapes, deliberately not fitted or normalized by peak.
    a=np.sqrt(np.pi)*w0*(.55*np.exp(-(k*w0)**2/4)
        +(.16+.08j)*(np.exp(-((k-G)*w0)**2/4)+np.exp(-((k+G)*w0)**2/4))
        +(.09-.03j)*(np.exp(-((k-3*G)*w0)**2/4)+np.exp(-((k+3*G)*w0)**2/4)))
    pin=w0*math.sqrt(np.pi/2)/(2*Z0)
    s=np.linspace(-4,4,121)
    feats=b_features(k,a,pin,selected_orders=[1,3,5],stop_halfwidth_rad_um=.12)
    y=scan_from_features(feats,s,include_zero_pair=True)
    aw=angular_weights(k,pin);mask=stop_cell_field(k,.12);band=nearest_order(k)
    direct=[]
    allowed=np.isin(np.abs(band),[0,1,3,5])
    for ss in s:
        shifted=a*mask*np.exp(-1j*band*G*ss*1e-6)*allowed
        direct.append(.5*np.sum(aw*np.abs(shifted-shifted[::-1])**2))
    direct_error=float(np.max(np.abs(y-direct)))
    # Amplitude normalization and power normalization must cancel together.
    scale=.73-.41j
    fn=b_features(k,a/scale,pin/abs(scale)**2,selected_orders=[1,3,5],stop_halfwidth_rad_um=.12)
    norm_error=float(np.max(np.abs(scan_from_features(fn,s)-scan_from_features(feats,s))))
    # A synthetic global coefficient is recovered across marks, with a held-out
    # mark; this verifies the API, not a fit to any real numerical target.
    records=[]
    for j in range(4):
        nz=(.01+.007*j)*(1-np.cos(4*np.pi*s/16))
        pz=.1+.03*j
        records.append(dict(id=f"synthetic_{j}",channel="SUM",nonzero=nz,zero_power=pz,
            target=nz+zero_coefficient("SUM")*pz))
    fit=fit_global_parameters(records[:3])
    hold=evaluate_global_parameters(records[3:],fit["parameters"])
    checks=dict(pair_features_vs_full_angular_norm_max_error=direct_error,
        joint_field_power_normalization_max_error=norm_error,
        k0_sample_power_retained=feats.k0_sample_power,
        global_transmission_recovery_error=abs(fit["parameters"]["zero_transmission"]-T0),
        held_out_max_error=hold["records"][0]["max_abs_error"])
    a0=np.sqrt(np.pi)*w0*.55*np.exp(-(k*w0)**2/4)
    modal=b_features(k,a,pin,selected_orders=[1,3,5],stop_halfwidth_rad_um=.12,A_zero_ref_V=a0)
    accounting=zero_accounting(modal)
    checks["modal_stopped_power_closure_error"]=abs(modal.zero_pass_power+modal.nonzero_pass_power+modal.zero_nonzero_cross_pass_power-modal.collected_power)
    checks["double_counted_minus_scalar_equals_passed_error"]=abs(accounting["double_counted_zero"]-accounting["scalar_zero_once"]-modal.zero_pass_power)
    both=[]
    for r in records:
        for channel in ("SUM","DIFF"):
            copy=dict(r,channel=channel)
            copy["target"]=copy["nonzero"]+zero_coefficient(channel)*copy["zero_power"]
            both.append(copy)
    try:
        fit_global_parameters(both,free=("gain","zero_transmission","leak_sum","leak_diff"))
        raise AssertionError("A physically degenerate four-parameter fit was accepted")
    except ValueError as e:
        checks["four_parameter_degeneracy_rejected"]=str(e)
    assert direct_error < 1e-12 and norm_error < 1e-12
    assert checks["global_transmission_recovery_error"] < 1e-10
    audit=dict(source_file=str(source.relative_to(root)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_priority="User-supplied corrected formula, not the old illustrative Q network",
        T0=T0,L_SUM=LEAKAGE["SUM"],L_DIFF=LEAKAGE["DIFF"],
        c_SUM=zero_coefficient("SUM"),c_DIFF=zero_coefficient("DIFF"),
        K_confirmed=list(K_CONFIRMED),K_provisional=list(K_PROVISIONAL),K_missing=[],
        K_user_authorized_interpolation={str(m):K_VALUES[m] for m in K_INTERPOLATED},
        T_m_mapping_status="User confirmed T_m=|r_m|^2 in subsequent 2026-11-05 message",
        source_note_status="Original source file preserves pre-confirmation uncertainty; this audit records subsequent explicit user clarification",
        no_implicit_Tm_mapping=True,no_per_mark_offsets=True,no_free_per_order_weights=True,
        checks=checks,synthetic_only=True,
        zero_accounting_equations=dict(physical="P0_pass + L_chi*(P0_raw-P0_pass)",
            scalar="[T0 + L_chi*(1-T0)] * P0_raw",
            intentionally_repeated="P0_pass + [T0+L_chi*(1-T0)]*P0_raw"),
        normalization=dict(fourier="A(k)=integral E(x) exp(-ikx) dx",
            electric_spectrum_weight="gamma(k) dk / (4 pi Z0 Pin)",
            instrument_field_intensity_mapping="K_m/gamma_m in addition to electromagnetic gamma(k)",
            instrument_efficiency_mapping="K_m in addition to electromagnetic gamma(k)",
            k_zero="Single zero quadrature node is retained in all power norms; dark-pair amplitude is zero at that node"),
        finite_gaussian_controls=dict(waist_radius_um=19.,period_um=16.,periods=9,
            field_at_mark_edge=math.exp(-(72/19)**2),power_outside_mark=math.erfc(math.sqrt(2)*72/19),
            zero_power_outside_opaque_stop_halfwidth_0p2_rad_um=math.erfc(.2*19/math.sqrt(2)),
            qualification="Ideal Gaussian local-response estimates, not measured FDTD mode powers"),
        limitations=["DIFF interference sign not determined: SUM/DIFF vary leakage only here",
            "T_m=|r_m|^2 is user confirmed; callers must retain consistent E-field normalization",
            "B0 angular band is not the full physical zeroth-order mode",
            "A physical angular stop cannot generally implement only one modal attenuation",
            "gamma(k) electromagnetic weight and K_m empirical instrument weight are separate",
            "Bounds do not resolve tau/L_SUM/L_DIFF degeneracy",
            "Positive-k dark-pair prefactor is a model convention, not verified hardware throughput"])
    (out/"zero_order_math_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    return audit


if __name__ == "__main__":
    print(json.dumps(write_math_audit(),ensure_ascii=False,indent=2))
