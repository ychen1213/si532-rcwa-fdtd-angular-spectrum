"""TE RCWA with total reflected fields and passive complex substrate admittance.

Coordinates and phasors: exp(-i omega t), air above y=h, substrate below y=0,
downward incident exp(i q x-i ky y), upward reflected exp(i kx x+i ky y).
The periodic mother cell remains an infinite array. Its Gaussian illumination
is a coherent angular-spectrum superposition, not a finite-nine-period solve.
"""
from __future__ import annotations

import os
for _option in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_option, "2")

import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import eig, eigh, solve

ROOT = Path(__file__).resolve().parent
Z0 = 376.730313668


@dataclass(frozen=True)
class Config:
    wavelength_m: float = 643e-9
    period_m: float = 16e-6
    height_m: float = .3e-6
    ridge_width_m: float = 1.6e-6
    ridge_centers_m: tuple = (-3.2e-6, 0., 3.2e-6)
    air_n: complex = 1.+0j
    ridge_n: complex = 2.+0j
    substrate_n: complex = 4.+0j
    material_label: str = "controlled_n4_lossless_not_production"

    @property
    def k0(self): return 2*np.pi/self.wavelength_m

    @property
    def G(self): return 2*np.pi/self.period_m

    def metadata(self):
        result = asdict(self)
        for key in ("air_n", "ridge_n", "substrate_n"):
            result[key] = {"real": complex(getattr(self,key)).real,
                           "imag": complex(getattr(self,key)).imag}
        return result


def load_config(path):
    source = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    for key in ("air_n", "ridge_n", "substrate_n"):
        if key in source:
            n = source[key]
            source[key] = complex(n["real"], n.get("imag",0.)) if isinstance(n,dict) else complex(n)
    if "ridge_centers_m" in source:
        source["ridge_centers_m"] = tuple(source["ridge_centers_m"])
    return Config(**source)


def dump(path, obj):
    def serial(value):
        if isinstance(value, np.generic): return value.item()
        if isinstance(value, np.ndarray): return value.tolist()
        if isinstance(value, Path): return str(value)
        if isinstance(value, complex): return {"real":value.real,"imag":value.imag}
        raise TypeError(type(value).__name__)
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=serial), encoding="utf-8")


def outgoing_sqrt(value):
    """Passive outgoing branch: Im ky>=0, and Re ky>=0 when Im ky=0."""
    root = np.sqrt(np.asarray(value,dtype=complex))
    flip = (root.imag < 0) | ((root.imag == 0) & (root.real < 0))
    return np.where(flip,-root,root)


def epsilon_convolution(max_order, config=Config()):
    n = np.arange(-max_order,max_order+1)
    difference = n[:,None]-n[None,:]
    centers = np.asarray(config.ridge_centers_m)
    fill = config.ridge_width_m/config.period_m
    phase = np.exp(-1j*2*np.pi*difference[:,:,None]*centers[None,None,:]/config.period_m).sum(axis=2)
    eps = (config.ridge_n**2-config.air_n**2)*fill*np.sinc(difference*fill)*phase
    eps[difference==0] += config.air_n**2
    if np.max(abs(eps.imag)) < 1e-14:
        eps = eps.real
    return n,eps


def solve_angle(q, max_order=160, convolution=None, config=Config(),integrate_absorption=False):
    """Unit Ez incidence at y=0; interface fluxes normalized to incident flux.

    For a semi-infinite absorbing substrate, T_interface is power entering that
    substrate at y=0. It is ultimately absorbed there; it is NOT far-depth T.
    Evanescent y=0 r fields are represented by zero because y=0 intersects the
    patterned layer. The homogeneous top-plane coefficients retain all orders.
    """
    n,eps = epsilon_convolution(max_order,config) if convolution is None else convolution
    kx = q+n*config.G
    ya = outgoing_sqrt(config.air_n**2-(kx/config.k0)**2)
    ys = outgoing_sqrt(config.substrate_n**2-(kx/config.k0)**2)
    yi = outgoing_sqrt(config.air_n**2-(q/config.k0)**2).item()
    if yi.real <= 0 or abs(yi.imag)>1e-12:
        raise ValueError("The incident q must be propagating in a lossless superstrate")
    equation = np.asarray(eps).copy()
    equation[np.diag_indices_from(equation)] -= (kx/config.k0)**2
    hermitian = np.max(abs(equation-equation.conj().T)) < 1e-12
    if hermitian:
        eigenvalues,W = eigh(equation,check_finite=False,driver="evr")
        Winv = W.conj().T
    else:
        eigenvalues,W = eig(equation,check_finite=False)
        Winv = solve(W,np.eye(len(n)),check_finite=False)
    mode = outgoing_sqrt(eigenvalues)
    phase = config.k0*config.height_m*mode
    tiny = abs(phase)<1e-7
    ac = np.empty_like(mode);bc = np.empty_like(mode)
    ac[~tiny] = mode[~tiny]/np.tan(phase[~tiny])
    bc[~tiny] = mode[~tiny]/np.sin(phase[~tiny])
    ac[tiny] = 1/(config.k0*config.height_m)-(config.k0*config.height_m)*mode[tiny]**2/3
    bc[tiny] = 1/(config.k0*config.height_m)+(config.k0*config.height_m)*mode[tiny]**2/6
    A = (W*ac)@Winv
    B = (W*bc)@Winv
    lower = 1j*A.copy()
    lower[np.diag_indices_from(lower)] += ys
    transfer = solve(lower,1j*B,check_finite=False)
    Yin = 1j*A-1j*B@transfer
    boundary = Yin.copy();boundary[np.diag_indices_from(boundary)] += ya
    inc = np.zeros(len(n),complex);inc[n==0] = 1
    reflected = solve(boundary,ya*inc-Yin@inc,check_finite=False)
    transmitted = transfer@(inc+reflected)
    r_top = reflected*np.exp(-1j*yi*config.k0*config.height_m)
    r_y0 = reflected*np.exp(-1j*(ya+yi)*config.k0*config.height_m)
    t_y0 = transmitted*np.exp(-1j*yi*config.k0*config.height_m)
    ys0 = outgoing_sqrt(config.substrate_n**2-(q/config.k0)**2).item()
    r_fresnel = (yi-ys0)/(yi+ys0)
    background = np.zeros(len(n),complex);background[n==0] = r_fresnel
    propagating = (abs(ya.imag)<1e-12) & (ya.real>0)
    r_y0[~propagating] = 0
    R_orders = np.maximum(ya.real,0)*abs(reflected)**2/yi.real
    T_orders = np.maximum(ys.real,0)*abs(t_y0)**2/yi.real
    R = float(np.sum(R_orders));T = float(np.sum(T_orders))
    layer_absorption = 1-R-T
    lossless_layer = abs(config.ridge_n.imag)<1e-14 and abs(config.air_n.imag)<1e-14
    result=dict(order_n=n,kx_rad_m=kx,r_total_y0=r_y0,r_flat_y0=background,
                r_contrast_y0=r_y0-background,transmitted_y0=t_y0,
                r_total_layer_top_global_incident=r_top,
                reflected_order_power_fraction=R_orders,transmitted_interface_order_power_fraction=T_orders,
                reflected_power_fraction=R,transmitted_interface_power_fraction=T,
                layer_absorption_fraction=float(layer_absorption),
                semi_infinite_substrate_absorption_fraction=T if config.substrate_n.imag>0 else 0.,
                energy_error=abs(layer_absorption) if lossless_layer else max(-layer_absorption,0.),
                q_rad_m=float(q),substrate_ky_rad_m=config.k0*ys)
    if integrate_absorption:
        # Integrate omega*eps0/2 E* Im(eps) E independently through the layer.
        # Fourier orthogonality supplies the lateral per-period average.
        points,ww=np.polynomial.legendre.leggauss(64)
        y=(points+1)*config.height_m/2;ww=ww*config.height_m/2
        alpha=Winv@(inc+reflected);beta=Winv@transmitted
        upper=np.empty((len(mode),len(y)),complex);lower_part=np.empty_like(upper)
        upper[~tiny]=np.sin(mode[~tiny,None]*config.k0*y[None,:])/np.sin(phase[~tiny,None])
        lower_part[~tiny]=np.sin(mode[~tiny,None]*config.k0*(config.height_m-y[None,:]))/np.sin(phase[~tiny,None])
        upper[tiny]=y[None,:]/config.height_m;lower_part[tiny]=1-y[None,:]/config.height_m
        internal=W@(alpha[:,None]*upper+beta[:,None]*lower_part)
        eps_loss=(np.asarray(eps)-np.asarray(eps).conj().T)/(2j)
        density=np.einsum("iy,iy->y",internal.conj(),eps_loss@internal).real
        absorption=float(config.k0*np.dot(ww,density)/yi.real)
        result.update(layer_absorption_volume_integral_fraction=absorption,
                      independent_flux_volume_balance_abs_error=abs(R+T+absorption-1))
    return result


def build_library(config,max_order=320,q_max=450000.,q_step=2500.,folder=None):
    fingerprint=hashlib.sha256(json.dumps(config.metadata(),sort_keys=True).encode()).hexdigest()[:12]
    folder=Path(folder) if folder else ROOT/f"angle_library_N{max_order}_dq{int(q_step)}_{fingerprint}"
    folder.mkdir(parents=True,exist_ok=True)
    output=folder/"angle_library.npz"
    if output.exists() and (folder/"angle_library.csv").exists() and (folder/"metadata.json").exists():
        with np.load(output) as lib:
            if str(lib["config_sha256_12"])==fingerprint and float(lib["q_max_rad_m"])==q_max and float(lib["q_step_rad_m"])==q_step:
                return output
    centered_symmetric = np.allclose(np.sort(config.ridge_centers_m),-np.sort(config.ridge_centers_m)[::-1],atol=1e-15,rtol=0)
    if not centered_symmetric: raise ValueError("Negative-q symmetry requires a centered symmetric mother cell")
    qplus=np.arange(0,q_max+q_step*.5,q_step)
    n,eps=epsilon_convolution(max_order,config)
    rows=[];errors=[];started=time.perf_counter()
    cache_dir=ROOT/"point_cache"/f"{fingerprint}_N{max_order}"
    cache_dir.mkdir(parents=True,exist_ok=True)
    for i,q in enumerate(qplus):
        cache_file=cache_dir/f"q_{int(round(q))}.npz"
        if cache_file.exists():
            with np.load(cache_file) as cached:
                row={key:cached[key] for key in cached.files}
        else:
            row=solve_angle(q,max_order,(n,eps),config)
            np.savez_compressed(cache_file,**row)
        rows.append(row);errors.append(row["energy_error"])
        if i%10==0 or i==len(qplus)-1:
            print(f"N={max_order} q={i+1}/{len(qplus)} elapsed={time.perf_counter()-started:.1f}s balance={max(errors):.3g}",flush=True)
    q=np.r_[-qplus[:0:-1],qplus]
    values={}
    for key in ("r_total_y0","r_flat_y0","r_contrast_y0","transmitted_y0","r_total_layer_top_global_incident",
                "reflected_order_power_fraction","transmitted_interface_order_power_fraction","substrate_ky_rad_m"):
        positive=np.array([r[key] for r in rows])
        values[key]=np.concatenate([positive[:0:-1,::-1],positive])
    values.update(q_rad_m=q,order_n=n,q_max_rad_m=q_max,q_step_rad_m=q_step,
                  wavelength_m=config.wavelength_m,mother_period_m=config.period_m,layer_height_m=config.height_m,
                  phase_reference_y_m=0.,transverse_reference_x_m=0.,max_order=max_order,config_sha256_12=fingerprint,
                  energy_error_positive_q=np.array(errors),R_positive_q=np.array([r["reflected_power_fraction"] for r in rows]),
                  T_interface_positive_q=np.array([r["transmitted_interface_power_fraction"] for r in rows]),
                  substrate_n=config.substrate_n,ridge_n=config.ridge_n)
    np.savez_compressed(output,**values)
    with (folder/"angle_library.csv").open("w",newline="",encoding="utf-8-sig") as stream:
        writer=csv.writer(stream)
        writer.writerow(["q_rad_m","incident_theta_deg","order_m","output_kx_rad_m","output_propagating",
                         "r_total_re","r_total_im","r_flat_re","r_flat_im","r_contrast_re","r_contrast_im",
                         "r_total_layer_top_re","r_total_layer_top_im","t_interface_re","t_interface_im",
                         "R_order_fraction","T_interface_order_fraction","substrate_ky_re_rad_m","substrate_ky_im_rad_m"])
        for iq,qi in enumerate(q):
            for j,m in enumerate(n):
                rt=values["r_total_y0"][iq,j];rf=values["r_flat_y0"][iq,j];rc=values["r_contrast_y0"][iq,j]
                top=values["r_total_layer_top_global_incident"][iq,j];ts=values["transmitted_y0"][iq,j]
                ks=values["substrate_ky_rad_m"][iq,j];kx=qi+m*config.G
                writer.writerow([qi,np.degrees(np.arcsin(qi/config.k0)),m,kx,int(abs(kx)<config.k0),
                                 rt.real,rt.imag,rf.real,rf.imag,rc.real,rc.imag,top.real,top.imag,ts.real,ts.imag,
                                 values["reflected_order_power_fraction"][iq,j],values["transmitted_interface_order_power_fraction"][iq,j],ks.real,ks.imag])
    dump(folder/"metadata.json",dict(solver="Independent TE RCWA with passive complex Si substrate",config=config.metadata(),
         full_harmonic_count=2*max_order+1,total_q_count=len(q),independent_positive_q_solutions=len(qplus),
         elapsed_s=time.perf_counter()-started,passivity_violation_max=max(errors),
         layer_absorption_positive_q_min=float(min(r["layer_absorption_fraction"] for r in rows)),
         layer_absorption_positive_q_max=float(max(r["layer_absorption_fraction"] for r in rows)),
         physical_input_field="r_total: includes coherent flat-substrate reflection, especially m=0",
         subtraction_field="r_contrast is diagnostic only; cannot represent actual light before the zero-order stop",
         loss_scope="T_interface is inward y=0 flux. For an absorbing semi-infinite substrate it ultimately becomes substrate absorption.",
         evanescent_scope="All evanescent coefficients are retained at y=height, while y=0 far-field arrays store zero outside air light cone.",
         phase_convention="exp(-i omega t), upward exp(+i ky y), Ez incident at y=0 is unity",
         negative_q_reconstruction="r_m(-q)=r_-m(q) for centered symmetric cell"))
    return output


def gaussian_spectrum(q,w0_m,beam_center_m=0.):
    return np.sqrt(np.pi)*w0_m*np.exp(-(np.asarray(q)*w0_m)**2/4)*np.exp(-1j*np.asarray(q)*beam_center_m)


def library_spectrum(library,kx,w0_m=19e-6,beam_center_m=0.,shift_m=0.,field_kind="total",NA=.6,
                     measured_q=None,measured_D0=None,return_orders=False,frozen_angle=False):
    """Coherent A(k)=sum r_m(k-mG) D(k-mG) exp(-imG d); default TOTAL."""
    close=isinstance(library,(str,Path))
    lib=np.load(library) if close else library
    q=lib["q_rad_m"];n=lib["order_n"];r=lib[f"r_{field_kind}_y0"]
    wl=float(lib["wavelength_m"]);G=2*np.pi/float(lib["mother_period_m"]);k0=2*np.pi/wl
    useful=n[abs(n*G)<NA*k0+float(lib["q_max_rad_m"])]
    pieces=[]
    for m in useful:
        qi=np.asarray(kx)-m*G;j=int(np.flatnonzero(n==m)[0])
        coeff=(np.full(len(kx),r[len(q)//2,j],complex) if frozen_angle else
               np.interp(qi,q,r[:,j].real,left=0,right=0)+1j*np.interp(qi,q,r[:,j].imag,left=0,right=0))
        D=(gaussian_spectrum(qi,w0_m,beam_center_m) if measured_q is None else
           np.interp(qi,measured_q,np.real(measured_D0),left=0,right=0)+1j*np.interp(qi,measured_q,np.imag(measured_D0),left=0,right=0))
        pieces.append(coeff*D*np.exp(-1j*m*G*shift_m))
    pieces=np.asarray(pieces)
    if close: lib.close()
    return (pieces.sum(axis=0),useful,pieces) if return_orders else pieces.sum(axis=0)


def reconstruct_field(library,x,w0_m=19e-6,beam_center_m=0.,shift_m=0.,monitor_y_m=1.65e-6,
                      field_kind="total",measured_q=None,measured_D0=None):
    """Direct q integration at a homogeneous air plane, with all evanescent orders."""
    with np.load(library) as lib:
        h=float(lib["layer_height_m"]);k0=2*np.pi/float(lib["wavelength_m"]);G=2*np.pi/float(lib["mother_period_m"])
        if monitor_y_m<h: raise ValueError("Physical monitor must lie above patterned layer")
        q=lib["q_rad_m"];n=lib["order_n"];top=lib["r_total_layer_top_global_incident"];flat=lib["r_flat_y0"]
        D=(gaussian_spectrum(q,w0_m,beam_center_m) if measured_q is None else
           np.interp(q,measured_q,np.real(measured_D0),left=0,right=0)+1j*np.interp(q,measured_q,np.imag(measured_D0),left=0,right=0))
        weights=np.ones(len(q))*float(lib["q_step_rad_m"]);weights[[0,-1]]*=.5
        eq=np.exp(1j*q[:,None]*np.asarray(x)[None,:]);field=np.zeros(len(x),complex)
        for j,m in enumerate(n):
            kx=q+m*G;ky=outgoing_sqrt(k0**2-kx**2)
            rr=top[:,j]*np.exp(1j*ky*(monitor_y_m-h))
            ff=flat[:,j]*np.exp(1j*ky*monitor_y_m)
            if field_kind=="contrast": rr=rr-ff
            elif field_kind=="flat": rr=ff
            elif field_kind!="total": raise ValueError(field_kind)
            c=rr*D*weights/(2*np.pi)*np.exp(-1j*m*G*shift_m)
            if np.sum(abs(c))<1e-15:continue
            field+=(c@eq)*np.exp(1j*m*G*np.asarray(x))
        return field


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--config",required=True)
    parser.add_argument("--max-order",type=int,default=320)
    parser.add_argument("--q-max",type=float,default=450000.)
    parser.add_argument("--q-step",type=float,default=2500.)
    parser.add_argument("--folder")
    args=parser.parse_args()
    print(build_library(load_config(args.config),args.max_order,args.q_max,args.q_step,args.folder),flush=True)


if __name__=="__main__":main()
