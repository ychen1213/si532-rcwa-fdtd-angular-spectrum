"""TE RCWA for a Si relief with optional planar SiO2 groove fill and cap.

Phasors exp(-i omega t), y points upwards, Ez scalar polarization.
Incident electric amplitude at y=0 is one; reported r_y0 shares that plane.
The original, previously validated single-layer solver is retained in the
provenance directory; the added uniform cap is another exact modal layer.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '2'
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eig, eigh, solve


@dataclass(frozen=True)
class Config:
    wavelength_m: float = 532e-9
    period_m: float = 16e-6
    height_m: float = 300e-9
    ridge_width_m: float = 1.6e-6
    ridge_centers_m: tuple = (-3.2e-6, 0., 3.2e-6)
    cap_m: float = 0.
    air_n: complex = 1.+0j
    ridge_n: complex = 4.15+.05j
    substrate_n: complex = 4.15+.05j
    fill_n: complex = 1.+0j
    cap_n: complex = 1.46+0j

    @property
    def k0(self): return 2*np.pi/self.wavelength_m
    @property
    def G(self): return 2*np.pi/self.period_m
    @property
    def top_m(self): return self.height_m+self.cap_m

    def metadata(self):
        d=asdict(self)
        for k in ('air_n','ridge_n','substrate_n','fill_n','cap_n'):
            z=complex(d[k]);d[k]={'real':z.real,'imag':z.imag}
        return d


def outgoing_sqrt(value):
    root=np.sqrt(np.asarray(value,complex))
    return np.where((root.imag<0)|((root.imag==0)&(root.real<0)),-root,root)


def epsilon_convolution(max_order,config):
    n=np.arange(-max_order,max_order+1)
    d=n[:,None]-n[None,:]
    f=config.ridge_width_m/config.period_m
    centers=np.asarray(config.ridge_centers_m)
    phase=np.exp(-1j*2*np.pi*d[:,:,None]*centers[None,None,:]/config.period_m).sum(axis=2)
    eps=(config.ridge_n**2-config.fill_n**2)*f*np.sinc(d*f)*phase
    eps[d==0]+=config.fill_n**2
    return n,eps


def modal_cotangent(mode,thickness_k0):
    phase=thickness_k0*mode
    tiny=abs(phase)<1e-7
    a=np.empty_like(mode);b=np.empty_like(mode)
    a[~tiny]=mode[~tiny]/np.tan(phase[~tiny])
    b[~tiny]=mode[~tiny]/np.sin(phase[~tiny])
    a[tiny]=1/thickness_k0-thickness_k0*mode[tiny]**2/3
    b[tiny]=1/thickness_k0+thickness_k0*mode[tiny]**2/6
    return a,b


def solve_angle(q,max_order,config):
    n,eps=epsilon_convolution(max_order,config)
    size=len(n);ix=np.diag_indices(size)
    kx=q+n*config.G
    ya=outgoing_sqrt(config.air_n**2-(kx/config.k0)**2)
    ys=outgoing_sqrt(config.substrate_n**2-(kx/config.k0)**2)
    yi=outgoing_sqrt(config.air_n**2-(q/config.k0)**2).item()
    if yi.real<=0 or abs(yi.imag)>1e-12: raise ValueError('Incident wave must propagate.')
    equation=eps.copy();equation[ix]-=(kx/config.k0)**2
    if np.max(abs(equation-equation.conj().T))<1e-12:
        ev,W=eigh(equation,check_finite=False,driver='evr');Wi=W.conj().T
    else:
        ev,W=eig(equation,check_finite=False);Wi=solve(W,np.eye(size),check_finite=False)
    mode=outgoing_sqrt(ev)
    ac,bc=modal_cotangent(mode,config.k0*config.height_m)
    A=(W*ac)@Wi;B=(W*bc)@Wi
    lower=1j*A.copy();lower[ix]+=ys
    transfer_pattern=solve(lower,1j*B,check_finite=False)
    Ypattern=1j*A-1j*B@transfer_pattern
    if config.cap_m>0:
        yc=outgoing_sqrt(config.cap_n**2-(kx/config.k0)**2)
        ca,cb=modal_cotangent(yc,config.k0*config.cap_m)
        lower=Ypattern.copy();lower[ix]+=1j*ca
        transfer_cap=solve(lower,np.diag(1j*cb),check_finite=False)
        Yin=np.diag(1j*ca)-(1j*cb[:,None])*transfer_cap
    else:
        transfer_cap=None;Yin=Ypattern
    boundary=Yin.copy();boundary[ix]+=ya
    inc=np.zeros(size,complex);inc[n==0]=1
    reflected=solve(boundary,ya*inc-Yin@inc,check_finite=False)
    field_pattern_top=(inc+reflected) if transfer_cap is None else transfer_cap@(inc+reflected)
    transmitted=transfer_pattern@field_pattern_top
    global_phase=np.exp(-1j*yi*config.k0*config.top_m)
    r_top=reflected*global_phase
    r_y0=reflected*np.exp(-1j*(ya+yi)*config.k0*config.top_m)
    t_y0=transmitted*global_phase
    propagating=(abs(ya.imag)<1e-12)&(ya.real>0)
    r_y0[~propagating]=0
    Rorders=np.maximum(ya.real,0)*abs(reflected)**2/yi.real
    Torders=np.maximum(ys.real,0)*abs(t_y0)**2/yi.real
    R=float(Rorders.sum());T=float(Torders.sum())
    return dict(order_n=n,q_rad_m=float(q),kx_rad_m=kx,r_total_y0=r_y0,
        r_total_top_global_incident=r_top,transmitted_y0=t_y0,
        reflected_order_power_fraction=Rorders,transmitted_interface_order_power_fraction=Torders,
        R=R,T_interface=T,finite_layer_absorption=1-R-T,
        max_order=max_order,full_harmonic_count=size,wavelength_m=config.wavelength_m,
        period_m=config.period_m,reference_y_m=0.,top_y_m=config.top_m)


def dump(path,obj):
    def serial(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,complex):return {'real':v.real,'imag':v.imag}
        if isinstance(v,Path):return str(v)
        raise TypeError(type(v).__name__)
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=serial),encoding='utf8')
