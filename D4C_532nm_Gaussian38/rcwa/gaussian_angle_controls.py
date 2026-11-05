"""Small-angle TE controls with a shared patterned-layer response per cap family.

N160 angular diagnostics compare r(q) against r(0) at the SAME truncation.
They quantify the frozen-q0 approximation, not finite-nine-period edge error.
The positive-q controls reuse one patterned-layer eigensystem across all three
uniform cap thicknesses. Negative q is filled by exact mirror symmetry.
"""
from pathlib import Path
from dataclasses import replace
import math,json,time
import numpy as np
from scipy.linalg import eig,eigh,solve
from scipy.special import erfcinv
from te_planarized_rcwa import epsilon_convolution,outgoing_sqrt,modal_cotangent,dump
from run_reference import make_config,csvwrite,T0,sha

ROOT=Path(__file__).resolve().parents[1];HERE=ROOT/'rcwa';N=160


def pattern(c,q):
    n,eps=epsilon_convolution(N,c);sz=len(n);ix=np.diag_indices(sz)
    kx=q+n*c.G;ya=outgoing_sqrt(c.air_n**2-(kx/c.k0)**2)
    ys=outgoing_sqrt(c.substrate_n**2-(kx/c.k0)**2)
    yi=outgoing_sqrt(c.air_n**2-(q/c.k0)**2).item()
    M=eps.copy();M[ix]-=(kx/c.k0)**2
    if np.max(abs(M-M.conj().T))<1e-12:
        ev,W=eigh(M,check_finite=False,driver='evr');Wi=W.conj().T
    else:
        ev,W=eig(M,check_finite=False);Wi=solve(W,np.eye(sz),check_finite=False)
    ac,bc=modal_cotangent(outgoing_sqrt(ev),c.k0*c.height_m)
    A=(W*ac)@Wi;B=(W*bc)@Wi;lower=1j*A.copy();lower[ix]+=ys
    transfer=solve(lower,1j*B,check_finite=False)
    return n,kx,ya,yi,1j*A-1j*B@transfer


def finish(c,q,p):
    n,kx,ya,yi,Ypattern=p;sz=len(n);ix=np.diag_indices(sz)
    if c.cap_m:
        yc=outgoing_sqrt(c.cap_n**2-(kx/c.k0)**2)
        ca,cb=modal_cotangent(yc,c.k0*c.cap_m)
        lower=Ypattern.copy();lower[ix]+=1j*ca
        trans=solve(lower,np.diag(1j*cb),check_finite=False)
        Yin=np.diag(1j*ca)-(1j*cb[:,None])*trans
    else:Yin=Ypattern
    inc=np.zeros(sz,complex);inc[n==0]=1
    bnd=Yin.copy();bnd[ix]+=ya
    rr=solve(bnd,ya*inc-Yin@inc,check_finite=False)
    rr*=np.exp(-1j*(ya+yi)*c.k0*c.top_m)
    rr[(abs(ya.imag)>1e-12)|(ya.real<=0)]=0
    return rr


def main():
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf8'))
    mat=json.loads((ROOT/'fdtd/materials.json').read_text(encoding='utf8'))
    configs={v['id']:make_config(v,mat) for v in plan['variants']}
    qp=np.arange(0,200000.1,25000.);order=np.arange(-N,N+1)
    results={k:[] for k in configs};timing=[]
    cached=all((HERE/'cases'/k/'angular_controls_N160.npz').exists() for k in configs)
    if cached:
        previous=HERE/'gaussian_angle_controls.json'
        if previous.exists():timing=json.loads(previous.read_text(encoding='utf8')).get('timings',[])
        for name in configs:
            with np.load(HERE/'cases'/name/'angular_controls_N160.npz') as z:
                if not np.allclose(z['q_rad_m'],np.r_[-qp[:0:-1],qp]):raise ValueError('Cached angular grid mismatch.')
                results[name]=z['r_total_y0'][len(qp)-1:]
    else:
        for geometry in ('ah53','binary50'):
            for fill in (False,True):
                names=[k for k in configs if k.startswith(geometry+'_') and bool(configs[k].cap_m)==fill]
                for q in qp:
                    start=time.perf_counter();p=pattern(configs[names[0]],q)
                    for name in names:results[name].append(finish(configs[name],q,p))
                    elapsed=time.perf_counter()-start
                    timing.append(dict(geometry=geometry,fill=fill,q_rad_m=float(q),cap_variants=len(names),elapsed_s=elapsed))
                    print(f'Angle control {geometry} fill={fill} q={q:g} {elapsed:.2f}s',flush=True)
    q=np.r_[-qp[:0:-1],qp];w0=19e-6
    # Incoming power measure gamma_i(q)|D(q)|² dq. Overall constants cancel.
    weight=np.exp(-.5*(q*w0)**2)*np.sqrt(1-(q/configs[next(iter(configs))].k0)**2)
    rows=[]
    for name,rr in results.items():
        positive=np.asarray(rr);r=np.concatenate([positive[:0:-1,::-1],positive])
        center=len(q)//2;base=r[center]
        selected=abs(order)<=18
        residual=np.sum(abs(r[:,selected]-base[None,selected])**2,axis=1)
        denom=np.sum(abs(r[:,selected])**2,axis=1)
        gL2=float(np.sqrt(np.trapezoid(weight*residual,q)/np.trapezoid(weight*denom,q)))
        nonzero=selected&(order!=0)
        nzres=np.sum(abs(r[:,nonzero]-base[None,nonzero])**2,axis=1)
        nzden=np.sum(abs(r[:,nonzero])**2,axis=1)
        nzL2=float(np.sqrt(np.trapezoid(weight*nzres,q)/np.trapezoid(weight*nzden,q)))
        j0=int(np.flatnonzero(order==0)[0]);I0avg=float(np.trapezoid(weight*abs(r[:,j0])**2,q)/np.trapezoid(weight,q))
        I0q0=float(abs(base[j0])**2)
        folder=HERE/'cases'/name
        np.savez_compressed(folder/'angular_controls_N160.npz',q_rad_m=q,order_n=order,r_total_y0=r,
            Gaussian_power_weight=weight,w0_m=w0,max_order=N,
            scope='Small-angle frozen-q0 diagnostic at N160; not a convergence-certified full Gaussian scattering library.')
        rows.append(dict(variant=name,max_order=N,positive_q_solves=len(qp),negative_q_from_mirror=True,
            Gaussian_weighted_complex_L2_frozen_q0=gL2,Gaussian_weighted_nonzero_complex_L2_frozen_q0=nzL2,
            I0_q0=I0q0,I0_Gaussian_angle_average=I0avg,
            I0_angle_average_relative_change=(I0avg-I0q0)/I0q0,
            excluded_Gaussian_incident_power_fraction=math.erfc(200000*w0/math.sqrt(2))))
    csvwrite(HERE/'gaussian_angle_controls.csv',rows)
    optical=dict(T0=float(T0),w0_m=w0,
        matched_opaque_stop_halfwidth_rad_m=float(math.sqrt(2)*erfcinv(T0)/w0),
        old_stop_halfwidth_rad_m=200000.,old_stop_Gaussian_tail_power_fraction=math.erfc(200000*w0/math.sqrt(2)),
        m0_nearest_bin_halfwidth_rad_m=math.pi/16e-6,
        definition='Ideal isolated Gaussian zero-order power fraction outside an opaque symmetric stop. Not the measured finite FDTD zero spectrum.')
    dump(HERE/'gaussian_zero_stop_control.json',optical)
    dump(HERE/'gaussian_angle_controls.json',dict(structures=rows,optical=optical,timings=timing,
        materials_sha256=sha(ROOT/'fdtd/materials.json'),
        limitations=['N160 is a same-grid angular-sensitivity diagnostic. Primary q0 has independent higher-harmonic convergence.',
            '25k rad/m q grid, ±200k support, trapezoid. No angular sampling convergence claim.',
            'Gaussian spectrum belongs to an infinite periodic array, not finite nine-period FDTD.',
            'Coefficient L2 uses incoming gamma_i-weighted Gaussian power; the I0 moment is the m0 angular power ratio because k_out=q.',
            'Negative q is reconstructed from exact mirror symmetry of the centered materials geometry.']))


if __name__=='__main__':main()
