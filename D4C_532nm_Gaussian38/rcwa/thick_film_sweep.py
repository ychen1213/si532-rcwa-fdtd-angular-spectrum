"""Extra held-out weak-scan mark; preserve the original eight-model contract.

Screen cap300..1000nm on the ORIGINAL h300 AH53 by reusing its patterned
layer input admittance. A uniform-cap scattering recurrence needs one vector
linear solve per thickness. The extra selected h180/cap400 mark is then solved
independently through N640 with the main RCWA solver and identical D4C fields.
"""
from pathlib import Path
import argparse
from dataclasses import replace
import json,time
import numpy as np
from scipy.linalg import solve
import gaussian_angle_controls as angular
from te_planarized_rcwa import outgoing_sqrt,dump
from run_reference import ROOT,HERE,make_config,solve_cached,compare,csvwrite,sha,K_VALUES,K_PROVISIONAL,T0,L

M=np.arange(1,19);K=np.array([{**K_VALUES,**K_PROVISIONAL}[int(m)] for m in M])
S=np.arange(-800,801)*.005
KERNEL=1-np.cos(4*np.pi*M[:,None]*S[None,:]/16.)


def cap_setup(c,N):
    angular.N=N
    n,kx,ya,yi,Y=angular.pattern(c,0.)
    yc=outgoing_sqrt(c.cap_n**2-(kx/c.k0)**2)
    Rfilm=solve(np.diag(yc)+Y,np.diag(yc)-Y,check_finite=False)
    rac=(ya-yc)/(ya+yc);tac=2*ya/(ya+yc);tca=2*yc/(ya+yc)
    inc=np.zeros(len(n),complex);inc[n==0]=1
    return n,kx,ya,yi,yc,Rfilm,rac,tac,tca,inc


def cap_response(c,p):
    n,kx,ya,yi,yc,Rfilm,rac,tac,tca,inc=p
    P=np.exp(1j*c.k0*yc*c.cap_m)
    Q=P[:,None]*Rfilm*P[None,:]
    forward=solve(np.eye(len(n))+rac[:,None]*Q,tac*inc,check_finite=False)
    rb=rac*inc+tca*(Q@forward)
    rr=rb*np.exp(-1j*(ya+yi)*c.k0*c.top_m)
    rr[(abs(ya.imag)>1e-12)|(ya.real<=0)]=0
    return n,rr,float(np.sum(np.maximum(ya.real,0)*abs(rb)**2)/yi.real)


def target_stats(n,r,R,c):
    j0=int(np.flatnonzero(n==0)[0]);rp=r[np.searchsorted(n,M)];rn=r[np.searchsorted(n,-M)]
    I0=float(abs(r[j0])**2);phi=np.angle(rp*np.conj(rn))
    # All geometries remain mirror symmetric. Keep actual pair phase in final data.
    ac=(2*K*abs(rp)**2)@(1-np.cos(4*np.pi*M[:,None]*S[None,:]/16.+phi[:,None]))
    zero=(T0+L['SUM']*(1-T0))*I0
    return dict(height_nm=c.height_m*1e9,cap_nm=c.cap_m*1e9,I0=I0,R=R,
        S_SUM_min=float(ac.min()+zero),S_SUM_peak=float(ac.max()+zero),
        S_SUM_mean=float(ac[:-1].mean()+zero),zero_SUM=zero)


def extra_target(folder,v,c,r):
    n=r['order_n'];rr=r['r_total_y0'];rp=rr[np.searchsorted(n,M)];rn=rr[np.searchsorted(n,-M)]
    I0=float(abs(rr[np.searchsorted(n,0)])**2);tm=abs(rp)**2;phi=np.angle(rp*np.conj(rn))
    cm=2*tm[:,None]*(1-np.cos(4*np.pi*M[:,None]*S[None,:]/16.+phi[:,None]))
    known=sum(K_VALUES[int(m)]*cm[j] for j,m in enumerate(M) if m in K_VALUES)
    ac=known+.4947*cm[1]+.49135*cm[3]
    zs=(T0+L['SUM']*(1-T0))*I0;zd=(T0+L['DIFF']*(1-T0))*I0
    np.savez_compressed(folder/'d4c_targets.npz',s_um=S,order_positive=M,r_plus=rp,r_minus=rn,T_m=tm,T_m_working_definition=tm,
        phi_m_rad=phi,I0=I0,T0=T0,AC=ac,S_SUM=ac+zs,S_DIFF_same_kernel=ac+zd,
        AC_known_K=known,coefficient_K2=cm[1],coefficient_K4=cm[3],AC_provisional_interpolated_K2K4=ac,
        zero_SUM=zs,zero_DIFF=zd,S_SUM_known_K=known+zs,S_DIFF_same_kernel_known_K=known+zd,
        S_SUM_provisional_interpolated_K2K4=ac+zs,S_DIFF_same_kernel_provisional_interpolated_K2K4=ac+zd)
    csvwrite(folder/'d4c_targets.csv',[dict(s_um=float(s),AC=float(a),S_SUM=float(a+zs),S_DIFF_same_kernel=float(a+zd),
        AC_known_K=float(kn),coefficient_K2=float(c2),coefficient_K4=float(c4),zero_SUM=zs,zero_DIFF=zd,
        S_SUM_known_K=float(kn+zs),S_DIFF_same_kernel_known_K=float(kn+zd),
        S_SUM_provisional_interpolated_K2K4=float(a+zs),S_DIFF_same_kernel_provisional_interpolated_K2K4=float(a+zd))
        for s,a,kn,c2,c4 in zip(S,ac,known,cm[1],cm[3])])
    np.savez_compressed(folder/'reference.npz',**r,
        r_common_y085=r['r_total_y0']*np.exp(1j*np.sqrt(np.maximum(0,c.k0**2-r['kx_rad_m']**2))*.85e-6),
        config_json=json.dumps(c.metadata(),sort_keys=True),material_source_sha256=sha(ROOT/'fdtd/materials.json'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--candidate-only',action='store_true');args=parser.parse_args()
    materials=json.loads((ROOT/'fdtd/materials.json').read_text(encoding='utf8'))
    original=json.loads((ROOT/'experiment.json').read_text(encoding='utf8'))['variants'][0]
    v=dict(original,id='ah53_h180_cap400',ridge_height_nm=180.,cap_nm=400.)
    dump(ROOT/'inputs/extra_weak_mark.json',dict(variant=v,
        objective='User means scan vertical-axis peak .00x, not I0 or total reflection R.',
        split='Extra held-out structure only; original cap0/100 fits and original eight variants remain fixed.',
        physical_change='Both Si relief height300->180nm and planar SiO2 cap400nm; not a film-only effect.',
        fit_policy='Apply the previously frozen shared gain/T0. Do not include this mark in any calibration.',
        material_source='fdtd/materials.json',same_source_and_monitor_geometry=True))
    c=make_config(v,materials);history=[];previous=None
    # The selected extra mark is prioritized before the optional thickness sweep.
    for N in (160,240,320,480,640):
        r,meta=solve_cached(v['id'],c,N)
        item=dict(max_order=N,full_harmonic_count=2*N+1,**target_stats(r['order_n'],r['r_total_y0'],r['R'],c))
        if previous is not None:item.update(compare(r,previous))
        history.append(item);previous=r
    folder=HERE/'cases'/v['id'];extra_target(folder,v,c,r)
    csvwrite(folder/'convergence.csv',history)
    dump(HERE/'thick_film_candidate.json',dict(variant=v,selected_reference=(folder/'reference.npz').relative_to(ROOT).as_posix(),
        d4c_targets=(folder/'d4c_targets.npz').relative_to(ROOT).as_posix(),stats=history[-1],convergence=history,
        T0=T0,K2=.4947,K4=.49135,split='extra frozen-parameter heldout; not part of original eight fits',
        materials_sha256=sha(ROOT/'fdtd/materials.json')))
    print('CANDIDATE_FINAL',json.dumps(history[-1]),flush=True)
    if args.candidate_only and (HERE/'thick_film_sweep.json').exists():
        old=json.loads((HERE/'thick_film_sweep.json').read_text(encoding='utf8'))
        old['selected_candidate_max_order']=N
        old['limitation']=f'Thickness-grid survey, not a proof of the global continuous minimum; N160 screening, N{N} selected extra mark.'
        dump(HERE/'thick_film_sweep.json',old)
        return
    # Original h300: shared-pattern sweep, covering all requested cap thicknesses.
    base=dict(original,cap_nm=300.)
    bc=make_config(base,materials);p=cap_setup(bc,160)
    rows=[];started=time.perf_counter()
    for cap in np.arange(300.,1000.01,10.):
        cc=replace(bc,cap_m=cap*1e-9)
        n,rr,R=cap_response(cc,p);rows.append(target_stats(n,rr,R,cc))
    csvwrite(HERE/'thick_film_sweep.csv',rows)
    # Independent recurrence checks at three already computed cap thicknesses.
    checks=[]
    for cap in (100.,150.,200.):
        cc=replace(bc,cap_m=cap*1e-9);n,rr,R=cap_response(cc,p)
        with np.load(HERE/'cases'/f'ah53_cap{int(cap):03d}'/'N160_q0.npz') as z:
            err=float(np.max(abs(rr-z['r_total_y0'])))
            checks.append(dict(cap_nm=cap,max_abs_complex_recurrence_vs_modal=err,passed=err<1e-10))
    if not all(q['passed'] for q in checks):raise RuntimeError('Cap recurrence did not match the full solver')
    out=dict(original_height_nm=300,cap_range_nm=[300,1000],cap_step_nm=10,max_order=160,
        elapsed_s=time.perf_counter()-started,min_scan_peak=min(rows,key=lambda x:x['S_SUM_peak']),
        min_I0=min(rows,key=lambda x:x['I0']),min_total_R=min(rows,key=lambda x:x['R']),
        meets_user_scan_peak_001_to_009=any(.001<=x['S_SUM_peak']<=.009 for x in rows),
        all_sweep_rows=rows,recurrence_checks=checks,
        limitation='Thickness-grid survey, not a proof of the global continuous minimum; N160 screening, N640 selected extra mark.',
        original_eight_configuration_unchanged=True)
    dump(HERE/'thick_film_sweep.json',out)
    print('ORIGINAL_H300_SWEEP',json.dumps({k:v for k,v in out.items() if k not in ('all_sweep_rows','recurrence_checks')}),flush=True)


if __name__=='__main__':main()
