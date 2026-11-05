"""Same-N angular-sensitivity diagnostic of the ninth weak-scan mark.

It does NOT replace the primary normal-incidence D4C target. The symmetric
Gaussian parity-paired packet has r_m(q) r_-m(-q)* = |r_m(q)|², so the recorded
coefficient moment is useful for diagnosing angular illumination differences.
No collection stop or extra empirical NA weighting is introduced in this test.
"""
from pathlib import Path
import json,math
import numpy as np
from run_reference import ROOT,HERE,make_config,solve_cached,csvwrite,T0,L,sha
from thick_film_sweep import M,K,S
from te_planarized_rcwa import dump

N=160;ID='ah53_h180_cap400'


def diagnostic(q,r,n,c,label):
    w0=19e-6;weight=np.exp(-.5*(q*w0)**2)*np.sqrt(1-(q/c.k0)**2)
    selected=np.searchsorted(n,M);j0=int(np.flatnonzero(n==0)[0]);center=int(np.flatnonzero(q==0)[0])
    tm=abs(r[:,selected])**2;tm0=tm[center]
    mean=np.trapezoid(weight[:,None]*tm,q,axis=0)/np.trapezoid(weight,q)
    I0=float(abs(r[center,j0])**2);Iavg=float(np.trapezoid(weight*abs(r[:,j0])**2,q)/np.trapezoid(weight,q))
    kernel=1-np.cos(4*np.pi*M[:,None]*S[None,:]/16.)
    y0=(2*K*tm0)@kernel;ya=(2*K*mean)@kernel
    dev=ya-y0;yc=y0[:-1]-y0[:-1].mean();ac=ya[:-1]-ya[:-1].mean()
    stats=dict(label=label,max_order=N,q_max_rad_m=float(q.max()),q_nodes=len(q),
        q0_AC_peak=float(y0.max()),angle_averaged_AC_peak=float(ya.max()),
        AC_relative_L2=float(np.linalg.norm(dev)/np.linalg.norm(y0)),
        centered_AC_relative_L2=float(np.linalg.norm(ac-yc)/np.linalg.norm(yc)),
        AC_span_relative_change=float(np.ptp(ya)/np.ptp(y0)-1),
        AC_peak_relative_change=float(ya.max()/y0.max()-1),
        I0_q0=I0,I0_averaged=Iavg,I0_relative_change=Iavg/I0-1,
        power_outside_input_q_range=math.erfc(q.max()*w0/math.sqrt(2)),
        input_power_weight='gamma_in(q)*exp[-(q*w0)^2/2]',
        comparison='Frozen q0 and angular average use identical N160; neither is the primary N640 target.')
    rows=[dict(m=int(m),K=float(kk),T_q0=float(t0),T_Gaussian_mean=float(tavg),
        relative_change=float(tavg/t0-1) if t0 else None,absolute_change=float(tavg-t0),
        weighted_AC_DC_q0=float(2*kk*t0),weighted_AC_DC_mean=float(2*kk*tavg)) for m,kk,t0,tavg in zip(M,K,tm0,mean)]
    return stats,rows,dict(s_um=S,AC_frozen_q0=y0,AC_Gaussian_average=ya,
        S_SUM_frozen_q0=y0+(T0+L['SUM']*(1-T0))*I0,
        S_SUM_Gaussian_average=ya+(T0+L['SUM']*(1-T0))*Iavg)


def main():
    v=json.loads((ROOT/'inputs/extra_weak_mark.json').read_text(encoding='utf8'))['variant']
    materials=json.loads((ROOT/'fdtd/materials.json').read_text(encoding='utf8'));c=make_config(v,materials)
    positive={};n=None
    for q in (0.,50000.,100000.):
        r,_=solve_cached(ID,c,N,q);positive[q]=r['r_total_y0'];n=r['order_n']
    q=np.array([-100000.,-50000.,0.,50000.,100000.])
    r=np.array([positive[abs(qq)][::-1] if qq<0 else positive[qq] for qq in q])
    stats,_,_=diagnostic(q,r,n,c,'PRELIMINARY truncated ±100k coarse angular average')
    print('QUICK_SAME_N',json.dumps(stats),flush=True)
    for q in np.arange(0,250000.1,25000.):
        if q not in positive:
            r,_=solve_cached(ID,c,N,float(q));positive[q]=r['r_total_y0']
    q=np.arange(-250000.,250000.1,25000.)
    r=np.array([positive[abs(qq)][::-1] if qq<0 else positive[qq] for qq in q])
    stats,rows,curves=diagnostic(q,r,n,c,'Same-N160 Gaussian sensitivity; ±250k, 25k step')
    folder=HERE/'cases'/ID
    np.savez_compressed(folder/'Gaussian_angle_diagnostic_N160.npz',q_rad_m=q,order_n=n,r_total_y0=r,**curves)
    csvwrite(folder/'Gaussian_angle_diagnostic_orders.csv',rows)
    csvwrite(folder/'Gaussian_angle_diagnostic_curves.csv',[{k:float(v[j]) for k,v in curves.items()} for j in range(len(S))])
    stats.update(variant=ID,materials_sha256=sha(ROOT/'fdtd/materials.json'),
        per_order=rows,primary_q0_reference_unchanged=True,
        limitations=['N160 same-truncation diagnostic, not N640-minus-N160.',
            'No physical NA clipping/stop applied; user K already in the D4C coefficient mapping.',
            'Gaussian intensity coefficient average diagnoses finite angular illumination; it is not a newly fitted target.',
            'Negative q reconstructed with exact centered-material mirror symmetry.'])
    dump(HERE/'extra_gaussian_angle_diagnostic.json',stats)
    print('FINAL_SAME_N',json.dumps({k:v for k,v in stats.items() if k!='per_order'}),flush=True)


if __name__=='__main__':main()
