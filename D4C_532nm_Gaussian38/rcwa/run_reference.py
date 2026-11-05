"""Compute actual matched-material RCWA and the supplied 2026-11-05 D4C model.

The user subsequently confirmed T_m=|r_m|² and authorized K2/K4 linear
interpolation. The complete target is primary; symbolic K sensitivity remains.
"""
from __future__ import annotations
from pathlib import Path
import csv,hashlib,json,math,time
from dataclasses import replace
import numpy as np
from te_planarized_rcwa import Config,solve_angle,dump

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'rcwa'
K_VALUES={1:.4957,3:.4937,5:.4890,6:.4860,7:.4875,8:.4783,9:.4704,
    10:.4708,11:.4629,12:.4540,13:.4478,14:.4380,15:.4293,16:.4200,17:.4056,18:.2480}
K_PROVISIONAL={2:.4947,4:.49135}
K_CONFIRMED={3,6,7,9,10,11,12,14,15}
T0=math.erfc((5e-4*math.sqrt(2))/(.778*3.76e-4))
L={'SUM':.00138,'DIFF':.000431}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csvwrite(path,rows):
    with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
        columns=list(dict.fromkeys(k for row in rows for k in row))
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)


def make_config(v,materials):
    def index(name):
        z=materials['materials'][name]['FDTD_fit_index'];return complex(z['n'],z['k'])
    si=index('Si');oxide=index('SiO2');cap=v['cap_nm']*1e-9
    return Config(wavelength_m=532e-9,period_m=v['period_um']*1e-6,
        height_m=v['ridge_height_nm']*1e-9,ridge_width_m=v['ridge_width_um']*1e-6,
        ridge_centers_m=tuple(x*1e-6 for x in v['ridge_centers_local_um']),
        cap_m=cap,ridge_n=si,substrate_n=si,fill_n=oxide if cap else 1+0j,cap_n=oxide)


def solve_cached(variant,c,N,q=0.):
    folder=HERE/'cases'/variant;folder.mkdir(parents=True,exist_ok=True)
    qname='q0' if not q else f'q{int(round(q)):+d}'
    path=folder/f'N{N}_{qname}.npz';meta=path.with_suffix('.json')
    fingerprint=hashlib.sha256(json.dumps(c.metadata(),sort_keys=True).encode()).hexdigest()
    if path.exists() and meta.exists():
        m=json.loads(meta.read_text(encoding='utf8'))
        if m.get('config_sha256')==fingerprint and m.get('solver_sha256')==sha(HERE/'te_planarized_rcwa.py'):
            with np.load(path) as z:return {k:z[k] for k in z.files},m
    started=time.perf_counter();r=solve_angle(q,N,c)
    np.savez_compressed(path,**r)
    m=dict(variant=variant,config=c.metadata(),config_sha256=fingerprint,
        solver_sha256=sha(HERE/'te_planarized_rcwa.py'),max_order=N,q_rad_m=q,
        elapsed_s=time.perf_counter()-started,R=r['R'],T_interface=r['T_interface'],
        finite_layer_absorption=r['finite_layer_absorption'],npz=path.relative_to(ROOT).as_posix())
    dump(meta,m);print(f"RCWA {variant} N{N} {qname} {m['elapsed_s']:.2f}s R={r['R']:.9f}",flush=True)
    return r,m


def compare(new,old):
    order=np.arange(-18,19)
    def get(r):return r['r_total_y0'][np.searchsorted(r['order_n'],order)]
    a=get(new);b=get(old);delta=a-b
    return dict(max_complex_abs_change=float(abs(delta).max()),
        relative_complex_L2=float(np.linalg.norm(delta)/np.linalg.norm(a)),
        zero_intensity_abs_change=float(abs(abs(a[18])**2-abs(b[18])**2)),
        nonzero_relative_complex_L2=float(np.linalg.norm(np.delete(delta,18))/np.linalg.norm(np.delete(a,18))))


def main():
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf8'))
    materials=json.loads((ROOT/'fdtd/materials.json').read_text(encoding='utf8'))
    s=np.arange(-800,801)*.005;mpos=np.arange(1,19)
    all_orders=[];conv=[];summaries=[];finals={};selected=[]
    source_sha=dict(experiment=sha(ROOT/'experiment.json'),materials=sha(ROOT/'fdtd/materials.json'),
        user_formula=sha(ROOT/'inputs/D4C_FDTD_user_formula_20261105.md'))
    for v in plan['variants']:
        c=make_config(v,materials);previous=None;history=[]
        for N in (160,240,320,480,640):
            r,meta=solve_cached(v['id'],c,N)
            item=dict(variant=v['id'],max_order=N,full_harmonic_count=2*N+1,**{k:meta[k] for k in ('elapsed_s','R','T_interface','finite_layer_absorption')})
            if previous is not None:item.update(compare(r,previous))
            conv.append(item);history.append(item)
            # Require at least N320, then stop when the collected complex vector is stable to 0.05%.
            if N>=320 and item['relative_complex_L2']<.0005 and item['nonzero_relative_complex_L2']<.001:
                break
            previous=r
        finals[v['id']]=(r,c);selected.append(dict(variant=v['id'],max_order=N,npz=meta['npz']))
        order=r['order_n'];rr=r['r_total_y0'];j0=int(np.flatnonzero(order==0)[0]);I0=float(abs(rr[j0])**2)
        pp=np.searchsorted(order,mpos);nn=np.searchsorted(order,-mpos)
        rp=rr[pp];rn=rr[nn];Tm=abs(rp)**2;phase=np.angle(rp*np.conj(rn))
        cs=2*Tm[:,None]*(1-np.cos(2*mpos[:,None]*c.G*s[None,:]*1e-6+phase[:,None]))
        known=sum(K_VALUES[int(m)]*cs[j] for j,m in enumerate(mpos) if m in K_VALUES)
        missing2=cs[1];missing4=cs[3]
        interp=known+K_PROVISIONAL[2]*missing2+K_PROVISIONAL[4]*missing4
        zsum=(T0+L['SUM']*(1-T0))*I0;zdiff=(T0+L['DIFF']*(1-T0))*I0
        folder=HERE/'cases'/v['id']
        np.savez_compressed(folder/'d4c_targets.npz',s_um=s,order_positive=mpos,r_plus=rp,r_minus=rn,
            T_m=Tm,T_m_working_definition=Tm,phi_m_rad=phase,I0=I0,T0=T0,
            AC=interp,S_SUM=interp+zsum,S_DIFF_same_kernel=interp+zdiff,
            AC_known_K=known,coefficient_K2=missing2,coefficient_K4=missing4,
            AC_provisional_interpolated_K2K4=interp,
            zero_SUM=zsum,zero_DIFF=zdiff,
            S_SUM_known_K=known+zsum,S_DIFF_same_kernel_known_K=known+zdiff,
            S_SUM_provisional_interpolated_K2K4=interp+zsum,
            S_DIFF_same_kernel_provisional_interpolated_K2K4=interp+zdiff)
        rows=[dict(s_um=float(ss),AC=float(p),S_SUM=float(p+zsum),S_DIFF_same_kernel=float(p+zdiff),
            AC_known_K=float(a),coefficient_K2=float(b),coefficient_K4=float(d),
            zero_SUM=zsum,zero_DIFF=zdiff,S_SUM_known_K=float(a+zsum),S_DIFF_same_kernel_known_K=float(a+zdiff),
            S_SUM_provisional_interpolated_K2K4=float(p+zsum),S_DIFF_same_kernel_provisional_interpolated_K2K4=float(p+zdiff))
            for ss,a,b,d,p in zip(s,known,missing2,missing4,interp)]
        csvwrite(folder/'d4c_targets.csv',rows)
        np.savez_compressed(folder/'reference.npz',**r,r_common_y085=r['r_total_y0']*np.exp(1j*np.sqrt(np.maximum(0,c.k0**2-r['kx_rad_m']**2))*.85e-6),
            config_json=json.dumps(c.metadata(),sort_keys=True),material_source_sha256=source_sha['materials'])
        for j,m in enumerate(order):
            z=rr[j];all_orders.append(dict(variant=v['id'],max_order=N,m=int(m),kx_rad_um=float(r['kx_rad_m'][j]*1e-6),
                propagating=bool(abs(r['kx_rad_m'][j])<c.k0),r_y0_re=float(z.real),r_y0_im=float(z.imag),r_abs=float(abs(z)),
                r_phase_rad=float(np.angle(z)),intensity_abs_r_squared=float(abs(z)**2),R_order_fraction=float(r['reflected_order_power_fraction'][j])))
        summaries.append(dict(variant=v['id'],geometry=v['geometry'],cap_nm=v['cap_nm'],max_order=N,
            I0=I0,r0_re=float(rr[j0].real),r0_im=float(rr[j0].imag),T0=T0,
            residual_zero_per_I0=T0,blocked_leak_SUM_per_I0=L['SUM']*(1-T0),blocked_leak_DIFF_per_I0=L['DIFF']*(1-T0),
            zero_SUM=zsum,zero_DIFF=zdiff,R=float(r['R']),T_interface=float(r['T_interface']),
            finite_layer_absorption=float(r['finite_layer_absorption']),
            mirror_max_complex_error=float(np.max(abs(rr-rr[::-1]))),
            max_positive_negative_pair_phase_rad=float(abs(phase).max()),
            AC_period_mean=float(interp[:-1].mean()),AC_provisional_period_mean=float(interp[:-1].mean()),
            missing_K2_coefficient_period_mean=float(missing2[:-1].mean()),
            missing_K4_coefficient_period_mean=float(missing4[:-1].mean())))
    csvwrite(HERE/'all_orders.csv',all_orders);csvwrite(HERE/'convergence.csv',conv)
    csvwrite(HERE/'structure_summary.csv',summaries)
    csvwrite(HERE/'K_table.csv',[dict(m=m,K=K_VALUES.get(m,K_PROVISIONAL.get(m)),status='user_authorized_linear_interpolation' if m in (2,4) else 'user_corrected' if m in K_CONFIRMED else 'prior_reading',
        sensitivity_interpolation=K_PROVISIONAL.get(m,'')) for m in range(1,19)])
    dump(HERE/'reference_summary.json',dict(T0=T0,L=L,zero_coefficient_per_I0={k:T0+x*(1-T0) for k,x in L.items()},
        selected=selected,structures=summaries,source_sha256=source_sha,
        working_definition='T_m=abs(r_+m)^2 confirmed by user 2026-11-05; I0=abs(r0)^2; phi=arg(r_+m conj(r_-m)).',
        formula='2 sum_m K_m T_m [1-cos(2mGs+phi_m)] + [T0+L_chi(1-T0)] I0',
        K_interpolation_policy='User authorized linear interpolation: K2=.4947 and K4=.49135. Primary S_SUM/S_DIFF_same_kernel/AC contain all m=1..18. Symbolic coefficients remain for sensitivity. Legacy provisional-named arrays are identical compatibility aliases, not a second target.',
        channel_scope='DIFF uses the same supplied dark interference kernel solely to compare L_DIFF; actual DIFF optical sign remains unspecified.',
        phase_convention='exp(-i omega t), Ez; air incident -y, reflection +y; incident amplitude at y0=1; lateral origin x0 at symmetric center.',
        field_model_scope='Infinite periodic RCWA at q0. Finite nine-period FDTD is compared separately. Gaussian frozen-q0 template is an approximation, not a finite FDTD solve.'))
    print('RCWA all 8 structures and D4C symbolic/interpolated targets complete.',flush=True)


if __name__=='__main__':main()
