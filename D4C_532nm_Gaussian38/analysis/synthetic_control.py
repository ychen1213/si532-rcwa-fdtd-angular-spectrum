"""Frozen-q0 RCWA + finite Gaussian controls, without any new solver run.

Each structural coefficient comes from saved RCWA. The SAME Gaussian source,
stop widths, leakage and K convention are used for all eight marks. This is a
model-control calculation, never labelled as a genuine finite-mark FDTD result.
"""
from __future__ import annotations
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

import numpy as np
from scipy.special import erfcinv

from zero_order_models import (Z0,T0,KNOWN_ORDERS,K_VALUES,angular_weights,
    b_features,candidate_scan,curve_harmonics,d4c_target,nearest_order,
    scan_from_features,zero_accounting,zero_coefficient,instrument_weights)

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
OUT=HERE/'synthetic_control'
K0=2*np.pi/532e-9
G=2*np.pi/16e-6
YREF=.85e-6
WAIST=19e-6


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def complex_interp(x,a,q):
    return np.interp(q,x,a.real,left=0.,right=0.)+1j*np.interp(q,x,a.imag,left=0.,right=0.)


def measured_source(path):
    """Independent source E/H Fourier projection propagated to ridge plane 0."""
    with np.load(path) as z:
        x=z['x_m'];e=z['Ez'];h=z['Hx'];y=float(z['monitor_y_m'])
    dx=float(np.median(np.diff(x)));n=2**int(np.ceil(np.log2(64*len(x))))
    tw=np.ones(len(x));tw[[0,-1]]=.5
    k=2*np.pi*np.fft.fftshift(np.fft.fftfreq(n,dx))
    phase=np.exp(-1j*k*x[0]);a=dx*np.fft.fftshift(np.fft.fft(tw*e,n))*phase
    ah=dx*np.fft.fftshift(np.fft.fft(tw*h,n))*phase
    sel=np.abs(k)<.98*K0;k=k[sel];a=a[sel];ah=ah[sel]
    gamma=np.sqrt(1-(k/K0)**2)
    down=.5*(a-Z0*ah/gamma)
    source0=down*np.exp(1j*K0*gamma*y)
    pin=float(np.sum(angular_weights(k,1.,NA=1.)*np.abs(down)**2))
    return k,source0,pin,dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),
        Pin_W_m=pin,reference_y_m=0.,definition='air INC E/H downward component, propagated to y0')


def source_grid(kind,k,measured=None):
    if kind=='ideal':
        B=np.sqrt(np.pi)*WAIST*np.exp(-(k*WAIST)**2/4)
        pin=float(np.sum(angular_weights(k,1.,NA=1.)*abs(B)**2))
        return B,pin,lambda q:np.sqrt(np.pi)*WAIST*np.exp(-(q*WAIST)**2/4)
    x,a,pin,_=measured
    B=complex_interp(x,a,k)
    return B,pin,lambda q:complex_interp(x,a,q)


def frozen_spectrum(reference,k,source_at):
    orders=reference['order_n'];r=reference['r_total_y0']
    phase=np.exp(1j*np.sqrt(K0*K0-k*k)*YREF)
    components={int(m):r[j]*source_at(k-m*G)*phase for j,m in enumerate(orders) if abs(m)<=20}
    return sum(components.values()),components


def componentwise_scan(components,k,pin,s,selected,NA=.6):
    """Separately integrate each true +/-m pair, deliberately no cross-m terms.

    Entire positive pupil is used, instead of nearest-order cutting. This is
    the closest finite-template analogue of D4C's independent order sum. NA
    truncation and continuous gamma(k) remain physical finite-spectrum effects.
    """
    w=angular_weights(k,pin,NA=NA);pos=k>0;angle=-4j*np.pi*np.outer(selected,s)/16
    dc=[];cross=[]
    for m in selected:
        plus=components[int(m)];minus_mirror=components[-int(m)][::-1]
        dc.append(float(np.sum(w[pos]*(abs(plus[pos])**2+abs(minus_mirror[pos])**2))))
        cross.append(np.sum(w[pos]*plus[pos]*np.conj(minus_mirror[pos])))
    weights=instrument_weights(selected,convention='K_over_gamma')
    return np.sum(weights[:,None]*(np.array(dc)[:,None]-2*np.real(np.array(cross)[:,None]*np.exp(angle))),axis=0)


def metrics(y,target,s,selected):
    error=np.asarray(y)-target
    hh=curve_harmonics(s,error,selected)
    return dict(max_abs_error=float(max(abs(error))),rmse=float(np.sqrt(np.mean(error**2))),
        relative_l2=float(np.linalg.norm(error)/max(np.linalg.norm(target),1e-30)),
        dc_error=hh['dc'],maximum_cosine_error=float(max(abs(np.asarray(hh['cosine'])))),
        maximum_sine_error=float(max(abs(np.asarray(hh['sine'])))),
        residual_harmonic_reconstruction_error=hh['max_reconstruction_error'])


def main(dk=500.,allow_partial=False):
    OUT.mkdir(exist_ok=True)
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf-8'))
    n=int(np.ceil(.98*K0/dk));k=np.arange(-n,n+1)*dk
    k=k[np.abs(k)<.98*K0]
    s=np.arange(-800,801)*.005
    sourcefile=ROOT/'fdtd/cases/air/nearfield_INC.npz'
    measured=measured_source(sourcefile) if sourcefile.exists() else None
    sources=['ideal']+(['measured'] if measured is not None else [])
    h_T0=float(math.sqrt(2)*erfcinv(T0)/WAIST*1e-6)
    stops={'no_stop':0.,'erfc_matched':h_T0,'old_0p2':.2}
    rows=[];power_rows=[];curve_columns={'s_um':s};files={}
    # A restricted-known-order control isolates the m=18 NA-boundary question.
    selections={'known_all':np.array(KNOWN_ORDERS),'known_le17':np.array([m for m in KNOWN_ORDERS if m<18])}
    for variant in plan['variants']:
        vid=variant['id'];rp=ROOT/'rcwa/cases'/vid/'reference.npz'
        if not rp.exists():
            if allow_partial:continue
            raise FileNotFoundError(f'Need all eight independent RCWA references: {rp}')
        files[vid]=dict(path=rp.relative_to(ROOT).as_posix(),sha256=sha(rp))
        with np.load(rp) as z:reference={name:z[name] for name in z.files}
        rr=reference['r_total_y0'];mm=reference['order_n'];r0=rr[mm==0][0]
        for source in sources:
            B,pin,source_at=source_grid(source,k,measured)
            total,components=frozen_spectrum(reference,k,source_at)
            true_zero=components[0]
            for selection,selected in selections.items():
                plus=rr[np.searchsorted(mm,selected)]
                minus=rr[np.searchsorted(mm,-selected)]
                # The symmetric reference has phi~0; explicit opposite signs
                # are retained in the saved phase audit rather than concealed.
                phi=np.angle(plus*np.conj(minus))
                target=d4c_target(s,selected,abs(plus)**2,abs(r0)**2,phi_m=phi,
                    T_m_definition='User-confirmed |r_+m|^2; I0=|r0|^2')
                comp=componentwise_scan(components,k,pin,s,selected)+zero_coefficient()*float(np.sum(angular_weights(k,pin)*abs(true_zero)**2))
                row=dict(variant=vid,source=source,selection=selection,zero_definition='physical_mode',stop='no_stop',model='componentwise_scalar_zero',**metrics(comp,target,s,selected))
                rows.append(row)
                if selection=='known_all':
                    comp65=componentwise_scan(components,k,pin,s,selected,NA=.65)+zero_coefficient()*float(np.sum(angular_weights(k,pin,NA=.65)*abs(true_zero)**2))
                    rows.append(dict(variant=vid,source=source,selection=selection,zero_definition='physical_mode',stop='no_stop',model='componentwise_scalar_zero_NA065',**metrics(comp65,target,s,selected)))
                    for zeroname,a0 in [('band',None),('physical_mode',true_zero)]:
                        f65=b_features(k,total,pin,selected_orders=selected,NA=.65,
                            stop_halfwidth_rad_um=h_T0,A_zero_ref_V=a0)
                        for model in ('scalar_zero_once','physical_zero_once'):
                            yy=candidate_scan(f65,s,model=model,convention='K_over_gamma')
                            rows.append(dict(variant=vid,source=source,selection=selection,zero_definition=zeroname,
                                stop='erfc_matched',model=model+'_NA065',**metrics(yy,target,s,selected)))
                            if source=='measured':curve_columns[f'{vid}__{zeroname}__{model}_NA065']=yy
                for stopname,width in stops.items():
                    for zeroname,a0 in [('band',None),('physical_mode',true_zero)]:
                        f=b_features(k,total,pin,selected_orders=selected,
                            stop_halfwidth_rad_um=width,A_zero_ref_V=a0)
                        account=zero_accounting(f)
                        power_rows.append(dict(variant=vid,source=source,selection=selection,stop=stopname,
                            halfwidth_rad_um=width,zero_definition=zeroname,
                            raw_power=f.zero_raw_power,passed_power=f.zero_pass_power,
                            blocked_power=f.zero_blocked_power,tau_effective=account['measured_effective_transmission'],
                            physical_zero_once=account['physical_zero_once'],scalar_zero_once=account['scalar_zero_once'],
                            intentionally_repeated_zero=account['double_counted_zero'],
                            cross_power=f.zero_nonzero_cross_pass_power,
                            B0_power=f.zero_band_power,k0_node_power=f.k0_sample_power))
                        # No-stop physical leakage means almost all zero passes;
                        # retain it as an explicit falsification branch.
                        for model in ('scalar_zero_once','physical_zero_once','double_counted_zero'):
                            y=candidate_scan(f,s,model=model,convention='K_over_gamma')
                            rows.append(dict(variant=vid,source=source,selection=selection,zero_definition=zeroname,
                                stop=stopname,model=model,**metrics(y,target,s,selected)))
                            if source=='measured' and selection=='known_all' and stopname=='erfc_matched':
                                curve_columns[f'{vid}__{zeroname}__{model}']=y
                if source=='measured' and selection=='known_all':
                    curve_columns[f'{vid}__target']=target
                    curve_columns[f'{vid}__componentwise']=comp
                    curve_columns[f'{vid}__componentwise_NA065']=comp65
                print(f'Synthetic {vid} {source} {selection} done',flush=True)
    def table(path,data):
        with path.open('w',encoding='utf-8-sig',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    table(OUT/'absolute_errors.csv',rows);table(OUT/'zero_powers.csv',power_rows)
    np.savez_compressed(OUT/'selected_scan_curves.npz',**curve_columns)
    convergence=stop_quadrature_convergence(measured,h_T0)
    zero_only=zero_only_harmonics(k,s,measured,files,h_T0)
    table(OUT/'zero_only_spurious_harmonics.csv',zero_only)
    table(OUT/'stop_quadrature_convergence.csv',convergence)
    summary=dict(kind='Frozen q0 RCWA coefficients times finite incident Gaussian spectrum, not FDTD mark solves',
        normalization='Raw E spectrum and actual source Pin; no mark-specific gain/offset or fitted per-order coefficient',
        source=measured[3] if measured is not None else 'Ideal Gaussian only',reference_files=files,
        all_eight_structures_complete=len(files)==8,k_step_rad_m=dk,k_samples=len(k),
        common_parameters=dict(waist_radius_um=19.,T0=T0,L_SUM=.00138,L_DIFF=.000431,
            erfc_matched_stop_halfwidth_rad_um=h_T0,old_stop_halfwidth_rad_um=.2,K_mapping='K_m/gamma_m'),
        target_mapping='User confirmed T_m=|r_+m|^2; I0=|r0|^2; K2=.4947 and K4=.49135 are user-authorized interpolation',
        measured_source_resampling_power_ratio=float(np.sum(angular_weights(k,measured[2],NA=1.)*abs(complex_interp(measured[0],measured[1],k))**2)) if measured is not None else None,
        stop_quadrature='Centered hard-stop cell overlap; piecewise-constant density in boundary cells, not an optical soft edge',
        stop_quadrature_convergence=convergence,
        zero_only_spurious_harmonics=zero_only,
        NA_boundary_control='Suffix NA065 collects the whole selected m18 lobe; its center sin(theta)=18*.532/16=.5985 lies just inside NA=.6. K18=.248 is retained unchanged. Whether K18 already includes edge loss is NOT established.',
        phase_convention='FDTD B phase negative, target phi=arg(r+ conj r-); symmetric marks have phi~0 so the unresolved sign has no material effect here',
        interpretive_limits=[
            'At m18 the NA=.6 edge cuts part of the finite lobe; K18 may already include an unknown instrument-edge weighting',
            'Band A0 includes some nonzero tails and omits physical zero tails; physical A0 comes from this synthetic construction only',
            'Components of total E are coherent; separate D4C scalar-power channels intentionally omit modal cross power',
            'Physical-stop and scalar-T0 curves can differ even with identical RCWA coefficients; their difference is not an electromagnetic solver error',
            'Synthetic q0 freezing omits angle-dependent RCWA response; the RCWA angle-control study is separate'],
        rows=rows,power_rows=power_rows)
    (OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    make_figures(summary)
    print(f'Wrote {len(rows)} model comparisons for {len(files)} structures',flush=True)
    return summary


def zero_only_harmonics(k,s,measured,files,halfwidth):
    """A zero-only field must have no scan harmonic; B-bin relabelling can add it."""
    rows=[]
    for source in ['ideal']+(['measured'] if measured is not None else []):
        B,pin,_=source_grid(source,k,measured)
        for vid,provenance in files.items():
            with np.load(ROOT/provenance['path']) as ref:
                rr=ref['r_total_y0'][ref['order_n']==0][0]
            a0=rr*B*np.exp(1j*np.sqrt(K0*K0-k*k)*YREF)
            for stop,width in [('no_stop',0.),('erfc_matched',halfwidth)]:
                f=b_features(k,a0,pin,selected_orders=KNOWN_ORDERS,stop_halfwidth_rad_um=width)
                band_curve=scan_from_features(f,s,convention='K_over_gamma')
                physical=b_features(k,a0,pin,selected_orders=KNOWN_ORDERS,
                    stop_halfwidth_rad_um=width,A_zero_ref_V=a0)
                physical_curve=scan_from_features(physical,s,convention='K_over_gamma')
                h=curve_harmonics(s,band_curve,KNOWN_ORDERS)
                rows.append(dict(variant=vid,source=source,stop=stop,
                    band_only_false_scan_range=float(np.ptp(band_curve)),
                    band_only_false_DC=h['dc'],band_only_false_m1_cosine=h['cosine'][0],
                    physical_zero_false_scan_range=float(np.ptp(physical_curve)),
                    true_zero_raw_power=physical.zero_raw_power,
                    false_peak_to_D4C_zero_baseline=float(np.max(band_curve)/(zero_coefficient()*physical.zero_raw_power)),
                    explanation='Only physical zero was supplied. A nonzero scan harmonic is caused solely by nearest-band relabelling of its tails.'))
    return rows


def stop_quadrature_convergence(measured,halfwidth):
    """Separate binary-edge sampling error from the physical Gaussian model."""
    from scipy.integrate import quad
    lower=halfwidth*1e6*WAIST
    f=lambda u:math.sqrt(1-(u/(K0*WAIST))**2)*math.exp(-u*u/2)
    exact_ideal=quad(f,lower,12,epsabs=1e-14)[0]/quad(f,0,12,epsabs=1e-14)[0]
    rows=[]
    for source in ['ideal']+(['measured'] if measured is not None else []):
        for step in (2000.,1000.,500.,250.,125.,62.5):
            n=int(np.floor(.98*K0/step));k=np.arange(-n,n+1)*step
            B,pin,_=source_grid(source,k,measured)
            for cells in (False,True):
                feat=b_features(k,B,pin,selected_orders=[1],A_zero_ref_V=B,
                    stop_halfwidth_rad_um=halfwidth,stop_boundary_cells=cells)
                tau=feat.zero_pass_power/feat.zero_raw_power
                rows.append(dict(source=source,dk_rad_m=step,boundary_cells=cells,
                    measured_tau=tau,ideal_TE_flux_reference=exact_ideal,
                    error_vs_ideal_TE_reference=tau-exact_ideal if source=='ideal' else '',
                    supplied_erfc_T0=T0,error_vs_supplied_erfc=tau-T0,
                    delta_vs_finest_cell_result=0.))
        finest=[r['measured_tau'] for r in rows if r['source']==source and r['dk_rad_m']==62.5 and r['boundary_cells']][0]
        for row in rows:
            if row['source']==source:row['delta_vs_finest_cell_result']=row['measured_tau']-finest
    return rows


def make_figures(summary):
    """Scientific plots of controls, labelled synthetic instead of FDTD marks."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    ids=list(summary['reference_files'])
    labels=[x.replace('ah53_cap','AH').replace('binary50_cap','B') for x in ids]
    source='measured' if isinstance(summary['source'],dict) else 'ideal'
    powers=[r for r in summary['power_rows'] if r['source']==source and r['selection']=='known_all' and r['stop']=='erfc_matched']
    fig,ax=plt.subplots(figsize=(9,4.6),layout='constrained')
    for definition,style,label in [('band','o-','B0 angular band'),('physical_mode','s-','Known physical zero mode')]:
        rr={r['variant']:r for r in powers if r['zero_definition']==definition}
        ax.plot(labels,[rr[x]['tau_effective'] for x in ids],style,label=label)
    ax.axhline(T0,color='black',linestyle='--',label='Supplied scalar T0')
    ax.set(ylabel='Passed / raw zero power',xlabel='Structure (same source and same stop)',
           title='Synthetic RCWA + Gaussian: band and physical zero definitions')
    ax.grid(alpha=.2);ax.legend(loc='upper left')
    fig.savefig(OUT/'zero_transmission_definition.png',dpi=180)
    fig.savefig(OUT/'zero_transmission_definition.svg');plt.close(fig)
    cases=[('band','scalar_zero_once','erfc_matched','B0 + scalar zero'),
           ('physical_mode','scalar_zero_once','erfc_matched','Mode A0 + scalar zero'),
           ('physical_mode','physical_zero_once','erfc_matched','Mode A0 + physical zero'),
           ('physical_mode','componentwise_scalar_zero','no_stop','Independent physical pairs')]
    fig,axes=plt.subplots(2,1,figsize=(9,7.3),sharex=True,layout='constrained')
    x=np.arange(len(ids));width=.19
    for j,(zd,model,stop,label) in enumerate(cases):
        rr={r['variant']:r for r in summary['rows'] if r['source']==source and r['selection']=='known_all' and r['zero_definition']==zd and r['model']==model and r['stop']==stop}
        axes[0].bar(x+(j-1.5)*width,[rr[v]['relative_l2']*100 for v in ids],width,label=label)
        axes[1].bar(x+(j-1.5)*width,[rr[v]['dc_error'] for v in ids],width,label=label)
    axes[0].set(ylabel='Curve relative L2 error (%)',title='Synthetic control against D4C: all 18 orders')
    axes[1].set(ylabel='DC error / incident power',xlabel='Structure',xticks=x,xticklabels=labels)
    axes[0].legend(loc='upper left',ncols=2,fontsize=9)
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.savefig(OUT/'absolute_curve_errors.png',dpi=180)
    fig.savefig(OUT/'absolute_curve_errors.svg');plt.close(fig)
    descriptions=dict(zero_transmission_definition=dict(
        purpose='Test whether one fixed angular stop implements one universal zero transmission across structural variants',
        conclusion='With known physical A0 the effective transmission is universal; defining zero as a nearest angular band introduces structural contamination'),
        absolute_curve_errors=dict(
        purpose='Separate D4C-to-finite-spectrum representation error from electromagnetic solver error',
        conclusion='Finite-spectrum controls have nonzero model residuals even before FDTD; no mark-specific DC correction is applied'))
    (OUT/'figure_purposes.json').write_text(json.dumps(descriptions,ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dk',type=float,default=500.);p.add_argument('--allow-partial',action='store_true')
    args=p.parse_args();main(args.dk,args.allow_partial)
