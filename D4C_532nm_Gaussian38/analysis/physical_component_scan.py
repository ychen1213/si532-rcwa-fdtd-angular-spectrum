"""Independent source-projected physical-wavepacket B scans from real FDTD.

No RCWA amplitudes enter the field decomposition. Saved D4C target curves are
read ONLY for comparison and shared training/holdout parameter estimation.
The polynomial source projection is an approximate representation, not a unique
physical modal measurement; its field residual is retained in the audit.
"""
from __future__ import annotations
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
from scipy.optimize import least_squares

from zero_order_models import (K_VALUES,T0,angular_weights,b_features,
    evaluate_global_parameters,fit_global_parameters,instrument_weights,
    scan_from_features,zero_coefficient)

ROOT=Path(__file__).resolve().parents[1]
K0=2*np.pi/532e-9
G=2*np.pi/16e-6
SELECTED=np.arange(1,19)


def dump(path,obj):Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf8')


def table(path,rows):
    if not rows:return
    with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def interpolate(k,a,q):
    return np.interp(q,k,a.real,left=0,right=0)+1j*np.interp(q,k,a.imag,left=0,right=0)


def reconstruct_components(data,metadata,degree=2):
    """Reconstruct the SAME shifted-source polynomial basis saved by root."""
    k=data['kx_rad_m'];source=data['Ainc_y0_V'];projection=metadata[f'projection_degree{degree}']
    orders=np.asarray(projection['order'],int)
    coeff=np.asarray(projection['fitted_coefficients_re'])+1j*np.asarray(projection['fitted_coefficients_im'])
    if coeff.shape!=(len(orders),degree+1):raise ValueError('Saved coefficient shape differs from degree/order schema')
    yref=float(data['reference_y_m']);gamma=np.sqrt(1-(k/K0)**2)
    scale=float(projection['incident_q_scale_rad_m']);cutoff=float(projection['template_q_cutoff_rad_m'])
    components={}
    for j,m in enumerate(orders):
        q=k-m*G
        base=interpolate(k,source,q)*(abs(q)<cutoff)*np.exp(1j*K0*gamma*yref)
        response=sum(coeff[j,d]*(q/scale)**d for d in range(degree+1))
        components[int(m)]=base*response
    return components


def component_features(k,components,pin,NA=.6):
    """Integrate complete +m / reflected -m packets over positive k once.

    The saved |q|<G source cutoff guarantees positive m>=1 packet support does
    not extend below k=0 (m=1 starts at 0), within THIS representation. Real
    Gaussian and Maxwell fields have tails beyond that artificial cutoff.
    Those tails are a small-model approximation, never asserted to vanish.
    """
    w=angular_weights(k,pin,NA=NA);positive=k>0
    dc=[];cross=[];power_rows=[]
    for m in SELECTED:
        plus=components[int(m)];minus_reflected=components[-int(m)][::-1]
        pp=float(np.sum(w[positive]*abs(plus[positive])**2))
        pm=float(np.sum(w[positive]*abs(minus_reflected[positive])**2))
        cc=np.sum(w[positive]*plus[positive]*np.conj(minus_reflected[positive]))
        dc.append(pp+pm);cross.append(cc)
        outside=float(np.sum(w[~positive]*abs(plus[~positive])**2))
        power_rows.append(dict(m=int(m),P_plus_over_Pin=pp,P_minus_over_Pin=pm,
            cross_re=float(cc.real),cross_im=float(cc.imag),
            plus_packet_power_nonpositive_k=outside,
            visibility=2*float(abs(cc))/max(pp+pm,1e-30)))
    return np.asarray(dc),np.asarray(cross),power_rows


def packet_scan(dc,cross,s_um):
    weights=instrument_weights(SELECTED,convention='K_over_gamma',K=K_VALUES)
    phase=np.exp(-4j*np.pi*np.outer(SELECTED,np.asarray(s_um))/16.)
    return np.sum(weights[:,None]*(dc[:,None]-2*np.real(cross[:,None]*phase)),axis=0)


def metrics(target,y,s0_index=None):
    t=np.asarray(target)[:-1];v=np.asarray(y)[:-1];e=v-t
    centered=(v-v.mean())-(t-t.mean())
    j=len(target)//2 if s0_index is None else int(s0_index)
    baseline_target=float(np.asarray(target)[j]);baseline_pred=float(np.asarray(y)[j])
    return dict(relative_L2_pct=float(100*np.linalg.norm(e)/max(np.linalg.norm(t),1e-30)),
        centered_relative_L2_pct=float(100*np.linalg.norm(centered)/max(np.linalg.norm(t-t.mean()),1e-30)),
        RMSE=float(np.sqrt(np.mean(e*e))),max_abs_error=float(max(abs(e))),
        DC_error=float(e.mean()),span_ratio=float(np.ptp(v)/max(np.ptp(t),1e-30)),
        baseline_target=baseline_target,baseline_pred=baseline_pred,
        baseline_relative_pct=100*(baseline_pred/baseline_target-1) if baseline_target else None)


def fit_nonzero_gain_parameters(records,*,free=('gain_ac','zero_transmission'),fixed=None,baseline_aware=False):
    """Shared gain on the nonzero scan only; scalar zero is not gain-scaled.

    This is a separate low-dimensional readout candidate, not a physical
    calibration claim. All records share one gain_ac and one T0. Leakage values
    are fixed to the supplied constants; no mark-specific baseline is fitted.
    Records here have already excluded the duplicate periodic endpoint.
    """
    defaults=dict(gain_ac=1.,zero_transmission=T0,leak_sum=.00138,leak_diff=.000431)
    if fixed:defaults.update(fixed)
    if not set(free)<=set(('gain_ac','zero_transmission')):raise ValueError('Only shared gain_ac/T0 may be fitted')
    if not records:raise ValueError('No records for shared nonzero-gain fit')
    bounds={'gain_ac':(.5,1.5),'zero_transmission':(0.,.02)}
    def unpack(x):return {**defaults,**dict(zip(free,map(float,x)))}
    def predict(r,p):
        leak=p['leak_'+r['channel'].lower()]
        c=zero_coefficient(r['channel'],p['zero_transmission'],leak)
        return p['gain_ac']*np.asarray(r['nonzero'])+c*float(r['zero_power'])
    scales=[max(float(np.sqrt(np.mean(np.asarray(r['target'])**2))),1e-12) for r in records]
    def residual(x):
        p=unpack(x)
        blocks=[]
        for r,scale in zip(records,scales):
            error=predict(r,p)-r['target']
            blocks.append(error/(scale*np.sqrt(len(r['target']))))
            if baseline_aware:
                j=int(r['s0_index']);baseline=float(r['target'][j])
                if baseline<=0:raise ValueError('Relative s0 objective requires positive D4C baseline')
                blocks.append(np.array([error[j]/baseline]))
        return np.concatenate(blocks)
    def jacobian(x):
        p=unpack(x);blocks=[]
        for r,scale in zip(records,scales):
            nz=np.asarray(r['nonzero']);leak=p['leak_'+r['channel'].lower()]
            derivatives=dict(gain_ac=nz,zero_transmission=np.full_like(nz,(1-leak)*float(r['zero_power'])))
            blocks.append(np.column_stack([derivatives[name] for name in free])/(scale*np.sqrt(len(nz))))
            if baseline_aware:
                j=int(r['s0_index']);baseline=float(r['target'][j])
                blocks.append(np.array([[derivatives[name][j]/baseline for name in free]]))
        return np.vstack(blocks)
    if free:
        result=least_squares(residual,[defaults[name]for name in free],jac=jacobian,
            bounds=([bounds[name][0]for name in free],[bounds[name][1]for name in free]),
            xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=500)
        params=unpack(result.x);singular=np.linalg.svd(result.jac,compute_uv=False)
        rank=int(sum(singular>singular[0]*1e-7))
        if rank<len(free):raise ValueError(f'Nonzero-gain parameters not identifiable: rank {rank}/{len(free)}')
        success=bool(result.success)
    else:params=unpack([]);singular=np.array([]);rank=0;success=True
    rows=[]
    for r in records:
        predicted=predict(r,params);error=predicted-r['target'];j=int(r.get('s0_index',len(error)//2))
        baseline=float(r['target'][j]);bp=float(predicted[j])
        rows.append(dict(id=r['id'],channel=r['channel'],
            rmse=float(np.sqrt(np.mean(error**2))),max_abs_error=float(max(abs(error))),
            dc_error=float(error.mean()),relative_l2=float(np.linalg.norm(error)/max(np.linalg.norm(r['target']),1e-30)),
            baseline_target=baseline,baseline_pred=bp,baseline_relative_pct=100*(bp/baseline-1) if baseline else None))
    return dict(parameters=params,free=list(free),success=success,gain_scope='nonzero_only',
        formula='gain_ac*nonzero + [T0+L_channel*(1-T0)]*zero_power',
        objective='curve_RMS_normalized_and_s0_relative_equal_weight' if baseline_aware else 'curve_RMS_normalized_only',
        s0_objective_weight=1. if baseline_aware else 0.,
        records=rows,jacobian_rank=rank,jacobian_singular_values=singular.tolist(),
        shared_across_marks=True,per_mark_offsets=False,leakage_fixed=True)


def audit_nonzero_gain_fit():
    """Synthetic algebra check of the shared-parameter API, not real validation."""
    s=np.arange(1600)*.005;truth=dict(gain_ac=1.013,zero_transmission=.00041)
    records=[]
    for j in range(8):
        nz=(.007+.003*j)*(1-np.cos(4*np.pi*s/16))+.002*(1-np.cos(12*np.pi*s/16))
        p0=.025+.037*j
        for channel in ('SUM','DIFF'):
            target=truth['gain_ac']*nz+zero_coefficient(channel,truth['zero_transmission'])*p0
            records.append(dict(id=f'synthetic_{j}',channel=channel,target=target,nonzero=nz,zero_power=p0,s0_index=0))
    training=[r for r in records if int(r['id'].split('_')[-1])<4]
    heldout=[r for r in records if int(r['id'].split('_')[-1])>=4]
    fit=fit_nonzero_gain_parameters(training)
    held=fit_nonzero_gain_parameters(heldout,free=(),fixed=fit['parameters'])
    eg=abs(fit['parameters']['gain_ac']-truth['gain_ac'])
    et=abs(fit['parameters']['zero_transmission']-truth['zero_transmission'])
    error=max(r['max_abs_error']for r in held['records'])
    assert eg<1e-10 and et<1e-10 and error<1e-12
    aware=fit_nonzero_gain_parameters(training,baseline_aware=True)
    held_aware=fit_nonzero_gain_parameters(heldout,free=(),fixed=aware['parameters'],baseline_aware=True)
    aware_error=max(r['max_abs_error']for r in held_aware['records'])
    assert aware_error<1e-12
    return dict(synthetic_only=True,training_structures=4,heldout_structures=4,
        known_gain_ac=truth['gain_ac'],known_zero_transmission=truth['zero_transmission'],
        recovered_parameters=fit['parameters'],gain_error=eg,transmission_error=et,
        heldout_max_absolute_error=error,jacobian_rank=fit['jacobian_rank'],
        baseline_aware_recovered_parameters=aware['parameters'],baseline_aware_heldout_max_absolute_error=aware_error,
        interpretation='Checks shared nonzero-gain algebra only; does not validate any real FDTD or instrument model')


def add_whole_fit_baseline_diagnostics(fit,records):
    """Add s0 diagnostics to the unchanged generic whole-gain fit output."""
    indexed={(r['id'],r['channel']):r for r in records};p=fit['parameters']
    for row in fit['records']:
        r=indexed[(row['id'],row['channel'])];j=int(r['s0_index'])
        baseline=float(r['target'][j])
        predicted=p['gain']*(float(r['nonzero'][j])+zero_coefficient(r['channel'],p['zero_transmission'])*float(r['zero_power']))
        row.update(baseline_target=baseline,baseline_pred=predicted,
            baseline_relative_pct=100*(predicted/baseline-1) if baseline else None)
    return fit


def field_diagnostics(k,data,components,NA):
    a=data['A_ref_V'];pred=sum(components.values());zero=components[0]
    savedzero=data['A_zero_degree2_V'];w=angular_weights(k,float(data['Pin_W_m']),NA=NA)
    norm=lambda v:float(np.sum(w*abs(v)**2))
    nonzero=a-zero;pn=pred-zero
    return dict(NA=NA,projection_total_relative_field_L2=float(np.sqrt(norm(a-pred)/max(norm(a),1e-30))),
        projection_nonzero_relative_field_L2=float(np.sqrt(norm(nonzero-pn)/max(norm(nonzero),1e-30))),
        reconstructed_A0_vs_saved_relative_field_L2=float(np.sqrt(norm(zero-savedzero)/max(norm(savedzero),1e-30))),
        actual_total_power=norm(a),projection_total_power=norm(pred),
        physical_zero_power=norm(zero),projection_residual_power=norm(a-pred),
        residual_coherent_cross_power=float(2*np.real(np.sum(w*pred*np.conj(a-pred)))))


def analyze(suffix='_fine',require_all_eight=False):
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf-8-sig'))
    out=ROOT/'outputs'/('physical_component'+(suffix or '_coarse'));out.mkdir(parents=True,exist_ok=True)
    rows=[];powers=[];projection=[];inputs=[];curves={};records_by_method={};available=[]
    for variant in plan['variants']:
        vid=variant['id'];path=ROOT/'outputs/spectra'/(vid+suffix+'.npz');meta_path=path.with_suffix('.json')
        if not path.exists() or not meta_path.exists():continue
        with np.load(path) as z:d={name:z[name] for name in z.files}
        meta=json.loads(meta_path.read_text(encoding='utf-8-sig'))
        targetpath=ROOT/'rcwa/cases'/vid/'d4c_targets.npz'
        with np.load(targetpath) as z:
            s=z['s_um'];targets={'SUM':z['S_SUM'],'DIFF':z['S_DIFF_same_kernel']}
        components=reconstruct_components(d,meta)
        pred=sum(components.values());zero=components[0];pin=float(d['Pin_W_m']);k=d['kx_rad_m']
        available.append(vid)
        inputs.append(dict(id=vid,spectrum_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            projection_metadata_sha256=hashlib.sha256(meta_path.read_bytes()).hexdigest(),
            target_curve_sha256=hashlib.sha256(targetpath.read_bytes()).hexdigest()))
        for NA in (.6,.65):
            tag=f'NA{NA:.2f}'
            projection.append(dict(id=vid,**field_diagnostics(k,d,components,NA)))
            dc,cross,pr=component_features(k,components,pin,NA)
            powers.extend(dict(id=vid,NA=NA,**r) for r in pr)
            packet=packet_scan(dc,cross,s)
            actual_f=b_features(k,d['A_ref_V'],pin,selected_orders=SELECTED,A_zero_ref_V=zero,NA=NA)
            projected_f=b_features(k,pred,pin,selected_orders=SELECTED,A_zero_ref_V=zero,NA=NA)
            nonzero={f'actual_nearest_band_{tag}':scan_from_features(actual_f,s,convention='K_over_gamma'),
                f'projected_nearest_band_{tag}':scan_from_features(projected_f,s,convention='K_over_gamma'),
                f'projected_physical_pairs_{tag}':packet}
            for method,nz in nonzero.items():
                records_by_method.setdefault(method,[])
                for channel in ('SUM','DIFF'):
                    y=nz+zero_coefficient(channel)*actual_f.zero_raw_power
                    rows.append(dict(id=vid,cap_nm=variant['cap_nm'],geometry=variant['geometry'],
                        channel=channel,method=method,**metrics(targets[channel],y)))
                    curves[f'{vid}__{channel}__{method}']=y
                    records_by_method[method].append(dict(id=vid,channel=channel,target=targets[channel],
                        nonzero=nz,zero_power=actual_f.zero_raw_power,cap_nm=variant['cap_nm'],s0_index=int(np.argmin(abs(s)))))
        curves[f'{vid}__SUM__D4C']=targets['SUM'];curves[f'{vid}__DIFF__D4C']=targets['DIFF']
        print(f'Physical packet control: {vid} complete',flush=True)
    if require_all_eight and len(available)!=8:raise RuntimeError(f'Only {len(available)} real spectra available; eight required')
    if not available:raise RuntimeError('No processed real FDTD spectra available')
    fits={};fitrows=[]
    for method,records in records_by_method.items():
        # Duplicated endpoints are excluded in both fitting and score metrics.
        train=[dict(r,target=r['target'][:-1],nonzero=r['nonzero'][:-1]) for r in records if r['cap_nm'] in (0,100)]
        held=[dict(r,target=r['target'][:-1],nonzero=r['nonzero'][:-1]) for r in records if r['cap_nm'] in (150,200)]
        if not train:continue
        for free in (('gain',),('zero_transmission',),('gain','zero_transmission')):
            key=method+'__'+'+'.join(free)
            fit=fit_global_parameters(train,free=free,bounds={'gain':(.5,1.5),'zero_transmission':(0.,.02)})
            add_whole_fit_baseline_diagnostics(fit,train)
            fit['gain_scope']='whole_signal'
            fit['objective']='curve_RMS_normalized_only'
            fit['heldout']=add_whole_fit_baseline_diagnostics(evaluate_global_parameters(held,fit['parameters']),held) if held else None
            fits[key]=fit;p=fit['parameters']
            for r in records:
                y=p['gain']*(r['nonzero']+zero_coefficient(r['channel'],p['zero_transmission'])*r['zero_power'])
                fitrows.append(dict(id=r['id'],channel=r['channel'],method=method,free='+'.join(free),gain_scope='whole_signal',objective='curve_only',
                    split='train' if r['cap_nm'] in (0,100) else 'heldout',**metrics(r['target'],y)))
                curves[f"{r['id']}__{r['channel']}__fit_{key}"]=y
        for free in (('gain_ac',),('zero_transmission',),('gain_ac','zero_transmission')):
            key=method+'__nonzero_only__'+'+'.join(free)
            fit=fit_nonzero_gain_parameters(train,free=free)
            fit['heldout']=fit_nonzero_gain_parameters(held,free=(),fixed=fit['parameters']) if held else None
            fits[key]=fit;p=fit['parameters']
            for r in records:
                y=p['gain_ac']*r['nonzero']+zero_coefficient(r['channel'],p['zero_transmission'])*r['zero_power']
                fitrows.append(dict(id=r['id'],channel=r['channel'],method=method,free='+'.join(free),gain_scope='nonzero_only',objective='curve_only',
                    split='train' if r['cap_nm'] in (0,100) else 'heldout',**metrics(r['target'],y)))
                curves[f"{r['id']}__{r['channel']}__fit_{key}"]=y
        free=('gain_ac','zero_transmission')
        key=method+'__nonzero_only_baseline_aware__'+'+'.join(free)
        fit=fit_nonzero_gain_parameters(train,free=free,baseline_aware=True)
        fit['heldout']=fit_nonzero_gain_parameters(held,free=(),fixed=fit['parameters'],baseline_aware=True) if held else None
        fits[key]=fit;p=fit['parameters']
        for r in records:
            y=p['gain_ac']*r['nonzero']+zero_coefficient(r['channel'],p['zero_transmission'])*r['zero_power']
            fitrows.append(dict(id=r['id'],channel=r['channel'],method=method,free='+'.join(free),gain_scope='nonzero_only',objective='curve_and_s0_equal_weight',
                split='train' if r['cap_nm'] in (0,100) else 'heldout',**metrics(r['target'],y)))
            curves[f"{r['id']}__{r['channel']}__fit_{key}"]=y
    table(out/'candidate_errors.csv',rows);table(out/'component_powers.csv',powers)
    table(out/'projection_audit.csv',projection);table(out/'shared_fit_errors.csv',fitrows)
    dump(out/'shared_fits.json',fits)
    main_keys=('actual_nearest_band_NA0.60__nonzero_only__gain_ac',
        'actual_nearest_band_NA0.60__nonzero_only_baseline_aware__gain_ac+zero_transmission')
    dump(out/'main_branch_readout_comparison.json',{key:dict(parameters=fits[key]['parameters'],
        heldout=fits[key]['heldout']['records'] if fits[key]['heldout']else None,
        objective=fits[key]['objective'])for key in main_keys if key in fits})
    fit_audit=audit_nonzero_gain_fit()
    dump(out/'nonzero_gain_math_audit.json',fit_audit)
    np.savez_compressed(out/'scan_curves.npz',s_um=s,**curves)
    summary=dict(available_structures=available,all_eight_complete=len(available)==8,inputs=inputs,
        training_structures=[v['id']for v in plan['variants']if v['id']in available and v['cap_nm']in(0,100)],
        heldout_structures=[v['id']for v in plan['variants']if v['id']in available and v['cap_nm']in(150,200)],
        reconstruction='Saved real-FDTD measured-source joint degree2 polynomial coefficients; no RCWA amplitudes used',
        readout='Each physical +/-m source packet gets a single exp(-imGs); independent-pair power sums; scalar zero exactly once',
        target_curve_role='Saved D4C curves only used after reconstruction for comparison and shared parameter training',
        selected_orders=SELECTED.tolist(),NA_controls=[.6,.65],instrument_weight='K_m/gamma_m with user-authorized K2/K4 interpolation',
        zero_term='[T0+L_channel*(1-T0)] times full-pupil reconstructed physical-zero power',
        training='both geometries cap0/cap100 when available',heldout='both geometries cap150/cap200 when available',
        fitting_candidates=['whole_signal: gain*(nonzero+c(T0)*P0)','nonzero_only: gain_ac*nonzero+c(T0)*P0; fixed leakage'],
        fitting_objectives=['curve-only, retained including failures','nonzero-gain joint fit with one equal-weight relative s0 block per curve'],
        baseline_protection='The baseline-aware branch jointly minimizes curve/(RMS*sqrt(N)) and (s0_pred-s0_target)/s0_target, weight1 for both blocks for every training record',
        nonzero_gain_math_audit=fit_audit,
        partial_run_warning='No eight-mark cross-validation claim unless all_eight_complete is true',
        limitations=[
            'Source projection is an approximate polynomial field representation, not unique Maxwell mode decomposition',
            'Independent physical-pair summation deliberately discards inter-order coherent cross terms, as the supplied D4C order sum does',
            'Projection residual is omitted by packet scans; projected-nearest-band control separately quantifies its removal',
            'The saved |q|<G source cutoff makes the modeled +m packets have no nonpositive-k power for m>=1; actual far tails need not vanish',
            'NA=.65 is a representation control for full selected18 packets, not a claim that the real optical NA differs from .6',
            'No per-mark DC offset and no per-order fitted gain is allowed',
            'SUM/DIFF here differ only in supplied leakage, with no assertion about unresolved hardware interference signs'],
        rerun=f'python analysis/physical_component_scan.py --suffix {suffix} --require-all-eight')
    dump(out/'summary.json',summary)
    print(f'Wrote controls for {len(available)} real FDTD structures to {out}',flush=True)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--suffix',default='_fine');parser.add_argument('--require-all-eight',action='store_true')
    args=parser.parse_args();analyze(args.suffix,args.require_all_eight)
