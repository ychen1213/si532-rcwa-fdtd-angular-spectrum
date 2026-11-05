"""Compare actual D4C target curves with independently solved FDTD B curves.

One set of readout parameters is fitted on caps 0/100 nm; 150/200 nm are
held out. No structure-specific baseline or per-order fitted equalizer.
"""
from __future__ import annotations
from pathlib import Path
import argparse,csv,json
import numpy as np
from zero_order_models import (b_features,candidate_scan,scan_from_features,
    fit_global_parameters,evaluate_global_parameters,K_VALUES,T0,zero_coefficient)

ROOT=Path(__file__).resolve().parents[1]

def dump(path,obj):Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf8')
def writecsv(path,rows):
    if not rows:return
    with Path(path).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def periodic_records(records):
    """Do not give the duplicated periodic endpoint extra fitting weight."""
    return [dict(r,target=np.asarray(r['target'])[:-1],nonzero=np.asarray(r['nonzero'])[:-1]) for r in records]

def metrics(target,y):
    # Exclude duplicated endpoint for periodic mean/harmonic metrics.
    t=np.asarray(target)[:-1];a=np.asarray(y)[:-1];e=a-t
    scale=max(float(np.linalg.norm(t)),1e-30);span=max(float(np.ptp(t)),1e-30)
    tc=t-t.mean();ac=a-a.mean()
    center=len(np.asarray(target))//2
    zt=float(np.asarray(target)[center]);za=float(np.asarray(y)[center])
    return dict(relative_L2_pct=float(100*np.linalg.norm(e)/scale),
        RMSE=float(np.sqrt(np.mean(e**2))),max_abs_error=float(np.max(np.abs(e))),
        mean_error=float(e.mean()),mean_error_over_target_span_pct=float(100*e.mean()/span),
        centered_relative_L2_pct=float(100*np.linalg.norm(ac-tc)/max(np.linalg.norm(tc),1e-30)),
        span_ratio=float(np.ptp(a)/span),minimum_target=float(t.min()),minimum_fdtd=float(a.min()),
        target_mean=float(t.mean()),fdtd_mean=float(a.mean()),
        baseline_s0_target=zt,baseline_s0_fdtd=za,baseline_s0_error=za-zt,
        baseline_s0_relative_error_pct=100*(za/zt-1) if zt else float('nan'))

def load(suffix):
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf-8-sig'))
    data=[]
    for v in plan['variants']:
        name=v['id'];p=ROOT/'outputs/spectra'/(name+suffix+'.npz')
        if not p.exists():continue
        with np.load(p) as z:d={k:z[k] for k in z.files}
        with np.load(ROOT/'rcwa/cases'/name/'d4c_targets.npz') as z:
            d.update(s_um=z['s_um'],target_SUM=z['S_SUM'],target_DIFF=z['S_DIFF_same_kernel'],I0=float(z['I0']))
        data.append(dict(v,**d))
    return data

def features(d,h=0.,decomposition='band',tau=0.,NA=.6):
    zero=None if decomposition=='band' else d[f'A_zero_degree{decomposition}_V']
    return b_features(d['kx_rad_m'],d['A_ref_V'],float(d['Pin_W_m']),selected_orders=np.arange(1,19),
        stop_halfwidth_rad_um=h,stop_power_transmission=tau,A_zero_ref_V=zero,NA=NA)

def records(data,h=0.,decomposition='band',convention='K_over_gamma',model='scalar_zero_once'):
    rec=[];feat={}
    for d in data:
        f=features(d,h,decomposition);feat[d['id']]=f
        nz=scan_from_features(f,d['s_um'],convention=convention,K=K_VALUES)
        for channel in ('SUM','DIFF'):
            # Physical branch includes its measured passed power in nonzero,
            # while the fixed leakage multiplies the blocked zero power.
            if model=='physical_zero_once':
                nonzero=nz+f.zero_pass_power
                pzero=f.zero_blocked_power
            else:nonzero=nz;pzero=f.zero_raw_power
            rec.append(dict(id=d['id'],channel=channel,target=d['target_'+channel],nonzero=nonzero,
                zero_power=pzero,cap_nm=d['cap_nm'],s_um=d['s_um']))
    return rec,feat

def analyze(suffix='_fine'):
    data=load(suffix)
    if not data:raise RuntimeError('No processed FDTD spectra; run prepare_spectra first')
    out=ROOT/'outputs'/('comparison'+(suffix or '_coarse'));out.mkdir(parents=True,exist_ok=True)
    models=[('raw_EM_oldstop','band',.2,'none','scalar_zero_once',.6),
        ('weighted_oldstop_scalar','band',.2,'K_over_gamma','scalar_zero_once',.6),
        ('weighted_matchedstop_scalar','band',.17991979895042,'K_over_gamma','scalar_zero_once',.6),
        ('weighted_nostop_band_scalar','band',0.,'K_over_gamma','scalar_zero_once',.6),
        ('weighted_nostop_mode0_scalar',0,0.,'K_over_gamma','scalar_zero_once',.6),
        ('weighted_nostop_mode2_scalar',2,0.,'K_over_gamma','scalar_zero_once',.6),
        ('weighted_nostop_mode2_scalar_fullband',2,0.,'K_over_gamma','scalar_zero_once',.65),
        ('weighted_matchedstop_mode2_scalar',2,.17991979895042,'K_over_gamma','scalar_zero_once',.6),
        ('weighted_matchedstop_mode2_physical',2,.17991979895042,'K_over_gamma','physical_zero_once',.6),
        ('weighted_matchedstop_mode2_double',2,.17991979895042,'K_over_gamma','double_counted_zero',.6)]
    rows=[];zeros=[];curves={};s=data[0]['s_um'];curvekeys=[]
    for label,dec,h,conv,model,NA in models:
        for d in data:
            f=features(d,h,dec,NA=NA)
            for ch in ('SUM','DIFF'):
                y=candidate_scan(f,s,model=model,channel=ch,convention=conv,K=K_VALUES)
                rows.append(dict(id=d['id'],channel=ch,model=label,**metrics(d['target_'+ch],y)))
                key=f"{d['id']}__{ch}__{label}";curves[key]=y;curvekeys.append(key)
            zeros.append(dict(id=d['id'],model=label,decomposition=dec,h_rad_um=h,
                P0_over_Pin=f.zero_raw_power,RCWA_I0=d['I0'],P0_relative_error_pct=100*(f.zero_raw_power/d['I0']-1),
                P0pass_over_Pin=f.zero_pass_power,tau_eff=f.zero_pass_power/f.zero_raw_power,
                P0_band_over_Pin=f.zero_band_power,zero_nonzero_cross_raw=f.zero_nonzero_cross_raw_power,
                zero_nonzero_cross_pass=f.zero_nonzero_cross_pass_power))
    fits={};fitrows=[];required=[]
    # Structured ablation: exact same datasets/weights for all shared fits.
    for dec in ('band',0,2):
        for conv in ('none','K_over_gamma'):
            rec,_=records(data,decomposition=dec,convention=conv)
            train=[r for r in rec if r['cap_nm'] in (0,100)];hold=[r for r in rec if r['cap_nm'] in (150,200)]
            if not train:continue
            for r in rec:
                leak=.00138 if r['channel']=='SUM' else .000431
                # Diagnostic only: value needed by one mark at gain=1. It is
                # never used to generate its prediction or held-out curve.
                needed=(np.mean((r['target']-r['nonzero'])[:-1])/r['zero_power']-leak)/(1-leak)
                required.append(dict(id=r['id'],channel=r['channel'],decomposition=dec,convention=conv,
                    needed_T0_at_gain1=float(needed),supplied_T0=T0,
                    meaning='per-mark diagnostic, not an implemented per-mark calibration'))
            for free in (('zero_transmission',),('gain',),('gain','zero_transmission')):
                key=f"dec{dec}_{conv}_"+'+'.join(free)
                fit=fit_global_parameters(periodic_records(train),free=free,bounds={'gain':(.1,2.),'zero_transmission':(0.,.02)})
                fit['heldout']=evaluate_global_parameters(periodic_records(hold),fit['parameters']) if hold else None
                fits[key]=fit
                p=fit['parameters']
                for r in rec:
                    y=p['gain']*(r['nonzero']+zero_coefficient(r['channel'],p['zero_transmission'])*r['zero_power'])
                    fitrows.append(dict(id=r['id'],channel=r['channel'],fit=key,split='train' if r['cap_nm'] in (0,100) else 'heldout',**metrics(r['target'],y)))
                    curves[f"{r['id']}__{r['channel']}__fit_{key}"]=y
    # A real opaque stop uses its measured zero transmission and has no second
    # independent T0. Search one common h; optimize only a common power gain.
    sweep=[];best=None
    for h in np.arange(.14,.22500001,.0025):
        rec,feat=records(data,h=float(h),decomposition=2,model='physical_zero_once')
        train=[r for r in rec if r['cap_nm'] in (0,100)]
        fit=fit_global_parameters(periodic_records(train),free=('gain',),fixed={'zero_transmission':0.},bounds={'gain':(.5,1.5)})
        score=float(np.mean([r['relative_l2']**2 for r in fit['records']]))
        p=fit['parameters'];errs=[]
        for r in rec:
            y=p['gain']*(r['nonzero']+zero_coefficient(r['channel'],0.)*r['zero_power'])
            errs.append(dict(id=r['id'],channel=r['channel'],**metrics(r['target'],y)))
        entry=dict(h_rad_um=float(h),gain=p['gain'],train_mean_relative_L2_pct=100*np.sqrt(score),
            heldout_mean_relative_L2_pct=float(np.mean([e['relative_L2_pct'] for e in errs if next(d['cap_nm'] for d in data if d['id']==e['id']) in (150,200)])) if any(d['cap_nm'] in (150,200) for d in data) else None)
        sweep.append(entry)
        if best is None or score<best[0]:best=(score,entry,errs)
    for d in data:
        curves[f"{d['id']}__SUM__D4C"]=d['target_SUM'];curves[f"{d['id']}__DIFF__D4C"]=d['target_DIFF']
    np.savez_compressed(out/'scan_curves.npz',s_um=s,**curves)
    writecsv(out/'candidate_errors.csv',rows);writecsv(out/'zero_power_audit.csv',zeros)
    writecsv(out/'shared_fit_errors.csv',fitrows);writecsv(out/'shared_stop_sweep.csv',sweep)
    writecsv(out/'required_T0_diagnostic.csv',required)
    dump(out/'shared_fits.json',fits)
    dump(out/'summary.json',dict(suffix=suffix,available_structures=[d['id'] for d in data],
        all_eight_complete=len(data)==8,train='both geometries, cap0/cap100',heldout='both geometries, cap150/cap200',
        best_shared_physical_stop=best[1],best_shared_physical_stop_errors=best[2],
        primary_unfitted_model='weighted_nostop_mode2_scalar',no_per_mark_offsets=True,
        no_per_order_fitted_gain=True,all_integrals='angular domain, infinite image plane',
        K2_K4='user-authorized linear interpolation',T_m='user-confirmed |r_m|^2',
        fitted_h_is_grid_search=True,stop_grid_spacing_rad_um=.0025,
        interpretive_limit='Fits compare numerical/readout models, not an independently calibrated physical instrument'))
    print(f'Analyzed {len(data)} real structures -> {out}',flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--suffix',default='_fine');args=p.parse_args();analyze(args.suffix)
