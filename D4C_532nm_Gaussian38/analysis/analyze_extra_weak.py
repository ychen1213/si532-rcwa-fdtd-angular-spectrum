"""Ninth weak-signal holdout, with all original calibration parameters frozen."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from run_analysis import ROOT,features,metrics
from zero_order_models import candidate_scan,scan_from_features,zero_coefficient

ID='ah53_h180_cap400'
def main():
    out=ROOT/'outputs/extra_weak';out.mkdir(exist_ok=True,parents=True)
    # 100fs failed the weak signal's OWN span criterion. The 200fs run is
    # the candidate production record, independently checked against 400fs.
    primary_suffix='_fine_t200'
    with np.load(ROOT/'outputs/spectra'/(ID+primary_suffix+'.npz')) as z:d={k:z[k] for k in z.files}
    targetfile=ROOT/'rcwa/cases'/ID/'d4c_targets.npz'
    with np.load(targetfile) as z:s=z['s_um'];targets={'SUM':z['S_SUM'],'DIFF':z['S_DIFF_same_kernel']}
    calibration=ROOT/'outputs/physical_component_fine/shared_fits.json'
    pfit=json.loads(calibration.read_text(encoding='utf8'))
    method='actual_nearest_band_NA0.60'
    good=pfit[method+'__nonzero_only__gain_ac']['parameters']
    aware=pfit[method+'__nonzero_only_baseline_aware__gain_ac+zero_transmission']['parameters']
    originalfits=json.loads((ROOT/'outputs/comparison_fine/shared_fits.json').read_text(encoding='utf8'))
    wrong=originalfits['dec2_K_over_gamma_gain+zero_transmission']['parameters']
    no=features(d,decomposition=2);stop=features(d,h=.17991979895041987,decomposition=2)
    band=features(d,decomposition='band')
    nz=scan_from_features(no,s,convention='K_over_gamma')
    rows=[];curves={};plain=[]
    for ch,t in targets.items():
        curves['D4C_'+ch]=t
        candidates={
            'B0_band_scalar':candidate_scan(band,s,channel=ch,model='scalar_zero_once',convention='K_over_gamma'),
            'fixed_T0_once':candidate_scan(no,s,channel=ch,model='scalar_zero_once',convention='K_over_gamma'),
            'physical_stop_once':candidate_scan(stop,s,channel=ch,model='physical_zero_once',convention='K_over_gamma'),
            'repeated_passed_zero':candidate_scan(stop,s,channel=ch,model='double_counted_zero',convention='K_over_gamma'),
            'frozen_nonzero_gain':good['gain_ac']*nz+zero_coefficient(ch,good['zero_transmission'])*no.zero_raw_power,
            'frozen_baseline_aware':aware['gain_ac']*nz+zero_coefficient(ch,aware['zero_transmission'])*no.zero_raw_power,
            'frozen_curve_only_gain_T0':wrong['gain']*(nz+zero_coefficient(ch,wrong['zero_transmission'])*no.zero_raw_power)}
        for model,y in candidates.items():
            rows.append(dict(id=ID,channel=ch,model=model,**metrics(t,y)))
            curves[model+'_'+ch]=y
    np.savez_compressed(out/'scan_curves.npz',s_um=s,**curves)
    with (out/'errors.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for j,ss in enumerate(s):plain.append(dict(s_um=float(ss),**{k:float(v[j]) for k,v in curves.items()}))
    with (out/'scan_curves.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(plain[0]));w.writeheader();w.writerows(plain)
    summary=dict(id=ID,extra_independent_holdout=True,training_additions=False,
        target_curve_file=str(targetfile.relative_to(ROOT)),target_curve_sha256=hashlib.sha256(targetfile.read_bytes()).hexdigest(),
        primary_case_suffix=primary_suffix,
        spectrum_file='outputs/spectra/'+ID+primary_suffix+'.npz',spectrum_sha256=hashlib.sha256((ROOT/'outputs/spectra'/(ID+primary_suffix+'.npz')).read_bytes()).hexdigest(),
        original_curve_only_calibration_file='outputs/comparison_fine/shared_fits.json',
        original_curve_only_calibration_sha256=hashlib.sha256((ROOT/'outputs/comparison_fine/shared_fits.json').read_bytes()).hexdigest(),
        structure=json.loads((ROOT/'inputs/extra_weak_mark.json').read_text(encoding='utf8')),
        target_summary={ch:dict(minimum=float(t.min()),maximum=float(t.max()),span=float(np.ptp(t))) for ch,t in targets.items()},
        frozen_calibration_file=str(calibration.relative_to(ROOT)),calibration_sha256=hashlib.sha256(calibration.read_bytes()).hexdigest(),
        frozen_nonzero_gain_parameters=good,frozen_baseline_aware_parameters=aware,frozen_failed_curve_only_parameters=wrong,
        extra_mark_was_not_used_for_any_fit=True,errors=rows,
        physical_zero_over_Pin=no.zero_raw_power,physical_zero_after_stop_over_Pin=stop.zero_pass_power,
        interpretive_limit='Weak scan ordinate is not the total or zeroth reflectance; this extra mark also changes ridge height to180nm')
    controls=[];control_curves={}
    for suffix in ('_fine','_fine_t200','_fine_t400'):
        path=ROOT/'outputs/spectra'/(ID+suffix+'.npz')
        if not path.exists():continue
        with np.load(path) as z:control_data={key:z[key] for key in z.files}
        cf=features(control_data,decomposition=2)
        done=json.loads((ROOT/'fdtd/cases'/(ID+suffix)/'DONE.json').read_text(encoding='utf8'))
        for ch,t in targets.items():
            cy=candidate_scan(cf,s,channel=ch,model='scalar_zero_once',convention='K_over_gamma')
            controls.append(dict(suffix=suffix,channel=ch,dy_nm=done['settings']['dy_m']*1e9,
                actual_time_fs=done['termination']['simulated_time_s']*1e15,**metrics(t,cy)))
            control_curves[suffix+'_'+ch]=cy
    summary['controls_complete']=len(controls)==6
    summary['time_mesh_control_rows']=controls
    if '_fine_t200_SUM' in control_curves:
        summary['weak_time_max_difference_over_own_span_pct']={ch:float(100*np.max(np.abs(control_curves['_fine_'+ch]-control_curves['_fine_t200_'+ch]))/np.ptp(t)) for ch,t in targets.items()}
    if '_fine_t400_SUM' in control_curves:
        summary['weak_200_vs_400_max_difference_over_own_span_pct']={ch:float(100*np.max(np.abs(control_curves['_fine_t200_'+ch]-control_curves['_fine_t400_'+ch]))/np.ptp(t)) for ch,t in targets.items()}
        summary['time_control_passed']=max(summary['weak_200_vs_400_max_difference_over_own_span_pct'].values())<.2
    np.savez_compressed(out/'control_curves.npz',s_um=s,**control_curves)
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    styles=[('fixed_T0_once','Fixed T0, zero once','#26847d','--'),
            ('B0_band_scalar','B0 angular-band zero','#8d6746',':'),
            ('frozen_nonzero_gain','Original frozen nonzero gain','#205b89',':'),
            ('repeated_passed_zero','Passed zero counted twice','#d77d27','-.'),
            ('frozen_curve_only_gain_T0','Original curve-only T0 fit','#9e4b87','--')]
    for j,ch in enumerate(('SUM','DIFF')):
        ax=axes[0,j];t=targets[ch];ax.plot(s,t,color='black',lw=2,label='D4C target')
        for model,label,col,ls in styles:ax.plot(s,curves[model+'_'+ch],color=col,ls=ls,lw=1.2,label=label)
        ax.set(title=ch+' leakage: extra weak scan, no refitting',xlabel='s / um',ylabel='Absolute normalized signal')
        ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.grid(alpha=.2)
        ax=axes[1,j];sel=np.abs(s)<.12;ax.plot(s[sel]*1000,t[sel],color='black',lw=2,label='D4C target')
        for model,label,col,ls in styles:ax.plot(s[sel]*1000,curves[model+'_'+ch][sel],color=col,ls=ls,lw=1.2,label=label)
        ax.set(title=f'Nonzero bottom: D4C min = {t.min():.7g}',xlabel='Local s / nm',ylabel='Absolute normalized signal')
        ax.grid(alpha=.2)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncols=3,fontsize=8)
    fig.suptitle('Ninth holdout: AH53 h180nm + planar SiO2 cap400nm; same 38um beam and original parameters',fontsize=12)
    for ext in ('png','svg'):fig.savefig(ROOT/'figures'/('17_extra_weak_scan.'+ext),dpi=190,bbox_inches='tight')
    plt.close(fig)
    if summary['controls_complete']:
        fig,axes=plt.subplots(2,2,figsize=(12,7.7),layout='constrained')
        for j,ch in enumerate(('SUM','DIFF')):
            t=targets[ch];axes[0,j].plot(s,t,'k-',lw=2,label='q0 D4C target')
            for suffix,col,label in [('_fine','#205b89','dy5nm / 100fs (failed)'),('_fine_t200','#d77d27','dy5nm / 200fs'),('_fine_t400','#26847d','dy5nm / 400fs')]:
                y=control_curves[suffix+'_'+ch]
                axes[0,j].plot(s,y,color=col,label=label,lw=1.1)
                axes[1,j].plot(s,100*(y-t)/np.ptp(t),color=col,label=label,lw=1.1)
            axes[0,j].set(title=ch+': fixed T0, no empirical gain',xlabel='s / um',ylabel='Absolute normalized signal',ylim=(0,None))
            axes[1,j].set(title='Residual normalized to this weak scan span',xlabel='s / um',ylabel='Error / weak D4C span (%)')
            for ax in axes[:,j]:ax.grid(alpha=.2)
        handles,labels=axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='outside lower center',ncols=4,fontsize=8)
        for ext in ('png','svg'):fig.savefig(ROOT/'figures'/('18_extra_weak_controls.'+ext),dpi=190,bbox_inches='tight')
        plt.close(fig)
    print(json.dumps(dict(targets=summary['target_summary'],controls_complete=summary['controls_complete'],
        errors=[{k:r[k] for k in ('channel','model','relative_L2_pct','baseline_s0_relative_error_pct')} for r in rows]),indent=2))

if __name__=='__main__':main()
