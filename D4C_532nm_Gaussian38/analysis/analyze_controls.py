"""Independent time and vertical-mesh controls with the same fixed readout."""
from pathlib import Path
import json
import numpy as np
from prepare_spectra import read_monitor,ROOT,K0,Z0
from run_analysis import metrics,features
from zero_order_models import candidate_scan

def main():
    out=ROOT/'outputs/controls';out.mkdir(parents=True,exist_ok=True)
    d={};rows=[];flat=[]
    with np.load(ROOT/'rcwa/cases/ah53_cap150/d4c_targets.npz') as z:
        s=z['s_um'];target=z['S_SUM']
    for suffix in ('_fine','_fine_t200','_finerY'):
        p=ROOT/'outputs/spectra'/('ah53_cap150'+suffix+'.npz')
        if not p.exists():raise RuntimeError(f'Control is not prepared: {p.name}')
        with np.load(p) as z:
            a={key:z[key] for key in z.files}
        f=features(a,decomposition=2)
        y=candidate_scan(f,s,model='scalar_zero_once',channel='SUM',convention='K_over_gamma')
        done=json.loads((ROOT/'fdtd/cases'/('ah53_cap150'+suffix)/'DONE.json').read_text(encoding='utf8'))
        rows.append(dict(suffix=suffix,dy_nm=done['settings']['dy_m']*1e9,dx_nm=done['settings']['dx_m']*1e9,
            actual_time_fs=done['termination']['simulated_time_s']*1e15,**metrics(target,y)))
        d[suffix]=y
    rf=(1-complex(4.149971915645845,.048064634885764825))/(1+complex(4.149971915645845,.048064634885764825))
    for suffix in ('','_fine','_finerY'):
        air=ROOT/'fdtd/cases'/('air'+suffix);case=ROOT/'fdtd/cases'/('flatSi'+suffix)
        i=read_monitor(air,'INC');r=read_monitor(case,'R');empty=read_monitor(air,'R')
        k=r['k'];g=r['gamma'];dk=k[1]-k[0]
        power=np.sum(g*np.abs(r['up']-empty['up'])**2)*dk/(4*np.pi*Z0)
        R=power/i['metadata']['down_power_W_m'];error=100*(R/abs(rf)**2-1)
        settings=json.loads((case/'settings.json').read_text(encoding='utf8'))
        flat.append(dict(suffix=suffix,R=float(R),Fresnel_R=abs(rf)**2,relative_error_pct=float(error),
            dy_nm=settings['dy_m']*1e9))
    y100=d['_fine'];y200=d['_fine_t200'];span=np.ptp(target)
    time=dict(max_curve_difference_over_D4C_span_pct=float(100*np.max(np.abs(y100-y200))/span),
        max_curve_difference=float(np.max(np.abs(y100-y200))),
        baseline_relative_change_pct=float(100*(y200[len(s)//2]/y100[len(s)//2]-1)))
    np.savez_compressed(out/'mesh_time_curves.npz',s_um=s,D4C_SUM=target,**d)
    result=dict(film='ah53 cap150',scan_rows=rows,flat_Si_rows=flat,time_comparison=time,
        fixed_readout='K_m/gamma_m, source-template A0, original T0/leakage, no fitted gain or offset',
        PML_limitation='finerY uses actual Custom128 versus fine Standard64; Y physical thickness retained, X thickness doubled; not a pure single-parameter PML test')
    (out/'mesh_time_summary.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
    fig,axes=plt.subplots(1,3,figsize=(14,4.4),layout='constrained')
    axes[0].bar(['20/10 nm','10/5 nm','10/2.5 nm'],[r['relative_error_pct'] for r in flat],color=['#c49b65','#205b89','#26847d'])
    axes[0].set(title='Flat-Si R: independent Fresnel control',ylabel='Relative reflected-power error (%)')
    for suffix,col,label in [('_fine','#205b89','dy5nm / 100fs'),('_fine_t200','#d77d27','dy5nm / 200fs'),('_finerY','#26847d','dy2.5nm / 100fs')]:
        axes[1].plot(s,100*(d[suffix]-target)/span,color=col,label=label)
    axes[1].set(title='AH53 cap150: same fixed readout',xlabel='s / um',ylabel='Residual / D4C span (%)');axes[1].legend(fontsize=8)
    x=np.arange(3);axes[2].bar(x,[r['relative_L2_pct'] for r in rows],color=['#205b89','#d77d27','#26847d'])
    axes[2].set(xticks=x,xticklabels=['dy5/100fs','dy5/200fs','dy2.5/100fs'],title='Whole curve before any gain fit',ylabel='Relative L2 error (%)')
    axes[2].tick_params(axis='x',rotation=15)
    for ext in ('png','svg'):fig.savefig(ROOT/'figures'/('15_film_mesh_and_time.'+ext),dpi=190,bbox_inches='tight')
    plt.close(fig)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
