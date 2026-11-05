"""Scientific plots from actual saved FDTD/D4C comparisons only."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures'
def readcsv(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=185,bbox_inches='tight');fig.savefig(OUT/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def main():
    comp=ROOT/'outputs/comparison_fine';summary=json.loads((comp/'summary.json').read_text(encoding='utf8'))
    ids=summary['available_structures'];curves=np.load(comp/'scan_curves.npz');s=curves['s_um']
    err=readcsv(comp/'candidate_errors.csv');fiterr=readcsv(comp/'shared_fit_errors.csv');zero=readcsv(comp/'zero_power_audit.csv')
    fits=json.loads((comp/'shared_fits.json').read_text(encoding='utf8'))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
    colors=['#205b89','#d77d27','#26847d','#9e4b87']
    keys=[('weighted_nostop_band_scalar','B0 band + scalar zero','#d77d27'),
          ('weighted_nostop_mode2_scalar','Mode A0 + scalar zero','#26847d'),
          ('weighted_matchedstop_mode2_double','Repeated passed zero','#9e4b87')]
    fig,axes=plt.subplots(4,2,figsize=(12,12),layout='constrained',sharex=True)
    for ax,name in zip(axes.T.flat,ids):
        ax.plot(s,curves[f'{name}__SUM__D4C'],color='#202937',lw=2,label='D4C target')
        for key,label,col in keys[:2]:ax.plot(s,curves[f'{name}__SUM__{key}'],color=col,lw=1,ls='--',label=label)
        ax.set_title(name.replace('_cap',' / cap ')+' nm');ax.set_ylabel('Signal / incident power')
    for ax in axes[-1]:ax.set_xlabel('Displacement s / um')
    for ax in list(axes.T.flat)[len(ids):]:ax.set_visible(False)
    axes[0,0].legend(fontsize=8);fig.suptitle('Actual finite-Gaussian FDTD vs D4C: absolute curves, no fitted gain/offset')
    save(fig,'06_all_absolute_scans')
    fig,axes=plt.subplots(4,2,figsize=(12,11),layout='constrained',sharex=True)
    for ax,name in zip(axes.T.flat,ids):
        t=curves[f'{name}__SUM__D4C'];span=np.ptp(t)
        for key,label,col in keys:ax.plot(s,100*(curves[f'{name}__SUM__{key}']-t)/span,color=col,label=label,lw=1.1)
        ax.axhline(0,color='black',lw=.6);ax.set_title(name.replace('_cap',' / cap ')+' nm');ax.set_ylabel('Residual / target span (%)')
    for ax in axes[-1]:ax.set_xlabel('Displacement s / um')
    for ax in list(axes.T.flat)[len(ids):]:ax.set_visible(False)
    axes[0,0].legend(fontsize=8);fig.suptitle('Residuals retain both baseline and waveform errors')
    save(fig,'07_scan_residuals')
    labels=[x.replace('ah53_cap','AH').replace('binary50_cap','B') for x in ids]
    variants=[('raw_EM_oldstop','Unweighted EM'),('weighted_oldstop_scalar','K/gamma, old stop'),
        ('weighted_nostop_band_scalar','B0 + scalar'),('weighted_nostop_mode2_scalar','Mode A0 + scalar'),
        ('weighted_matchedstop_mode2_physical','Mode A0 + physical once'),('weighted_matchedstop_mode2_double','Repeated passed zero')]
    matrix=np.array([[float(next(r['relative_L2_pct'] for r in err if r['id']==name and r['channel']=='SUM' and r['model']==m)) for name in ids] for m,_ in variants])
    fig,ax=plt.subplots(figsize=(12,5),layout='constrained')
    im=ax.imshow(np.log10(np.maximum(matrix,.001)),cmap='YlOrRd',aspect='auto')
    ax.set(xticks=np.arange(len(ids)),xticklabels=labels,yticks=np.arange(len(variants)),yticklabels=[l for _,l in variants],title='Absolute-curve relative L2 error (%); all models have gain = 1')
    for i in range(len(variants)):
        for j in range(len(ids)):ax.text(j,i,f'{matrix[i,j]:.3g}',ha='center',va='center',fontsize=10,color='black')
    fig.colorbar(im,ax=ax,label='log10(error in %)');save(fig,'08_model_comparison')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for dec,marker,label in [('band','o','B0 angular band'),('2','s','Source-only modal A0')]:
        rr=[next(r for r in zero if r['id']==name and r['model']==('weighted_nostop_band_scalar' if dec=='band' else 'weighted_nostop_mode2_scalar')) for name in ids]
        axes[0].plot(labels,[100*(float(r['P0_over_Pin'])/float(r['RCWA_I0'])-1) for r in rr],marker+'-',label=label)
        phys=[next(r for r in zero if r['id']==name and r['model']==('weighted_matchedstop_scalar' if dec=='band' else 'weighted_matchedstop_mode2_physical')) for name in ids]
        axes[1].plot(labels,[float(r['tau_eff']) for r in phys],marker+'-',label=label)
    axes[0].set(ylabel='Zero power relative to RCWA I0 (%)',title='Independent FDTD zero-power extraction')
    axes[1].axhline(.000629728295091705,color='black',ls='--',label='Supplied T0')
    axes[1].set(ylabel='Measured passed / raw zero power',title='One common opaque stop h = 0.17992 rad/um')
    for ax in axes:ax.tick_params(axis='x',rotation=25);ax.legend(fontsize=8)
    save(fig,'09_zero_power_and_transmission')
    selected=['decband_K_over_gamma_zero_transmission','decband_K_over_gamma_gain+zero_transmission',
              'dec2_K_over_gamma_zero_transmission','dec2_K_over_gamma_gain+zero_transmission']
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for fit,col in zip(selected,colors):
        rr=[next(r for r in fiterr if r['id']==name and r['channel']=='SUM' and r['fit']==fit) for name in ids]
        label=('B0' if 'decband' in fit else 'Mode A0')+(', gain + T0' if 'gain+' in fit else ', T0 only')
        axes[0].plot(labels,[float(r['relative_L2_pct']) for r in rr],'o-',color=col,label=label)
        axes[1].plot(labels,[float(r['mean_error_over_target_span_pct']) for r in rr],'o-',color=col,label=label)
    for ax in axes:
        ax.tick_params(axis='x',rotation=25);ax.legend(fontsize=8)
        for j,name in enumerate(ids):
            if name.endswith(('150','200')):ax.axvspan(j-.4,j+.4,color='#9abdb0',alpha=.13)
    axes[0].set(ylabel='Absolute-curve relative L2 error (%)',title='Shared parameters; shaded caps are held out')
    axes[1].set(ylabel='Mean residual / target span (%)',title='Absolute vertical displacement after shared calibration')
    save(fig,'10_global_fit_validation')
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    for ax,ch in zip(axes,('SUM','DIFF')):
        fixed=[next(r for r in err if r['id']==name and r['channel']==ch and r['model']=='weighted_nostop_mode2_scalar') for name in ids]
        calibrated=[next(r for r in fiterr if r['id']==name and r['channel']==ch and r['fit']=='dec2_K_over_gamma_gain+zero_transmission') for name in ids]
        repeated=[next(r for r in err if r['id']==name and r['channel']==ch and r['model']=='weighted_matchedstop_mode2_double') for name in ids]
        ax.plot(labels,[float(r['baseline_s0_relative_error_pct']) for r in fixed],'o-',color='#26847d',label='Fixed supplied T0')
        ax.plot(labels,[float(r['baseline_s0_relative_error_pct']) for r in calibrated],'s-',color='#9e4b87',label='Whole-curve gain + T0 fit')
        ax.plot(labels,[float(r['baseline_s0_relative_error_pct']) for r in repeated],'^-',color='#d77d27',label='Passed zero counted again')
        ax.axhline(0,color='black',lw=.7);ax.set(title=ch+' leakage version: actual s = 0 baseline',ylabel='Baseline relative error (%)')
        ax.tick_params(axis='x',rotation=25);ax.legend(fontsize=8)
    save(fig,'13_zero_baseline_validation')
    physical=ROOT/'outputs/physical_component_fine'
    if (physical/'shared_fit_errors.csv').exists():
        pr=readcsv(physical/'shared_fit_errors.csv')
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
        fixed=[next(r for r in err if r['id']==name and r['channel']=='SUM' and r['model']=='weighted_nostop_mode2_scalar') for name in ids]
        ac=[next(r for r in pr if r['id']==name and r['channel']=='SUM' and r['method']=='actual_nearest_band_NA0.60' and r['gain_scope']=='nonzero_only' and r['free']=='gain_ac' and r['objective']=='curve_only') for name in ids]
        aware={ch:[next(r for r in pr if r['id']==name and r['channel']==ch and r['method']=='actual_nearest_band_NA0.60' and r['gain_scope']=='nonzero_only' and r['free']=='gain_ac+zero_transmission' and r['objective']=='curve_and_s0_equal_weight') for name in ids] for ch in ('SUM','DIFF')}
        axes[0].plot(labels,[float(r['relative_L2_pct']) for r in fixed],'o-',label='Fixed all parameters',color='#26847d')
        axes[0].plot(labels,[float(r['relative_L2_pct']) for r in ac],'s-',label='One shared nonzero gain; fixed T0',color='#205b89')
        axes[0].plot(labels,[float(r['relative_L2_pct']) for r in aware['SUM']],'^--',label='Baseline-aware gain + T0',color='#9e4b87')
        axes[1].plot(labels,[float(r['baseline_s0_relative_error_pct']) for r in fixed],'s-',label='Fixed T0, with/without nonzero gain',color='#205b89')
        for ch,marker,col in [('SUM','o','#9e4b87'),('DIFF','^','#d77d27')]:
            axes[1].plot(labels,[float(r['baseline_relative_pct']) for r in aware[ch]],marker+'--',label='Baseline-aware '+ch,color=col)
        for ax in axes:
            ax.tick_params(axis='x',rotation=25);ax.legend(fontsize=8)
            for j,name in enumerate(ids):
                if name.endswith(('150','200')):ax.axvspan(j-.4,j+.4,color='#9abdb0',alpha=.13)
        axes[0].set(title='One common parameter across 8 structures',ylabel='Absolute-curve relative L2 error (%)')
        axes[1].set(title='Zero baseline remains an independent criterion',ylabel='s = 0 relative baseline error (%)')
        save(fig,'14_recommended_shared_readout')
    sweep=readcsv(comp/'shared_stop_sweep.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    h=[float(r['h_rad_um']) for r in sweep]
    axes[0].plot(h,[float(r['train_mean_relative_L2_pct']) for r in sweep],label='Training RMS')
    if all(r['heldout_mean_relative_L2_pct'] for r in sweep):axes[0].plot(h,[float(r['heldout_mean_relative_L2_pct']) for r in sweep],label='Held-out mean')
    axes[0].axvline(.17991979895,color='black',ls='--',label='Supplied T0 / ideal Gaussian')
    axes[0].set(xlabel='Common stop halfwidth h / rad per um',ylabel='Relative curve error (%)',title='Physical zero counted once; one shared gain')
    axes[0].legend(fontsize=8)
    axes[1].plot(h,[float(r['gain']) for r in sweep]);axes[1].set(xlabel='Common h / rad per um',ylabel='Fitted common power gain',title='Fitted only on cap 0 and cap 100 nm')
    save(fig,'11_shared_stop_width_sweep')
    # Mesh/time results are from independent real solves, not error bars on
    # inferred parameters. A control can be missing during an interim build.
    audit=[]
    for suff in ('','_fine','_yfine'):
        path=ROOT/'fdtd'/('baseline_flux_audit'+suff+'.json')
        if path.exists():audit.append(json.loads(path.read_text(encoding='utf8')))
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    controlpath=ROOT/'outputs/controls/mesh_time_summary.json'
    if controlpath.exists():
        controlled=json.loads(controlpath.read_text(encoding='utf8'))['flat_Si_rows']
        axes[0].bar(['20/10 nm','10/5 nm','10/2.5 nm'],[x['relative_error_pct'] for x in controlled],color='#205b89')
    elif audit:
        axes[0].bar([x.get('suffix') or 'coarse' for x in audit],[100*x['flat_R_relative_error'] for x in audit],color='#205b89')
    axes[0].set(ylabel='Flat-Si R error against Fresnel (%)',title='Absolute numerical accuracy, no fitted correction')
    for suff,col,label in [('',colors[1],'dx20 / dy10 nm'),('_fine',colors[0],'dx10 / dy5 nm')]:
        p=ROOT/'outputs'/('comparison'+(suff or '_coarse'))/'scan_curves.npz'
        if p.exists():
            with np.load(p) as z:
                name='ah53_cap000__SUM__weighted_nostop_mode2_scalar';tar='ah53_cap000__SUM__D4C'
                axes[1].plot(z['s_um'],100*(z[name]-z[tar])/np.ptp(z[tar]),color=col,label=label)
    axes[1].set(xlabel='s / um',ylabel='Residual / D4C span (%)',title='Same readout formula at two meshes');axes[1].legend(fontsize=8)
    save(fig,'04_numerical_accuracy')
    # Central angular spectrum from a true FDTD solve, including recovered A0.
    name='ah53_cap000'
    with np.load(ROOT/'outputs/spectra'/(name+'_fine.npz')) as z:
        k=z['kx_rad_m']*1e-6;A=z['A_ref_V'];A0=z['A_zero_degree2_V'];G=2*np.pi/16
        choose=np.abs(k)<1.05*G;scale=np.max(np.abs(A0))
        fig,axes=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
        for ax in axes:
            ax.semilogy(k[choose],np.abs(A[choose])/scale,color='#205b89',label='Actual total FDTD')
            ax.semilogy(k[choose],np.abs(A0[choose])/scale,color='#26847d',ls='--',label='Source-template A0 estimate')
            ax.semilogy(k[choose],np.abs((A-A0)[choose])/scale,color='#d77d27',ls=':',label='Nonzero remainder')
            ax.axvspan(-G/2,G/2,color='gray',alpha=.1,label='B0 angular band')
            ax.axvline(.17991979895,color='black',ls='--');ax.axvline(-.17991979895,color='black',ls='--')
            ax.set(xlabel='kx / rad per um',ylabel='Electric spectrum modulus / A0 peak')
        axes[0].set(title='True finite-Gaussian spectrum: AH53, bare Si',ylim=(1e-6,2))
        axes[1].set(title='Band boundary and stop boundary are different',xlim=(.13,.26),ylim=(1e-5,.1))
        axes[0].legend(fontsize=8);save(fig,'05_zero_angular_spectrum')
    from plot_d4c_minima import main as minima_plot
    minima_plot()
    print(f'Generated actual-data figures for {len(ids)} structures')

if __name__=='__main__':main()
