"""Traceable geometry, complete D4C targets and actual RCWA convergence plots."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from te_planarized_rcwa import dump
from run_reference import csvwrite

ROOT=Path(__file__).resolve().parents[1];FIG=ROOT/'figures';HERE=ROOT/'rcwa'
COLORS={0:'#245b8a',100:'#bd6329',150:'#36846a',200:'#9a4b7c'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
    'axes.labelsize':10,'legend.fontsize':9,'svg.fonttype':'none'})


def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=210,bbox_inches='tight')
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight');plt.close(fig)


def cross_section(ax,v,title):
    cap=v['cap_nm'];top=300+cap
    ax.add_patch(Rectangle((-8,-100),16,100,facecolor='#586b7a',edgecolor='none'))
    if cap:ax.add_patch(Rectangle((-8,0),16,top,facecolor='#a9d8ed',edgecolor='none'))
    for xc in v['ridge_centers_local_um']:
        ax.add_patch(Rectangle((xc-v['ridge_width_um']/2,0),v['ridge_width_um'],300,facecolor='#586b7a',edgecolor='none'))
    ax.axhline(0,color='#34495a',lw=.8);ax.set(xlim=(-8,8),ylim=(-100,550),xlabel='x / um',ylabel='y / nm',title=title)
    ax.set_yticks([0,300,450]);ax.set_xticks([-8,-4,0,4,8]);ax.grid(axis='y',alpha=.16)
    ax.text(.98,.96,'Air',ha='right',va='top',transform=ax.transAxes)
    ax.text(.02,.04,'Si substrate',ha='left',va='bottom',color='white',fontsize=9,transform=ax.transAxes)


def main():
    FIG.mkdir(exist_ok=True)
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf8'))
    variants=plan['variants'];byid={v['id']:v for v in variants}
    summary=json.loads((HERE/'reference_summary.json').read_text(encoding='utf8'))
    with (HERE/'convergence.csv').open(encoding='utf-8-sig',newline='') as f:conv=list(csv.DictReader(f))
    fig,aa=plt.subplots(2,2,figsize=(9.4,7.1),layout='constrained')
    cross_section(aa[0,0],byid['ah53_cap150'],'(a) AH53 mother cell; 150 nm cap example')
    cross_section(aa[0,1],byid['binary50_cap150'],'(b) Binary 50% mother cell; same film')
    x=np.linspace(-82,82,2001);intensity=np.exp(-2*(x/19)**2);ax=aa[1,0]
    for cell in range(-4,5):
        ax.axvline(cell*16-8,color='#8594a1',lw=.4,alpha=.5)
        for local in (-3.2,0,3.2):
            ax.add_patch(Rectangle((cell*16+local-.8,-.13),1.6,.10,facecolor='#586b7a',edgecolor='none'))
    ax.axvline(72,color='#8594a1',lw=.4);ax.plot(x,intensity,color='#bd6329',lw=1.8,label='Gaussian intensity')
    ax.axvline(-19,color='#bd6329',ls=':',lw=1);ax.axvline(19,color='#bd6329',ls=':',lw=1)
    ax.set(xlim=(-82,82),ylim=(-.16,1.06),xlabel='x / um',ylabel='I(x) / I(0)',title='(c) Finite 9 cells (144 um) and 38 um spot')
    ax.text(-78,.76,'38 um: intensity\n1/e² diameter',ha='left',fontsize=8.7)
    ax.text(-78,.55,'Gray bars: Si ridges\n(height schematic)',ha='left',fontsize=8)
    ax.grid(alpha=.16);ax.legend(loc='upper right')
    ax=aa[1,1];ax.axis('off');ax.set_title('(d) Two geometries × four film conditions = 8 marks')
    cells=[[str(cap),'Air' if cap==0 else 'SiO2',str(300+cap)] for cap in (0,100,150,200)]
    tab=ax.table(cellText=cells,colLabels=['Cap / nm','Groove fill','Top y / nm'],cellLoc='center',bbox=[.02,.33,.96,.57])
    tab.auto_set_font_size(False);tab.set_fontsize(10)
    ax.text(.02,.22,'Cap 0: no SiO2 anywhere.\nCap > 0: SiO2 fills grooves and covers Si ridges.\nSi relief height = 300 nm; mother period = 16 um.',va='top',fontsize=9.5,linespacing=1.5)
    fig.suptitle('Geometry and illumination used for the matched RCWA / finite FDTD study',fontsize=12)
    save(fig,'01_geometry')
    csvwrite(FIG/'01_geometry.csv',[dict(x_um=float(a),normalized_Gaussian_intensity=float(b)) for a,b in zip(x,intensity)])
    dump(FIG/'01_geometry.json',dict(variants=variants,waist_radius_um=19.,intensity_1e2_diameter_um=38.,
        geometry_scope='Panels a/b are one mother cell at cap150; panel c shows nine AH53 cells; lengths are actual but x and y have different units/scales.'))

    metrics=[];curves={};wide={};fig,aa=plt.subplots(2,2,figsize=(9.5,7.4),layout='constrained')
    for col,geometry in enumerate(('ah53','binary50')):
        for cap in (0,100,150,200):
            name=f'{geometry}_cap{cap:03d}'
            with np.load(HERE/'cases'/name/'d4c_targets.npz') as z:
                s=z['s_um'].copy();y=z['S_SUM'].copy();d=z['S_DIFF_same_kernel'].copy();ac=z['AC'].copy()
                zero=float(z['zero_SUM']);zd=float(z['zero_DIFF'])
            curves[name]=dict(s_um=s,S_SUM=y,S_DIFF_same_kernel=d,AC=ac)
            yy=y[:-1];dd=d[:-1]
            metrics.append(dict(variant=name,S_SUM_min=float(yy.min()),S_SUM_max=float(yy.max()),
                S_SUM_period_mean=float(yy.mean()),S_SUM_peak_to_peak=float(np.ptp(yy)),
                S_SUM_DOM=float(np.ptp(yy)/(yy.max()+yy.min())),
                S_DIFF_same_kernel_min=float(dd.min()),S_DIFF_same_kernel_max=float(dd.max()),
                zero_SUM=zero,zero_DIFF=zd,AC_period_mean=float(ac[:-1].mean())))
            for row in (0,1):
                aa[row,col].plot(s,y,color=COLORS[cap],lw=1.4,label=f'{cap} nm cap')
            wide[name+'_SUM']=y;wide[name+'_DIFF_same_kernel']=d
        aa[0,col].set(title=('(a) AH53' if col==0 else '(b) Binary 50%')+' — full scan',xlim=(-4,4),xlabel='s / um',ylabel='D4C model signal')
        local=abs(s)<=.12
        family=[v['S_SUM'] for name,v in curves.items() if name.startswith(geometry+'_')]
        aa[1,col].set(title=('(c) AH53' if col==0 else '(d) Binary 50%')+' — near s = 0',xlim=(-.12,.12),yscale='log',
            ylim=(.7*min(float(y.min()) for y in family),1.4*max(float(y[local].max()) for y in family)),
            xlabel='s / um',ylabel='D4C model signal (log)')
        for row in (0,1):aa[row,col].grid(alpha=.2,which='both')
    handles,labels=aa[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=4,fontsize=9)
    fig.suptitle('Complete 8-mark D4C targets: supplied dark kernel + SUM zero leakage\nT_m = |r_m|²; user-authorized K2/K4 interpolation; no fitted FDTD scale',fontsize=11.5)
    save(fig,'02_d4c_targets')
    csvwrite(FIG/'02_d4c_targets.csv',[dict(s_um=float(ss),**{k:float(v[j]) for k,v in wide.items()}) for j,ss in enumerate(s)])
    csvwrite(HERE/'d4c_target_metrics.csv',metrics)
    dump(FIG/'02_d4c_targets.json',dict(primary_channel='SUM leakage with supplied 1-cos interference kernel',
        source='rcwa/cases/<variant>/d4c_targets.npz::S_SUM',scan_samples=len(s),
        statistics_endpoint='First 1600 of 1601 points; duplicate final endpoint excluded for period statistics.',
        K2=.4947,K4=.49135,K_interpolation='Explicitly authorized by user in this task.',metrics=metrics))

    fig,aa=plt.subplots(1,3,figsize=(10,3.7),layout='constrained')
    for geometry in ('ah53','binary50'):
        for cap in (0,100,150,200):
            name=f'{geometry}_cap{cap:03d}';rows=[r for r in conv if r['variant']==name and r['relative_complex_L2']]
            xx=np.array([int(r['max_order']) for r in rows]);marker='o' if geometry=='ah53' else 's'
            ref=next(r for r in summary['structures'] if r['variant']==name)
            vals=[np.array([float(r['relative_complex_L2'])*100 for r in rows]),
                np.array([float(r['zero_intensity_abs_change'])/ref['I0']*100 for r in rows]),
                np.array([float(r['max_complex_abs_change']) for r in rows])]
            for ax,yy in zip(aa,vals):ax.semilogy(xx,yy,marker=marker,color=COLORS[cap],lw=1,label=f'{geometry} {cap} nm',ms=4)
    for ax,title,ylabel in zip(aa,['(a) Collected complex vector','(b) Zero-order intensity','(c) Largest coefficient change'],
        ['Relative L2 change / %','|delta I0| / final I0 / %','max |delta r_m|']):
        ax.set(title=title,xlabel='Retained maximum order N',ylabel=ylabel,xticks=[240,320]);ax.grid(alpha=.2,which='both')
    handles,labels=aa[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=4,fontsize=8.5)
    fig.suptitle('Actual RCWA truncation controls (N = 160, 240, 320; harmonics = 2N + 1)',fontsize=11.5)
    save(fig,'03_rcwa_convergence')
    csvwrite(FIG/'03_rcwa_convergence.csv',conv)
    dump(FIG/'03_rcwa_convergence.json',dict(source='rcwa/convergence.csv',orders_compared='m=-18..18',
        changes='N240 minus N160; N320 minus N240. Relative I0 uses final N320 I0.',
        final_max_relative_complex_L2=max(float(r['relative_complex_L2']) for r in conv if r['max_order']=='320'),
        final_max_abs_complex_change=max(float(r['max_complex_abs_change']) for r in conv if r['max_order']=='320')))
    np.savez_compressed(HERE/'d4c_all_targets.npz',variant=np.array(list(curves)),s_um=s,
        S_SUM=np.array([v['S_SUM'] for v in curves.values()]),
        S_DIFF_same_kernel=np.array([v['S_DIFF_same_kernel'] for v in curves.values()]),
        AC=np.array([v['AC'] for v in curves.values()]))
    print('Three reference figures, source tables, metrics and aggregate targets written.')


if __name__=='__main__':main()
