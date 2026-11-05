"""Show the strictly positive D4C leakage floors without full-scale masking."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from zero_order_models import K_VALUES

ROOT=Path(__file__).resolve().parents[1]
def main():
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf8'))
    fig,axes=plt.subplots(4,2,figsize=(12,11.5),layout='constrained')
    rows=[];minima=[]
    for ax,v in zip(axes.T.flat,plan['variants']):
        name=v['id']
        with np.load(ROOT/'rcwa/cases'/name/'d4c_targets.npz') as z:
            m=z['order_positive'];tm=z['T_m_working_definition'];phi=z['phi_m_rad']
            zs=float(z['zero_SUM']);zd=float(z['zero_DIFF'])
            stored_s=float(z['S_SUM'].min());stored_d=float(z['S_DIFF_same_kernel'].min())
        weights=np.array([K_VALUES[int(j)] for j in m])
        curvature=float(np.sum(weights*tm*(4*np.pi*m/16)**2))
        half=min(.1,np.sqrt(1.5*zd/max(curvature,1e-30)))
        s=np.linspace(-half,half,401)
        ac=2*np.sum((weights*tm)[:,None]*(1-np.cos(4*np.pi*m[:,None]*s[None,:]/16+phi[:,None])),axis=0)
        isum=ac+zs;diff=ac+zd
        assert zs>0 and zd>0 and abs(stored_s-zs)<1e-14 and abs(stored_d-zd)<1e-14
        ax.plot(s*1000,isum,color='#205b89',label='SUM: AC + residual + leakage')
        ax.plot(s*1000,diff,color='#d77d27',label='DIFF leakage: same supplied AC kernel')
        ax.axhline(zs,color='#205b89',ls=':',lw=.8);ax.axhline(zd,color='#d77d27',ls=':',lw=.8)
        ax.scatter([0,0],[zs,zd],color=['#205b89','#d77d27'],s=25,zorder=5)
        ax.set(title=name.replace('_cap',' / cap ')+' nm',xlabel='Local displacement s / nm',ylabel='Absolute normalized signal',ylim=(0,1.15*max(isum.max(),diff.max())))
        ax.text(.03,.97,f'min SUM = {zs:.6g}\nmin DIFF = {zd:.6g}',transform=ax.transAxes,va='top',fontsize=9,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2));ax.grid(alpha=.2)
        for ss,a,b in zip(s,isum,diff):rows.append(dict(id=name,s_um=float(ss),SUM=float(a),DIFF_same_kernel=float(b),zero_SUM=zs,zero_DIFF=zd))
        minima.append(dict(id=name,SUM_min=stored_s,DIFF_same_kernel_min=stored_d,zero_SUM=zs,zero_DIFF=zd,
            all_minima_positive=True,local_halfwidth_nm=half*1000))
    axes[0,0].legend(loc='lower left',fontsize=7)
    fig.suptitle('D4C minima are strictly positive: independent local scales show the retained leakage floors',fontsize=13)
    for ext in ('png','svg'):fig.savefig(ROOT/'figures'/('16_d4c_positive_minima.'+ext),dpi=190,bbox_inches='tight')
    plt.close(fig)
    out=ROOT/'figures/data';out.mkdir(exist_ok=True)
    for path,data in [('16_d4c_positive_minima.csv',rows),('16_d4c_minima_table.csv',minima)]:
        with (out/path).open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    (out/'16_d4c_positive_minima.json').write_text(json.dumps(dict(
        source='Same saved RCWA amplitudes and confirmed D4C formula, evaluated on finer local displacement grids',
        new_electromagnetic_solves=False,axis='Individual local x/y scales; every y axis starts at zero',
        purpose='Resolve whether leakage was omitted or merely hidden by full-scale linear plotting',
        conclusion='All 16 channel minima equal their strictly positive residual-plus-leakage term',minima=minima),indent=2),encoding='utf8')
    print(json.dumps(minima,indent=2))

if __name__=='__main__':main()
