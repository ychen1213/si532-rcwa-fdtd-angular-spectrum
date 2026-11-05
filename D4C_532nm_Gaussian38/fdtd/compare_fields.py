"""Independent saved-field time/grid comparison for the NA=.6 B scan.

Uses raw E minus measured air leakage; no RCWA fit or peak normalization.
The scan functions are the project's separately maintained formula contract.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'analysis'))
from prepare_spectra import read_monitor, interp_complex, Z0, K0, REF
from zero_order_models import b_features,candidate_scan

def field(case,air):
    R=read_monitor(ROOT/'cases'/case,'R');E=read_monitor(ROOT/'cases'/air,'R');I=read_monitor(ROOT/'cases'/air,'INC')
    k=R['k'];g=R['gamma'];back=np.exp(1j*K0*g*(REF-R['y']))
    A=(R['E']-E['E'])*back
    return k,A,I['metadata']['down_power_W_m']

def compare(baseline,test,baseair='air_fine',testair=None):
    testair=testair or baseair;k,a,pin=field(baseline,baseair);kt,at,pint=field(test,testair)
    if len(kt)!=len(k) or not np.allclose(kt,k,rtol=1e-11,atol=1e-3):at=interp_complex(kt,at,k)
    # Compare incident-normalized fields including their measured source phase.
    at=at*np.sqrt(pin/pint)
    g=np.sqrt(1-(k/K0)**2);dk=float(k[1]-k[0]);w=g*dk/(4*np.pi*Z0*pin);na=abs(k)<.6*K0;G=2*np.pi/16e-6
    norm=lambda v,sel:float(np.sum(w[sel]*abs(v[sel])**2))
    rel=lambda sel:float(np.sqrt(norm(at-a,sel)/norm(a,sel)))
    band=np.sign(k)*np.floor(abs(k/G)+.5+1e-12);nz=na&(band!=0)
    fs=[]
    for v in (a,at):fs.append(b_features(k,v,pin,selected_orders=np.arange(1,19)))
    s=np.linspace(-4,4,1601);scans={};scanmax=0.
    for ch in ('SUM','DIFF'):
        y0=candidate_scan(fs[0],s,channel=ch,convention='K_over_gamma');yt=candidate_scan(fs[1],s,channel=ch,convention='K_over_gamma')
        err=float(max(abs(yt-y0))/np.ptp(y0));scanmax=max(scanmax,err)
        scans[ch]=dict(relative_L2=float(np.linalg.norm(yt-y0)/np.linalg.norm(y0)),Linf_over_baseline_span=err,baseline_min=float(y0.min()),baseline_max=float(y0.max()),test_min=float(yt.min()),test_max=float(yt.max()),raw_max_absolute=float(max(abs(yt-y0))))
    out=dict(baseline=baseline,test=test,base_air_reference=baseair,test_air_reference=testair,comparison='source-leakage-subtracted raw E spectra, measured Pin, common reference .85um, no phase/gain fit',field_L2_NA06=rel(na),field_L2_nonzero_bands_NA06=rel(nz),baseline_power_NA06=norm(a,na),test_power_NA06=norm(at,na),nonzero_power_relative_change=norm(at,nz)/norm(a,nz)-1,zero_band_power_relative_change=norm(at,na&(band==0))/norm(a,na&(band==0))-1,order_bands=[dict(m=m,field_relative_L2=rel(na&(abs(band)==m)),power_relative_change=norm(at,na&(abs(band)==m))/norm(a,na&(abs(band)==m))-1) for m in range(0,19)],scans=scans,criterion=dict(max_field_relative=.001,max_scan_Linf_over_span=.002,passed_total_field_and_scan=rel(na)<.001 and scanmax<.002,passed_also_nonzero_field=rel(na)<.001 and rel(nz)<.001 and scanmax<.002))
    path=ROOT/('comparison_'+test+'_vs_'+baseline+'.json');path.write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2));return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('baseline');ap.add_argument('test');ap.add_argument('--base-air',default='air_fine');ap.add_argument('--test-air');args=ap.parse_args();compare(args.baseline,args.test,args.base_air,args.test_air)
