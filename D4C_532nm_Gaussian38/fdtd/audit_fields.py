"""Read-only field and absolute-flux audits; writes audit JSON, never a solver."""
from pathlib import Path
import argparse,json
import numpy as np

ROOT=Path(__file__).resolve().parent
Z0=376.730313668
K0=2*np.pi/532e-9
REF=.85e-6

def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def angular(folder,name):
    with np.load(ROOT/'cases'/folder/f'nearfield_{name}.npz') as z:
        x=z['x_m'];e=z['Ez'];h=z['Hx'];y=float(z['monitor_y_m']);dx=float(np.median(np.diff(x)));n=2**int(np.ceil(np.log2(len(x)*16)))
    weights=np.ones(len(x));weights[[0,-1]]=.5
    k=2*np.pi*np.fft.fftshift(np.fft.fftfreq(n,dx));phase=np.exp(-1j*k*x[0]);a=dx*np.fft.fftshift(np.fft.fft(weights*e,n))*phase;b=dx*np.fft.fftshift(np.fft.fft(weights*h,n))*phase
    good=abs(k)<.98*K0;k=k[good];a=a[good];b=b[good];gamma=np.sqrt(1-(k/K0)**2);dk=float(k[1]-k[0]);up=.5*(a+Z0*b/gamma);down=.5*(a-Z0*b/gamma)
    flux=lambda v:float(np.sum(gamma*abs(v)**2)*dk/(4*np.pi*Z0))
    return dict(k=k,gamma=gamma,up=up*np.exp(1j*K0*gamma*(REF-y)),down=down*np.exp(-1j*K0*gamma*(REF-y)),up_power=flux(up),down_power=flux(down),net_power=flux(up)-flux(down),direct_net_power=float(.5*np.real(np.trapezoid(e*np.conj(h),x))),dx=dx,y=y)

def audit(suffix='',mark_case='ah53_cap000'):
    air='air'+suffix;flat='flatSi'+suffix;ah=mark_case+suffix
    I=angular(air,'INC');R=angular(air,'R');R2=angular(air,'R2');F=angular(flat,'R');A=angular(ah,'R');A2=angular(ah,'R2')
    pin=I['down_power'];q=I['k'];sel=abs(I['down'])>max(abs(I['down']))*1e-3
    # Backpropagate reflection/incident ratio to y=0, then fit only scalar gain.
    rfresnel=(1-complex(4.149971915645845,.048064634885764825))/(1+complex(4.149971915645845,.048064634885764825))
    f0=F['up']*np.exp(-2j*K0*I['gamma']*REF);inc=I['down'];fit=np.vdot(inc[sel],f0[sel])/np.vdot(inc[sel],inc[sel])
    vals=dict(field_definition='raw six-component monitors; angular E/H directional split; trapezoid endpoints; |k|<.98k0; reference_y .85um',suffix=suffix,Pin_air_down_W_per_m=pin,air_R_up_over_Pin=R['up_power']/pin,air_R_down_over_Pin=R['down_power']/pin,air_R2_up_over_Pin=R2['up_power']/pin,flat_R_over_Pin=F['up_power']/pin,flat_fresnel_R=abs(rfresnel)**2,flat_R_relative_error=F['up_power']/pin/abs(rfresnel)**2-1,flat_fitted_r_y0=dict(re=float(fit.real),im=float(fit.imag)),flat_fresnel_r=dict(re=rfresnel.real,im=rfresnel.imag),flat_complex_relative_error=float(abs(fit-rfresnel)/abs(rfresnel)),ah53_R_over_Pin=A['up_power']/pin,ah53_R_down_over_up=A['down_power']/A['up_power'],ah53_R2_vs_R_complex_relative_L2=float(np.sqrt(np.sum(A['gamma']*abs(A2['up']-A['up'])**2)/np.sum(A['gamma']*abs(A['up'])**2))),direct_vs_spectral_flux={name:dict(direct_W_per_m=v['direct_net_power'],spectral_W_per_m=v['net_power'],relative_difference=(v['net_power']-v['direct_net_power'])/v['direct_net_power']) for name,v in [('airINC',I),('airR',R),('flatR',F),('ah53R',A)]})
    vals['mark_case']=mark_case
    dump(ROOT/('baseline_flux_audit'+suffix+'.json'),vals);print(json.dumps(vals,indent=2));return vals

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--suffix',default='');ap.add_argument('--mark-case',default='ah53_cap000');args=ap.parse_args();audit(args.suffix,args.mark_case)
