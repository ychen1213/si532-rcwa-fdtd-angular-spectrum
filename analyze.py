"""Reproducible complex angular-spectrum analysis of genuine saved FDTD fields.

No electromagnetic fields are generated here. Analytic window-only examples
are labelled synthetic and are not substitutes for the independently run FDTD.
"""
from pathlib import Path
import json, csv, math, re, textwrap
import numpy as np
from scipy.special import sici
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'analysis'; FIG=ROOT/'figures'
DATA.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
P=2.; WL=.532; G=2*np.pi/P; K0=2*np.pi/WL; M=np.arange(-3,4)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'figure.dpi':150,
    'axes.grid':True,'grid.alpha':.22,'savefig.bbox':'tight','axes.spines.top':False,
    'axes.spines.right':False,'lines.linewidth':1.5,'figure.constrained_layout.use':True,
    'axes.titlesize':11,'axes.labelsize':11})

def savefig(name):
    fig=plt.gcf()
    w,h=fig.get_size_inches()
    if w>9:fig.set_size_inches(9,h*9/w)
    for ax in fig.axes:
        if len(ax.get_title())>45:ax.set_title('\n'.join(textwrap.wrap(ax.get_title(),43)))
    for obj in fig.findobj(matplotlib.text.Text):
        if obj.get_fontsize()<11:obj.set_fontsize(11)
    plt.savefig(FIG/(name+'.png'),dpi=180)
    plt.savefig(FIG/(name+'.svg'));plt.close()

def writecsv(path,rows):
    if not rows:return
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def dump(path,obj):
    def default(o):
        if isinstance(o,np.ndarray):return o.tolist()
        if isinstance(o,np.generic):return o.item()
        if isinstance(o,complex):return {'real':o.real,'imag':o.imag}
        raise TypeError(type(o).__name__)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=default),encoding='utf-8')

def loadfield(path):
    with np.load(path) as z:
        return {k:np.asarray(z[k]).squeeze() for k in ['x_m','y_m','frequency_hz','Ez','Hx']}

def resolve_air(folder):
    ref=folder/'AIR_REFERENCE.json'
    if not ref.exists():return folder/'air',False
    info=json.loads(ref.read_text(encoding='utf-8'))
    base=ROOT/'fdtd'/'cases'
    path=(base/info.get('path_relative_to_cases',info.get('air_reference'))).resolve()
    assert path.is_relative_to(base.resolve()),'Air reference must be an archived case in this experiment.'
    return path,True

def crop_uniform(f,L,dx):
    # Exact integer-period endpoints; trapezoidal endpoint weighting avoids a
    # spurious extra pixel of aperture. Units here and below are micrometres.
    x=np.linspace(-L/2,L/2,int(round(L/dx))+1)
    xx=f['x_m']*1e6
    assert x[0]>=xx[0] and x[-1]<=xx[-1]
    interp=lambda a:np.interp(x,xx,a.real)+1j*np.interp(x,xx,a.imag)
    return x,interp(f['Ez']),interp(f['Hx'])

def spectrum(x,e,padding=16):
    dx=x[1]-x[0]; nfft=2**int(np.ceil(np.log2(len(x)*padding)))
    weighted=e.copy();weighted[[0,-1]]*=.5
    k=np.fft.fftshift(np.fft.fftfreq(nfft,dx))*2*np.pi
    A=dx*np.fft.fftshift(np.fft.fft(weighted,nfft))*np.exp(-1j*k*x[0])
    return k,A

def project(x,e,k):
    return np.trapezoid(e[None,:]*np.exp(-1j*k[:,None]*x[None,:]),x,axis=1)

def theory():
    rows=[]
    for n in [10,20,50,100,200,400]:
        opt=minimize_scalar(lambda t:-abs(np.sinc(n+t)),bounds=(-.999,-.001),method='bounded')
        leak=-opt.fun
        theta1=np.arcsin(WL/P)
        width=np.degrees(np.arcsin(WL/P+WL/(n*P))-np.arcsin(WL/P-WL/(n*P)))
        # integral sinc^2(pi v) from -b..b.
        def sinc_power(b):
            return float(2/np.pi*(sici(2*np.pi*b)[0]+(np.cos(2*np.pi*b)-1)/(2*np.pi*b)))
        rows.append(dict(N=n,L_um=n*P,mainlobe_width_deg=width,zero_order_at_exact_first=0.,
          zero_order_leak_in_first_mainlobe=leak,leak_dB=20*np.log10(leak),
          conservative_amplitude_bound=1/(np.pi*(n-1)),mainlobe_scalar_energy_fraction=sinc_power(1),
          nearest_order_band_scalar_energy_fraction=sinc_power(n/2)))
    writecsv(DATA/'analytic_window_metrics.csv',rows)
    mainpower=rows[0]['mainlobe_scalar_energy_fraction'];mainfield=2*sici(np.pi)[0]/np.pi
    fig,axs=plt.subplots(1,2,figsize=(11,3.8))
    u=np.linspace(.8,1.2,16001)
    for n in [10,20,50,100,200,400]:
        axs[0].plot(u,abs(np.sinc(n*(u-1))),label=f'N={n}')
    axs[0].set(xlabel=r'$u=k_x/G$',ylabel='Unit-amplitude isolated order',title='More periods narrow the peak');axs[0].legend(ncol=3,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.25))
    axs[1].loglog([r['N'] for r in rows],[r['zero_order_leak_in_first_mainlobe'] for r in rows],'o-',label='Max zero-order amplitude near order +1')
    axs[1].loglog([r['N'] for r in rows],[r['conservative_amplitude_bound'] for r in rows],'--',label=r'$1/[\pi(N-1)]$')
    axs[1].set(xlabel='Number of periods N',ylabel='Leak / zero-order peak',title='Envelope decreases; exact order center is zero');axs[1].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.25))
    savefig('03_sinc_and_leakage')
    fig,axs=plt.subplots(1,2,figsize=(10,3.7))
    v=np.linspace(-5,5,10001);f=np.sinc(v)
    axs[0].plot(v,f,label=r'$\mathrm{sinc}(\pi v)$');axs[0].fill_between(v,0,f,where=abs(v)<=1,alpha=.25)
    axs[0].set(xlabel=r'$v=(k_x-k_m)L/(2\pi)$',ylabel='Signed complex-field kernel',title=f'Main-lobe field integral = {mainfield:.6f}')
    axs[1].plot(v,f*f);axs[1].fill_between(v,0,f*f,where=abs(v)<=1,alpha=.25)
    axs[1].set(xlabel=r'$v=(k_x-k_m)L/(2\pi)$',ylabel='Scalar spectral-energy kernel',title=f'Scalar main-lobe fraction = {mainpower:.6f}')
    savefig('04_field_vs_power_integrals')
    fig,ax=plt.subplots(figsize=(10,3.8))
    for j in range(-4,4):
        ax.add_patch(patches.Rectangle((j*P+.5,-.1),1,.1,color='#476882'))
    ax.add_patch(patches.Rectangle((-7.5,-.5),15,.4,color='#476882'))
    ax.axhline(0,color='k',lw=.7);ax.axhline(1.,ls='--',color='#cf7936',label='R monitor: y=+1.0 um')
    ax.axhline(.7,ls='--',color='#4b8d60',label='Plane-wave source: y=+0.7 um')
    ax.annotate('Incident TE',xy=(0,.1),xytext=(0,.6),ha='center',arrowprops={'arrowstyle':'->'})
    ax.text(-6.5,-.35,'Si substrate',color='white');ax.text(-6.6,.18,'Air; native Ez = scalar TE field')
    ax.annotate('',xy=(-2,-.6),xytext=(0,-.6),arrowprops={'arrowstyle':'<->'});ax.text(-1,-.75,'pitch = 2 um',ha='center')
    ax.set(xlim=(-7.5,7.5),ylim=(-.9,1.5),xlabel='Transverse coordinate (um); schematic only',ylabel='Native y (um)',title='50% Si relief grating; height = 100 nm')
    ax.legend(loc='upper right',fontsize=8);savefig('01_geometry')
    return rows,mainpower,mainfield

def read_rcwa():
    # Agent exports a stable seven-row comparison table. This adapter is kept
    # deliberately explicit so that phase conventions remain reviewable.
    p=ROOT/'rcwa'/'reference_for_analysis.npz'
    z=np.load(p)
    assert np.array_equal(z['m'],M)
    return z['r_native_origin'].astype(complex)

def analyze_case(folder,rref):
    air,reused=resolve_air(folder);mark=folder/'mark'
    ad=json.loads((air/'DONE.json').read_text(encoding='utf-8'))
    md=json.loads((mark/'DONE.json').read_text(encoding='utf-8'));cfg=md['settings']
    observed_path=mark/'observed_settings.json'
    observed=json.loads(observed_path.read_text(encoding='utf-8')) if observed_path.exists() else {}
    N=cfg['N'];L=N*P;dx=cfg['dx_m']*1e6
    R=loadfield(mark/'fields_R.npz');Ra=loadfield(air/'fields_R.npz');I=loadfield(air/'fields_I.npz')
    x,e,h=crop_uniform(R,L,dx);_,ea,ha=crop_uniform(Ra,L,dx);_,g,gh=crop_uniform(I,L,dx)
    e=e-ea;h=h-ha
    g0=np.trapezoid(g,x)/L
    ym=float(R['y_m'])*1e6;kz=np.sqrt(K0*K0-(M*G)**2)
    Am=project(x,e,M*G)
    r=Am*np.exp(-1j*kz*ym)/(L*g0)
    phase=np.angle(r/rref,deg=True); rel=abs(r-rref)/abs(rref)
    # Compare two monitors without changing aperture length or normalization.
    R2=loadfield(mark/'fields_R2.npz');R2a=loadfield(air/'fields_R2.npz')
    _,e2,_=crop_uniform(R2,L,dx);_,e2a,_=crop_uniform(R2a,L,dx)
    r2=project(x,e2-e2a,M*G)*np.exp(-1j*kz*float(R2['y_m'])*1e6)/(L*g0)
    k,A=spectrum(x,e);keep=abs(k)<K0;k=k[keep];A=A[keep]
    u=k/G;a=A/(L*g0);gamma=np.sqrt(np.maximum(0,1-(k/K0)**2))
    deprop=a*np.exp(-1j*np.sqrt(np.maximum(0,K0*K0-k*k))*ym)
    # Synthetic windowed RCWA prediction is at the SAME monitor before FFT.
    emodel=g0*np.sum(rref[:,None]*np.exp(1j*(M*G)[:,None]*x[None,:]+1j*kz[:,None]*ym),axis=0)
    ks,As=spectrum(x,emodel);assert np.allclose(ks[keep],k)
    amodel=As[keep]/(L*g0)
    rows=[]
    Hm=project(x,h,M*G)
    for i,m in enumerate(M):
        band=abs(u-m)<=.5;main=abs(u-m)<=1/N
        # R= integral gamma |A|^2 dk /(2 pi L |g0|^2) in air.
        bandpow=np.trapezoid(gamma[band]*abs(A[band])**2,k[band])/(2*np.pi*L*abs(g0)**2)
        mainpow=np.trapezoid(gamma[main]*abs(A[main])**2,k[main])/(2*np.pi*L*abs(g0)**2)
        refpow=np.sqrt(1-(m*G/K0)**2)*abs(rref[i])**2
        rows.append(dict(case=folder.name,N=N,dx_nm=cfg['dx_m']*1e9,dy_nm=cfg['dy_m']*1e9,
          max_time_fs=cfg['t_max_s']*1e15,m=int(m),theta_deg=np.degrees(np.arcsin(m*G/K0)),
          rcwa_re=rref[i].real,rcwa_im=rref[i].imag,fdtd_re=r[i].real,fdtd_im=r[i].imag,
          rcwa_abs=abs(rref[i]),fdtd_abs=abs(r[i]),complex_absolute_error=abs(r[i]-rref[i]),
          complex_relative_error=rel[i],phase_error_deg=phase[i],rcwa_efficiency=refpow,
          peak_estimated_efficiency=np.sqrt(1-(m*G/K0)**2)*abs(r[i])**2,
          nearest_order_band_efficiency=bandpow,mainlobe_efficiency=mainpow,
          R2_re=r2[i].real,R2_im=r2[i].imag,R_R2_complex_abs_difference=abs(r[i]-r2[i]),
          EH_upgoing_relative_residual=abs(376.730313668*Hm[i]/np.sqrt(1-(m*G/K0)**2)-Am[i])/abs(Am[i])))
    np.savez_compressed(DATA/(folder.name+'_angular_spectrum.npz'),kx_rad_um=k,u=u,
      A_monitor_field_times_um=A,a_monitor=a,a_ridge_top=deprop,
      windowed_rcwa_a_monitor=amodel,window_x_um=x,reflected_E=e,incident_E=g,
      predicted_periodic_E=emodel,g0=g0,orders=M,r_fdtd=r,r_rcwa=rref,r_fdtd_R2=r2)
    # Near-uniform incident amplitude is measured, never merely assumed.
    ci=abs(x)<L/2-2*P
    incident_nonuniform=np.sqrt(np.mean(abs(g[ci]-g0)**2))/abs(g0) if ci.any() else np.nan
    logfile=mark/'engine_stdout.log'
    if not logfile.exists():logfile=mark/'engine_stdout_public.log'
    log=logfile.read_text(encoding='utf-8',errors='replace') if logfile.exists() else ''
    autos=re.findall(r'Auto Shutoff:\s*([\deE+.-]+)',log)
    fullR=float(np.trapezoid(gamma*abs(A)**2,k)/(2*np.pi*L*abs(g0)**2))
    directR=float(376.730313668*np.trapezoid(np.real(e*np.conj(h)),x)/(L*abs(g0)**2))
    actualIncident=float(-376.730313668*np.trapezoid(np.real(g*np.conj(gh)),x)/(L*abs(g0)**2))
    summary=dict(case=folder.name,N=N,dx_nm=cfg['dx_m']*1e9,dy_nm=cfg['dy_m']*1e9,
      time_fs=cfg['t_max_s']*1e15,pml_layers=int(observed.get('pml layers',cfg['pml_layers'])),
      requested_pml_layers=cfg['pml_layers'],requested_auto_shutoff=cfg['auto_shutoff'],
      mark_solve_s=md['solve_s'],air_solve_s=0. if reused else ad['solve_s'],
      air_reference_reused=reused,air_reference_original_solve_s=ad['solve_s'],
      g0=g0,incident_rms_nonuniformity_excluding_2_period_edges=incident_nonuniform,
      relative_complex_L2=np.linalg.norm(r-rref)/np.linalg.norm(rref),
      maximum_absolute_complex_error=float(max(abs(r-rref))),
      maximum_phase_error_deg=float(max(abs(phase))),
      maximum_monitor_difference=float(max(abs(r-r2))),
      total_rcwa_R=float(np.sum([t['rcwa_efficiency'] for t in rows])),
      total_peak_R=float(np.sum([t['peak_estimated_efficiency'] for t in rows])),
      total_band_R=float(np.sum([t['nearest_order_band_efficiency'] for t in rows])),
      full_propagating_window_R=fullR,
      direct_cropped_Poynting_R=directR,
      spectrum_vs_cropped_flux_relative_difference=fullR/directR-1,
      actual_incident_power_over_meanE_plane_power=actualIncident,
      mark_last_logged_auto_shutoff=float(autos[-1]) if autos else None)
    return rows,summary,(u,a,amodel,x,e,g,g0,r,r2)

def plots(results,rref):
    production=[v for v in results if v[1]['case']==f"N{v[1]['N']:03d}_dx20_dy10_t600"]
    production.sort(key=lambda v:v[1]['N']);best=max(production,key=lambda v:v[1]['N'])
    fig,axs=plt.subplots(1,2,figsize=(10,3.7))
    axs[0].stem(M,abs(rref));axs[0].set(xlabel='Order m',ylabel=r'$|r_m|$',title='RCWA complex electric-field amplitudes')
    axs[1].plot(M,np.angle(rref,deg=True),'o-');axs[1].set(xlabel='Order m',ylabel='Phase at ridge top (deg)',title='Same x origin as finite FDTD')
    savefig('02_rcwa_orders')
    fig,axs=plt.subplots(2,1,figsize=(11,7),sharex=True)
    for rows,s,v in production:
        u,a,*_=v;axs[0].plot(u,abs(a),label=f'N={s["N"]}');axs[1].plot(u,20*np.log10(np.maximum(abs(a),1e-8)),label=f'N={s["N"]}')
    for ax in axs:
        for m in M:ax.axvline(m,color='k',lw=.5,alpha=.25)
    axs[0].set(ylabel=r'$|A(k_x)/(L E_i)|$',title='Genuine finite-array FDTD: normalized complex angular spectrum')
    axs[1].set(xlabel=r'$u=k_x/G$ (angle: $\theta=\arcsin(u\lambda/\Lambda)$)',ylabel='Field magnitude (dB)',ylim=(-85,0))
    fig.legend(*axs[0].get_legend_handles_labels(),loc='outside upper center',ncol=6,frameon=False)
    savefig('05_fdtd_spectra_all_N')
    rows,s,v=best;u,a,amodel,x,e,g,g0,r,r2=v;N=s['N']
    fig,axs=plt.subplots(2,2,figsize=(11,7))
    for m,ax in zip([0,1,2,3],axs.ravel()):
        vloc=(u-m)*N;sel=abs(vloc)<6
        ax.plot(vloc[sel],abs(a[sel]),label='FDTD')
        ax.plot(vloc[sel],abs(amodel[sel]),'--',label='Windowed periodic RCWA')
        ax.set(xlabel=r'$(k_x-k_m)L/(2\pi)$',ylabel='Normalized magnitude',title=f'Order {m}; N={N}')
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside upper center',ncol=2,frameon=False)
    savefig('06_local_sinc_validation')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for i,m in enumerate([0,1,2,3]):
        rr=[next(row for row in t[0] if row['m']==m) for t in production]
        ns=[t['N'] for t in rr]
        axs[0].loglog(ns,[t['complex_absolute_error'] for t in rr],'o-',label=f'm={m}')
        axs[1].semilogx(ns,[t['phase_error_deg'] for t in rr],'o-',label=f'm={m}')
    axs[0].set(xlabel='N',ylabel=r'$|\hat r_m-r_m|$',title='Period convergence at fixed mesh')
    axs[1].set(xlabel='N',ylabel='Phase error (deg)',title='Residuals include discretization and edges')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='outside upper center',ncol=4,frameon=False)
    savefig('07_complex_order_convergence')
    fig,ax=plt.subplots(figsize=(6,5))
    ax.scatter(rref.real,rref.imag,marker='x',s=90,label='RCWA',color='k')
    for rr,ss,vv in production:
        rv=vv[-2];ax.plot(rv.real,rv.imag,'o',ms=4,label=f'N={ss["N"]}')
    for i,m in enumerate(M):
        if m>=0:ax.annotate(f'm={m}',(rref[i].real,rref[i].imag),xytext=(5,5),textcoords='offset points')
    ax.axhline(0,color='k',lw=.5);ax.axvline(0,color='k',lw=.5)
    ax.set(xlabel=r'$\Re(r_m)$',ylabel=r'$\Im(r_m)$',title='Complex-plane comparison; symmetric +/- orders overlap',aspect='equal');ax.legend(fontsize=8)
    savefig('08_complex_plane')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    pos=[row for row in rows if row['m']>=0];xx=np.arange(4);w=.2
    for off,key,label in [(-1.5,'rcwa_efficiency','RCWA'),(-.5,'peak_estimated_efficiency','FDTD center'),(.5,'nearest_order_band_efficiency','FDTD order band'),(1.5,'mainlobe_efficiency','FDTD main lobe')]:
        axs[0].bar(xx+off*w,[row[key] for row in pos],w,label=label)
    axs[0].set(xticks=xx,xticklabels=[0,1,2,3],xlabel='m',ylabel='Efficiency',title=f'Center extraction and power integration; N={N}');axs[0].legend(fontsize=7)
    for rr,ss,vv in production:
        pp=[row for row in rr if row['m']>=0]
        axs[1].plot([t['m'] for t in pp],[t['mainlobe_efficiency']/t['peak_estimated_efficiency'] for t in pp],'o-',label=f'N={ss["N"]}')
    axs[1].axhline(2/np.pi*(sici(2*np.pi)[0]),color='k',ls='--',label='Isolated sinc: 0.902823')
    axs[1].set(xlabel='m',ylabel='Main-lobe / inferred order power',title='Main-lobe fraction: isolated sinc versus actual fields');axs[1].legend(fontsize=7,ncol=3,loc='upper center',bbox_to_anchor=(.5,-.2))
    savefig('09_peak_and_integrated_power')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    axs[0].plot(x,abs(g/g0),label='Amplitude');axs[0].set(xlabel='x (um)',ylabel=r'$|E_i(x)/\overline{E_i}|$',title=f'Actual air-reference illumination, N={N}')
    ax2=axs[0].twinx();ax2.plot(x,np.angle(g/g0,deg=True),color='#d07a34',alpha=.8,label='Phase');ax2.set_ylabel('Incident phase relative to mean (deg)',color='#d07a34')
    axs[1].plot(M,abs(r-r2),'o-');axs[1].set(xlabel='m',ylabel=r'$|r_m(R)-r_m(R2)|$',title='After correcting each monitor propagation phase')
    savefig('10_incident_and_reference_plane_checks')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    axs[0].plot([t[1]['N'] for t in production],[t[1]['mark_solve_s'] for t in production],'o-',label='Mark FDTD')
    axs[0].plot([t[1]['N'] for t in production],[t[1]['air_solve_s'] for t in production],'o-',label='Air reference')
    axs[0].set(xlabel='N',ylabel='Measured wall time (s)',title='Actual solver cost');axs[0].legend(fontsize=8)
    # Every solved case is shown, with its numerical settings in the label.
    for rr,ss,vv in results:
        if ss['N']==20:
            axs[1].plot(M,[row['complex_absolute_error'] for row in rr],'o-',label=f'{ss["dx_nm"]:g}/{ss["dy_nm"]:g} nm; PML {ss["pml_layers"]}')
    axs[1].set(xlabel='m',ylabel='Absolute complex error',title='N=20: mesh and PML controls');axs[1].legend(fontsize=7)
    savefig('11_runtime_and_mesh_time_checks')
    selected=min(results,key=lambda v:v[1]['relative_complex_L2'])
    rr,ss,vv=selected
    rv=vv[-2]
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    axs[0].plot(M,rref.real,'kx',ms=9,label='RCWA real')
    axs[0].plot(M,rv.real,'o',mfc='none',label='FDTD real')
    axs[0].plot(M,rref.imag,'k+',ms=9,label='RCWA imaginary')
    axs[0].plot(M,rv.imag,'s',mfc='none',label='FDTD imaginary')
    axs[0].set(xlabel='m',ylabel='Complex electric-field coefficient',title='Best independently solved case: direct complex comparison')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='outside upper center',ncol=4,frameon=False)
    axs[1].plot(M,[t['complex_relative_error']*100 for t in rr],'o-',label='Complex relative error')
    axs[1].set(xlabel='m',ylabel='Complex relative error (%)',title=f'N={ss["N"]}; dx/dy={ss["dx_nm"]:g}/{ss["dy_nm"]:g} nm')
    savefig('12_best_case_complex_comparison')

def main():
    theory_rows,mainpower,mainfield=theory()
    rref=read_rcwa();results=[]
    for folder in sorted((ROOT/'fdtd'/'cases').iterdir()):
        air,_=resolve_air(folder)
        if (folder/'mark'/'DONE.json').exists() and (air/'DONE.json').exists():
            results.append(analyze_case(folder,rref))
    if not results:raise RuntimeError('No completed air/mark pair; do not fabricate FDTD results.')
    writecsv(DATA/'all_extracted_orders.csv',[r for v in results for r in v[0]])
    dump(DATA/'case_summary.json',[v[1] for v in results])
    writecsv(DATA/'case_summary.csv',[{k:v for k,v in t[1].items() if k!='g0'} for t in results])
    dump(DATA/'summary.json',{'theory':theory_rows,'mainlobe_scalar_energy_fraction':mainpower,
      'mainlobe_field_integral_factor':mainfield,'rcwa_r_native_origin':rref,
      'cases':[v[1] for v in results],
      'provenance':'Real Lumerical fields in fdtd/cases; analytic window examples separately identified.',
      'normalization':'A=integral E exp(-ikx) dx; units of x are um; divide by L(um) and complex incident mean; exp(-iky*y_R) propagation correction.'})
    plots(results,rref)
    print(json.dumps({'completed_case_pairs':len(results),'cases':[v[1]['case'] for v in results]},indent=2))

if __name__=='__main__':main()
