"""Raw FDTD fields -> common electric angular spectra and source-only modes.

No RCWA target is used in this conversion or in the modal field projection.
The scan calculation remains entirely in angular space.
"""
from __future__ import annotations
from pathlib import Path
import argparse, hashlib, json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
Z0=376.730313668
WAVE=532e-9
K0=2*np.pi/WAVE
REF=.85e-6

def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf8')

def read_monitor(folder,name):
    path=Path(folder)/f'nearfield_{name}.npz'
    with np.load(path) as z:
        x=z['x_m'].squeeze();ez=z['Ez'].squeeze();hx=z['Hx'].squeeze()
        y=float(z['monitor_y_m'])
    dx=float(np.median(np.diff(x)))
    # Same physical FFT span at both meshes. Padding refines numerical
    # quadrature; it does not increase the physical monitor aperture.
    nfft=2**int(round(np.log2(.00524288/dx)))
    assert nfft>=len(x)
    k=2*np.pi*np.fft.fftshift(np.fft.fftfreq(nfft,dx))
    phase=np.exp(-1j*k*x[0]);ends=np.ones(len(x));ends[[0,-1]]=.5
    e=dx*np.fft.fftshift(np.fft.fft(ends*ez,nfft))*phase
    h=dx*np.fft.fftshift(np.fft.fft(ends*hx,nfft))*phase
    good=np.abs(k)<.98*K0;k=k[good];e=e[good];h=h[good]
    g=np.sqrt(1-(k/K0)**2);dk=float(k[1]-k[0])
    up=.5*(e+Z0*h/g);down=.5*(e-Z0*h/g)
    power=lambda a:float(np.sum(g*np.abs(a)**2)*dk/(4*np.pi*Z0))
    return dict(k=k,E=e,H=h,up=up,down=down,y=y,gamma=g,
        metadata=dict(file=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
          dx_m=dx,nfft=nfft,native_points=len(x),native_x_min_m=float(x[0]),native_x_max_m=float(x[-1]),
          monitor_y_m=y,up_power_W_m=power(up),down_power_W_m=power(down),
          direct_net_power_W_m=float(.5*np.real(np.trapezoid(ez*np.conj(hx),x)))))

def interp_complex(x,a,target):
    return np.interp(target,x,a.real,left=0,right=0)+1j*np.interp(target,x,a.imag,left=0,right=0)

def source_projection(k,A,Ainc0,max_order=28,degree=2):
    """Joint shifted-source least squares, not fitted to any scan curve.

    Each order permits a quadratic response in its incident q; the source
    templates and common propagation are measured in the empty FDTD run.
    This is an approximate field decomposition, not a claim that finite
    angular bands uniquely identify Maxwell diffraction orders.
    """
    G=2*np.pi/16e-6;g=np.sqrt(1-(k/K0)**2)
    use=np.abs(k)<.94*K0
    qscale=1/19e-6
    orders=np.arange(-max_order,max_order+1)
    columns=[]
    for m in orders:
        q=k-m*G
        source=interp_complex(k,Ainc0,q)
        # Suppress machine-level/background spectral tails in polynomial
        # templates, without truncating appreciable Gaussian source power.
        source=source*(np.abs(q)<1.0*G)
        base=source*np.exp(1j*K0*g*REF)
        for d in range(degree+1):columns.append(base*(q/qscale)**d)
    design=np.column_stack(columns)
    weight=np.sqrt(g[use])
    coeff,_,rank,sing=np.linalg.lstsq(design[use]*weight[:,None],A[use]*weight,rcond=1e-10)
    pred=design@coeff
    j=int(np.flatnonzero(orders==0)[0])*(degree+1)
    zero=design[:,j:j+degree+1]@coeff[j:j+degree+1]
    norm=lambda a:float(np.sum(g[use]*np.abs(a[use])**2))
    return zero,dict(degree=degree,max_order=max_order,rank=int(rank),columns=len(coeff),
        condition=float(sing[0]/sing[-1]),relative_field_L2=float(np.sqrt(norm(A-pred)/norm(A))),
        fitted_r0_at_q0_re=float(coeff[j].real),fitted_r0_at_q0_im=float(coeff[j].imag),
        fitted_coefficients_re=coeff.real.reshape(-1,degree+1).tolist(),
        fitted_coefficients_im=coeff.imag.reshape(-1,degree+1).tolist(),
        order=orders.tolist(),incident_q_scale_rad_m=qscale,
        template_q_cutoff_rad_m=1.0*G,source='empty-domain FDTD incident field only; no RCWA data',
        limitation='finite-band polynomial source-template decomposition, not an exact modal observable')

def prepare(case,suffix='',project=True):
    folder=ROOT/'fdtd/cases'/(case+suffix)
    airref=folder/'AIR_REFERENCE.json'
    air=ROOT/'fdtd/cases'/(json.loads(airref.read_text(encoding='utf8'))['air_reference'] if airref.exists() else 'air'+suffix)
    R=read_monitor(folder,'R');R2=read_monitor(folder,'R2')
    E=read_monitor(air,'R');E2=read_monitor(air,'R2');I=read_monitor(air,'INC')
    k=R['k'];g=R['gamma'];dk=float(k[1]-k[0])
    for data in (R2,E,E2,I):assert np.allclose(k,data['k'],rtol=1e-12,atol=1e-4)
    # Main extraction: total E above source minus the identical empty-run E.
    # Avoid continuum E/H directional bias on a discrete Yee grid.
    A=(R['E']-E['E'])*np.exp(1j*K0*g*(REF-R['y']))
    A2=(R2['E']-E2['E'])*np.exp(1j*K0*g*(REF-R2['y']))
    AEH=(R['up']-E['up'])*np.exp(1j*K0*g*(REF-R['y']))
    Ain=I['down']*np.exp(-1j*K0*g*(REF-I['y']))
    Ain0=Ain*np.exp(1j*K0*g*REF)
    Pin=I['metadata']['down_power_W_m']
    power=lambda a,NA=.98:float(np.sum(g*np.abs(a)**2*(np.abs(k)<NA*K0))*dk/(4*np.pi*Z0))
    rel=lambda a,b:float(np.sqrt(power(a-b,.6)/power(b,.6)))
    out=ROOT/'outputs/spectra';out.mkdir(parents=True,exist_ok=True)
    arrays=dict(kx_rad_m=k,A_ref_V=A,A_ref_R2_V=A2,A_ref_EH_V=AEH,Ainc_ref_V=Ain,Ainc_y0_V=Ain0,
        Pin_W_m=Pin,wavelength_m=WAVE,reference_y_m=REF)
    meta=dict(case=case,suffix=suffix,method='empty-R-subtracted total E; upward continuum reference propagation',
        fourier='A(k)=integral E(x) exp(-ikx) dx',time='exp(-i omega t)',
        Pin_W_m=Pin,reflected_NA06_over_Pin=power(A,.6)/Pin,reflected_098_over_Pin=power(A)/Pin,
        R2_vs_R_field_L2_NA06=rel(A2,A),EH_vs_Eonly_field_L2_NA06=rel(AEH,A),
        Pin_definition='empty INC down-propagating E/H angular flux; no fitted gain',
        monitors={name:data['metadata'] for name,data in [('caseR',R),('caseR2',R2),('airR',E),('airINC',I)]})
    if project:
        for deg in (0,2):
            az,pm=source_projection(k,A,Ain0,degree=deg)
            arrays[f'A_zero_degree{deg}_V']=az
            pm['zero_power_over_Pin']=power(az,.6)/Pin
            meta[f'projection_degree{deg}']=pm
    np.savez_compressed(out/(case+suffix+'.npz'),**arrays)
    dump(out/(case+suffix+'.json'),meta)
    print(json.dumps({key:val for key,val in meta.items() if key!='monitors' and not key.startswith('projection')},ensure_ascii=False),flush=True)
    return arrays,meta

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--suffix',default='');p.add_argument('--no-project',action='store_true')
    args=p.parse_args();prepare(args.case,args.suffix,not args.no_project)
