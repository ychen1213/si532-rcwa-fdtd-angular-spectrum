"""Real Ansys Lumerical RCWA reference, single-frequency TE silicon relief mark.

Run with Lumerical's bundled Python or Python with numpy and lumapi installed.
Set LUMERICAL_ROOT if the installation differs from the default below.
Each launch creates a fresh hidden CAD session and closes only that session.
"""
from __future__ import annotations
import csv, json, os, sys, time
from pathlib import Path
from contextlib import contextmanager
import numpy as np

ROOT=Path(__file__).resolve().parent
INSTALL=Path(os.environ.get('LUMERICAL_ROOT',r'D:\Ansys Inc\v252\Lumerical'))
sys.path.insert(0,str(INSTALL/'api'/'python'))
WAVELENGTH=532e-9
PERIOD=2e-6
DEPTH=100e-9
REFERENCE=-200e-9
DATASETS=['grating_characterization','grating_power','grating_orders','field_amplitudes','total_energy','substrate','index']

def json_default(value):
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,np.generic): return value.item()
    if isinstance(value,complex): return {'real':value.real,'imag':value.imag}
    raise TypeError(type(value).__name__)

def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,default=json_default),encoding='utf-8')

def archive_dataset(dataset,out,name):
    arrays={k:v for k,v in dataset.items() if isinstance(v,np.ndarray)}
    metadata={k:v for k,v in dataset.items() if not isinstance(v,np.ndarray)}
    metadata['arrays']={k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in arrays.items()}
    np.savez_compressed(out/(name+'.npz'),**arrays)
    write_json(out/(name+'.json'),metadata)

@contextmanager
def session():
    import lumapi
    with lumapi.FDTD(hide=True,serverArgs={'threads':'4'}) as s:
        row=1
        props=['job launching preset','device type','processes','threads','capacity','active']
        old={p:s.getresource('FDTD',row,p) for p in props}
        try:
            for p,v in [('job launching preset','Local Computer'),('device type','CPU'),('processes','1'),('threads','4'),('capacity','1'),('active','1')]:
                s.setresource('FDTD',row,p,v)
            yield s
        finally:
            for p in props:
                try:s.setresource('FDTD',row,p,old[p])
                except Exception:pass

def lsf(nk,harmonics,dx_nm,flat=False):
    def q(v):
        if isinstance(v,str):return json.dumps(v)
        if isinstance(v,bool):return '1' if v else '0'
        if isinstance(v,(tuple,list)):
            if isinstance(v[0],(tuple,list)):return '['+';'.join(','.join(q(x) for x in row) for row in v)+']'
            return '['+','.join(q(x) for x in v)+']'
        return f'{float(v):.17g}'
    lines=['# SI units. canonical z increases into substrate; exp(-i omega t).','newproject;',
           'material_id=addmaterial("(n,k) Material");',
           'setmaterial(material_id,"name","Si532_fitted_fixed");',
           f'setmaterial("Si532_fitted_fixed","Refractive Index",{nk.real:.17g});',
           f'setmaterial("Si532_fitted_fixed","Imaginary Refractive Index",{nk.imag:.17g});']
    def props(d):
        for k,v in d.items():lines.append(f'set({q(k)},{q(v)});')
    lines.append('addrect;')
    props({'name':'substrate','x min':-2e-6,'x max':2e-6,'y min':-1e-6,'y max':1e-6,
           'z min':0 if flat else DEPTH,'z max':600e-9,'material':'Si532_fitted_fixed'})
    if not flat:
        lines.append('addrect;')
        props({'name':'ridge','x min':-0.5e-6,'x max':0.5e-6,'y min':-1e-6,'y max':1e-6,
               'z min':0,'z max':DEPTH,'material':'Si532_fitted_fixed'})
    lines.append('addrcwa;')
    props({'name':'RCWA','simulation region':'2D Y-normal','propagation axis':'z','propagation direction':'forward',
           'x min':-PERIOD/2,'x max':PERIOD/2,'y':0,'z min':-300e-9,'z max':400e-9,
           'background material':'<Object defined dielectric>','index':1,
           'interface position':'absolute','interface absolute positions':[-200e-9,0,100e-9,200e-9],
           'max number k vectors':harmonics,'use li factorization':False,'mesh type':'custom',
           'mesh cells x':round(PERIOD/(dx_nm*1e-9)),'mesh refinement':'staircase',
           'fit materials with multi-coefficient model':False,'incident angle':'table',
           'incident angle table':[[0,0]],'sample spacing':'custom','custom frequency samples':[299792458/WAVELENGTH],
           'report grating characterization':True,'report grating power':True,'report grating orders':True,
           'report field amplitudes':True,'report index':True})
    lines+=['save("model.fsp");','run;','save("model.fsp");']
    for name in DATASETS:lines.append(f'{name}=getresult("RCWA","{name}");')
    lines.append('matlabsave("native_datasets.mat",'+','.join(DATASETS)+');')
    return '\n'.join(lines)+'\n'

def extract(data,out):
    fa=data['field_amplitudes'];gc=data['grating_characterization'];gp=data['grating_power']
    # Metadata explicitly fixes axis order: lambda/f, theta/phi, n, m.
    groups=fa['Lumerical_dataset']['parameters']
    assert groups==[['lambda','f'],['theta','phi'],['n'],['m']],groups
    orders=np.asarray(fa['n']).ravel().astype(int)
    assert np.asarray(fa['m']).size==1
    k0=2*np.pi/WAVELENGTH;kx=2*np.pi/PERIOD*orders
    kz=np.sqrt((k0*k0-kx*kx).astype(complex))
    phase=np.exp(1j*(k0+kz)*REFERENCE)
    ey_raw=np.asarray(fa['eky_bs'])[0,0,:,0]
    # v252 normal-incidence native s input has Ey=-1. Convert to incident
    # canonical Ey=+1. A flat interface verifies this sign to <1e-13.
    ey=-ey_raw*phase
    # TE grating characterization is flux normalized; physical Ey uses sqrt(cos theta).
    rss=np.asarray(gc['Rss'])[0,0,:,0]
    power=(np.real(kz)/k0)*abs(ey)**2
    columns=['order','kx_per_m','kz_real_per_m','kz_imag_per_m','propagating','angle_deg','r_real','r_imag','r_abs','r_phase_deg','R_power','native_eky_bs_real','native_eky_bs_imag','native_Rss_real','native_Rss_imag']
    with (out/'reflection_orders.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader()
        for i,m in enumerate(orders):
            w.writerow(dict(zip(columns,[m,kx[i],kz[i].real,kz[i].imag,int(abs(kx[i])<k0),
                np.degrees(np.arcsin(kx[i]/k0)) if abs(kx[i])<k0 else '',ey[i].real,ey[i].imag,
                abs(ey[i]),np.degrees(np.angle(ey[i])),power[i],ey_raw[i].real,ey_raw[i].imag,rss[i].real,rss[i].imag])))
    keep=abs(orders)<=3
    np.savez_compressed(out/'reflection_propagating.npz',orders=orders[keep],r=ey[keep],R=power[keep],kx=kx[keep],kz=kz[keep])
    result={'reflection_total':float(power.sum()),'phase_reference_z_nm':0,
            'raw_reflection_reference_z_nm':REFERENCE*1e9,
            'raw_ey_key':'field_amplitudes.eky_bs','phase_deembedding':'r_m=-eky_bs_m*exp(i*(k0+kzm)*z_reference)',
            'TE_sign':'Ey incident = +1; validated against flat Fresnel at z=0',
            'orders':orders[keep].tolist(),'r':[complex(v) for v in ey[keep]],'R':power[keep].tolist()}
    return result

def one_run(s,nk,harmonics,dx_nm,flat=False):
    out=ROOT/('flat_calibration' if flat else f'harmonics_{harmonics:03d}_dx{dx_nm:g}nm')
    out.mkdir(parents=True,exist_ok=True)
    required=['summary.json','model.fsp','native_datasets.mat',
              'reflection_propagating.npz','field_amplitudes.npz']
    if all((out/name).exists() for name in required):
        print('REUSE '+out.name,flush=True)
        return json.loads((out/'summary.json').read_text())
    text=lsf(nk,harmonics,dx_nm,flat)
    (out/'build_run_export.lsf').write_text(text,encoding='utf-8')
    write_json(out/'input.json',{'wavelength_nm':532,'period_nm':2000,'duty_cycle':0.5,'depth_nm':0 if flat else 100,
        'Si_n':nk,'harmonic_max_order':harmonics,'expected_total_harmonics':2*harmonics+1,'dx_nm':dx_nm,
        'TE':'Ey','incident_direction':'+z','convention':'exp(-i omega t)','source_reference_z_nm':-200})
    s.cd(str(out));start=time.perf_counter();s.eval(text)
    # v252 lumapi cannot reshape this RCWA index rectilinear dataset when
    # dummy interfaces are present. Preserve it losslessly in native .mat;
    # the matrix datasets below translate correctly through lumapi.
    data={name:s.getresult('RCWA',name) for name in DATASETS if name!='index'}
    for name,dataset in data.items():archive_dataset(dataset,out,name)
    result=extract(data,out)
    result.update(elapsed_wall_seconds=time.perf_counter()-start,solver_version=s.version(),harmonics=harmonics,dx_nm=dx_nm)
    for name in ['num_k','simulation_run_time','total_threads','peak_memory']:
        try:result[name]=s.getresult('RCWA',name)
        except Exception as e:result[name+'_unavailable']=str(e)
    if flat:
        measured=result['r'][3]
        analytic=(1-nk)/(1+nk)
        result['fresnel_analytic']=analytic
        result['fresnel_complex_error']=abs(measured-analytic)
        assert result['fresnel_complex_error']<1e-8,(measured,analytic)
    write_json(out/'summary.json',result)
    print(out.name,json.dumps(result,default=json_default),flush=True)
    return result

def main():
    material=json.loads((ROOT/'material.json').read_text())
    nk=complex(material['getfdtdindex']['real'],material['getfdtdindex']['imag'])
    with session() as s:
        one_run(s,nk,20,5,True)
        for harmonics,dx in [(20,5),(40,5),(80,5),(80,2.5),(80,1.25),(160,1.25)]:one_run(s,nk,harmonics,dx)
    entries=[]
    for a,b in [('harmonics_020_dx5nm','harmonics_040_dx5nm'),('harmonics_040_dx5nm','harmonics_080_dx5nm'),('harmonics_080_dx5nm','harmonics_080_dx2.5nm')]:
        with np.load(ROOT/a/'reflection_propagating.npz') as p,np.load(ROOT/b/'reflection_propagating.npz') as q:
            entries.append({'from':a,'to':b,'max_complex_r_change':float(max(abs(p['r']-q['r']))),'max_power_R_change':float(max(abs(p['R']-q['R'])))})
    write_json(ROOT/'convergence.json',entries)
    print('CONVERGENCE',json.dumps(entries),flush=True)

if __name__=='__main__':main()
