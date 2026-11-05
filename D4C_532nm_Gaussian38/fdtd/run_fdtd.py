"""Genuine isolated all-PML Gaussian finite-mark solves for D4C zero-order audit.

The script writes only in this project's fdtd/ directory.  NPZ contains raw
six-component complex E/H; no flat-background subtraction is applied here.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, csv, hashlib, json, os, re, sys, time, traceback
import numpy as np
from runtime_support import hidden_fdtd, run_cpu_engine

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent
NATIVE=Path(os.environ.get('AH53_D4C_NATIVE_ROOT',r'D:\AH53_D4C_zero_native_20261105'))/'fdtd'
C0=299792458.
SI='Si (Silicon) - Palik'
OX='SiO2 (Glass) - Palik'
COMPS=('Ex','Ey','Ez','Hx','Hy','Hz')

def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,complex):return dict(re=x.real,im=x.imag)
    if isinstance(x,Path):return str(x)
    return x

def dump(path,data):
    Path(path).write_text(json.dumps(clean(data),ensure_ascii=False,indent=2),encoding='utf-8')

def quote(x):
    if isinstance(x,str):return json.dumps(x)
    if isinstance(x,np.ndarray):return '['+','.join(format(float(v),'.17g') for v in x.ravel())+']'
    return format(float(x),'.17g')

def probe_materials(s):
    freq=C0/532e-9
    result=dict(wavelength_nm=532.,fit_band_nm=[500.,600.],solver_version=s.version(),materials={})
    for key,name in [('Si',SI),('SiO2',OX)]:
        db=complex(np.asarray(s.getindex(name,freq)).squeeze())
        fit=complex(np.asarray(s.getfdtdindex(name,freq,C0/600e-9,C0/500e-9)).squeeze())
        result['materials'][key]=dict(name=name,database_index=dict(n=db.real,k=db.imag),FDTD_fit_index=dict(n=fit.real,k=fit.imag,epsilon_re=(fit*fit).real,epsilon_im=(fit*fit).imag))
    result['source']=dict(amplitude=1.,phase_deg=0.,waist_radius_um=19.,intensity_1e2_diameter_um=38.,waist_y_um=0.,source_y_um=1.35)
    result['phase_convention']='exp(-i omega t); incident -Y, reflected +Y; common reference y=0.85um; raw fields retain solver source phase'
    dump(ROOT/'materials.json',result)
    print('MATERIAL_PROBE',json.dumps(result),flush=True)
    return result

def parse_log(path):
    if not Path(path).exists():return {}
    txt=Path(path).read_text(encoding='utf-8',errors='replace');r={}
    for key,pat in [('iterations',r'Completed (\d+) iterations'),('simulated_time_s',r'or ([\deE.+-]+)s of Simulation Time'),('peak_memory_GiB',r'Peak CPU memory used in the simulation \(GiB\): ([\deE.+-]+)'),('wall_s',r'Overall wall time measurements in seconds: ([\deE.+-]+)'),('Nx',r'Total X spatial steps (\d+)'),('Ny',r'Total Y spatial steps (\d+)')]:
        values=re.findall(pat,txt)
        if values:r[key]=int(values[-1]) if key in ('iterations','Nx','Ny') else float(values[-1])
    vals=re.findall(r'Auto Shutoff:\s*([\deE.+-]+)',txt)
    if vals:r['final_auto_shutoff']=float(vals[-1])
    r['completed_successfully']='Simulation completed successfully' in txt
    r['terminated_early']='Early termination of simulation' in txt
    return r

def cfg(case,dx=20.,dy=10.,pml=32,t=600.,shut=1e-9):
    variants={v['id']:v for v in json.loads((PROJECT/'experiment.json').read_text(encoding='utf-8'))['variants']}
    extra_path=PROJECT/'inputs/extra_weak_mark.json'
    extra_ids=set()
    if extra_path.exists():
        extra=json.loads(extra_path.read_text(encoding='utf-8'))
        extra_variants=extra.get('variants',[extra.get('variant',extra)])
        for v in extra_variants:
            if not isinstance(v,dict) or 'id' not in v:continue
            if v['id'] in variants:raise RuntimeError('Extra holdout must not override the original8: '+v['id'])
            variants[v['id']]=v;extra_ids.add(v['id'])
    c=dict(case=case,wavelength_m=532e-9,dx_m=dx*1e-9,dy_m=dy*1e-9,pml_layers=pml,t_max_s=t*1e-15,auto_shutoff=shut,x_span_m=192e-6,monitor_span_m=184e-6,y_min_m=-.8e-6,y_max_m=2.5e-6,source_y_m=1.35e-6,source_span_m=200e-6,waist_radius_m=19e-6,reference_y_m=.85e-6,source_band_m=[500e-9,600e-9],monitor_y_m={'R':1.8e-6,'R2':2e-6,'INC':.8e-6,'T':-.6e-6},mesh_refinement='conformal variant 0')
    if case not in ('air','flatSi'):
        c['geometry']=variants[case]
        if case in extra_ids:c['variant_source']=dict(path='inputs/extra_weak_mark.json',geometry_sha256=hashlib.sha256(json.dumps(variants[case],sort_keys=True,separators=(',',':')).encode()).hexdigest(),held_out_from_original_eight=True)
    return c

def build(s,c):
    s.switchtolayout();s.deleteall();script=['# SI units; new project; finite Gaussian, all four boundaries PML.','newproject;']
    def add(kind,properties):
        getattr(s,kind)();script.append(kind+';')
        for k,v in properties.items():s.set(k,v);script.append('set('+quote(k)+','+quote(v)+');')
    add('addfdtd',{'dimension':'2D','x':0.,'x span':c['x_span_m'],'y min':c['y_min_m'],'y max':c['y_max_m'],'z':0.,'mesh type':'uniform','dx':c['dx_m'],'dy':c['dy_m'],'mesh refinement':c['mesh_refinement'],'simulation time':c['t_max_s'],'auto shutoff min':c['auto_shutoff'],'pml profile':1,'pml layers':c['pml_layers'],'x min bc':'PML','x max bc':'PML','y min bc':'PML','y max bc':'PML'})
    if c['pml_layers']>64:
        probe=json.loads((ROOT/'pml_128_probe.json').read_text(encoding='utf-8'))['success']
        if not probe:raise RuntimeError('Actual PML128 capability not established')
        for prop,value in [('pml profile',probe['profile']),(probe['max_layer_property'],c['pml_layers']),('pml layers',c['pml_layers'])]:
            s.set(prop,value);script.append('set('+quote(prop)+','+quote(value)+');')
        actual=int(np.max(np.asarray(s.getnamed('FDTD','pml layers'))))
        if actual!=c['pml_layers']:raise RuntimeError('PML layer request clamped: '+str(actual))
    # Rectangles extend through external PML; lower mesh order takes priority.
    def rect(name,x,xspan,ymin,ymax,mat,order=2):
        add('addrect',{'name':name,'x':x,'x span':xspan,'y min':ymin,'y max':ymax,'z':0.,'z span':1e-6,'material':mat,'override mesh order from material database':1,'mesh order':order})
    if c['case']!='air':
        rect('Si_substrate',0.,c['x_span_m']+4e-6,c['y_min_m']-2e-6,0.,SI,1)
        if 'geometry' in c:
            v=c['geometry'];cap=v['cap_nm']*1e-9;h=v['ridge_height_nm']*1e-9
            # cap000 is bare Si; positive cap fills all grooves and covers ridges.
            # The planarization film extends laterally into the PML.
            if cap>0:rect('SiO2_planarization',0.,c['x_span_m']+4e-6,0.,h+cap,OX,2)
            period=v['period_um']*1e-6
            for cell in range(v['periods']):
                cx=(cell-(v['periods']-1)/2)*period
                for j,local in enumerate(v['ridge_centers_local_um']):rect(f'Si_ridge_{cell:02d}_{j}',cx+local*1e-6,v['ridge_width_um']*1e-6,0.,h,SI,1)
    add('addgaussian',{'name':'source','injection axis':'y-axis','direction':'backward','use scalar approximation':1,'waist radius w0':c['waist_radius_m'],'distance from waist':-c['source_y_m'],'x':0.,'x span':c['source_span_m'],'y':c['source_y_m'],'z':0.,'angle theta':0.,'polarization angle':90.,'amplitude':1.,'phase':0.,'wavelength start':c['source_band_m'][0],'wavelength stop':c['source_band_m'][1]})
    for name,y in c['monitor_y_m'].items():add('addpower',{'name':name,'monitor type':'Linear X','x':0.,'x span':c['monitor_span_m'],'y':y,'z':0.,'spatial interpolation':'specified position','override global monitor settings':1,'sample spacing':'custom','custom frequency samples':np.array([C0/c['wavelength_m']])})
    return '\n'.join(script)+'\n'

def export(s,c,folder):
    meta={}
    for name in c['monitor_y_m']:
        freq=float(np.asarray(s.getdata(name,'f')).squeeze());x=np.asarray(s.getdata(name,'x')).real.squeeze();y=float(np.asarray(s.getdata(name,'y')).squeeze())
        d=dict(x_m=x,monitor_y_m=y,frequency_hz=freq,wavelength_m=C0/freq,**{f:np.asarray(s.getdata(name,f)).squeeze().astype(complex) for f in COMPS})
        for axis in ('x','y','z'):
            try:d['delta_'+axis]=np.asarray(s.getdata(name,'delta_'+axis)).squeeze()
            except Exception:pass
        np.savez_compressed(folder/f'nearfield_{name}.npz',**d)
        meta[name]=dict(points=len(x),x_min_m=x.min(),x_max_m=x.max(),monitor_y_m=y,wavelength_m=C0/freq,native_transmission_signed=float(np.asarray(s.transmission(name)).squeeze()))
    mesh={axis+'_m':np.asarray(s.getresult('FDTD',axis)).real.squeeze() for axis in ('x','y')}
    np.savez_compressed(folder/'mesh_coordinates.npz',**mesh)
    props=['dimension','x min bc','x max bc','y min bc','y max bc','pml profile','pml layers','mesh type','dx','dy','dt','simulation time','auto shutoff min','x span','y min','y max']
    observed={p:s.getnamed('FDTD',p) for p in props}
    observed['mesh_actual']={axis:dict(points=len(v),min_m=float(v.min()),max_m=float(v.max()),step_min_m=float(np.diff(v).min()),step_max_m=float(np.diff(v).max())) for axis,v in mesh.items()}
    observed['source']={p:s.getnamed('source',p) for p in ('x','x span','y','amplitude','phase','polarization angle','waist radius w0','distance from waist','wavelength start','wavelength stop')}
    observed['geometry']=dict(kind=c['case'],substrate_present=c['case']!='air',ridge_count=(c.get('geometry',{}).get('periods',0)*len(c.get('geometry',{}).get('ridge_centers_local_um',[]))),oxide_present=c.get('geometry',{}).get('cap_nm',0)>0)
    for name in ['Si_substrate']+(['SiO2_planarization'] if observed['geometry']['oxide_present'] else []):
        if c['case']!='air':observed['geometry'][name]={p:s.getnamed(name,p) for p in ('material','x','x span','y min','y max','mesh order')}
    if 'geometry' in c:
        observed['geometry']['Si_ridge_00_0']={p:s.getnamed('Si_ridge_00_0',p) for p in ('material','x','x span','y min','y max','mesh order')}
    dump(folder/'observed_settings.json',observed)
    assert all(str(observed[p])=='PML' for p in ('x min bc','x max bc','y min bc','y max bc'))
    for axis in ('x','y'):
        steps=np.diff(mesh[axis+'_m'])
        assert np.ptp(steps)<c['d'+axis+'_m']*1e-5, (axis,steps.min(),steps.max())
    return dict(monitors=meta,sourcepower_W_per_m=float(np.asarray(s.sourcepower(C0/c['wavelength_m'])).squeeze()),dt_s=float(np.asarray(observed['dt']).squeeze()),mesh=observed['mesh_actual'])

def sanitize_log(folder):
    p=folder/'engine_stdout.log'
    txt=p.read_text(encoding='utf-8',errors='replace')
    lines=[]
    for line in txt.splitlines():
        if any(w in line.lower() for w in ('license','running on host','running fsp file','@localhost','username','\\users\\')):continue
        line=line.replace(str(NATIVE),'[project]/fdtd').replace(str(ROOT),'[project]/fdtd')
        lines.append(line)
    (folder/'engine_stdout_public.log').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def run_one(s,c,tag):
    folder=ROOT/'cases'/tag;folder.mkdir(parents=True,exist_ok=True)
    native=NATIVE/'cases'/tag
    if (folder/'CANCELLED.json').exists():raise RuntimeError('Preserved cancelled case; use a new suffix rather than overwriting its diagnostic record: '+tag)
    if (folder/'DONE.json').exists():
        old=json.loads((folder/'DONE.json').read_text(encoding='utf-8'))
        if old['settings']!=clean(c):raise RuntimeError('Existing DONE settings mismatch: '+tag)
        print('SKIP_COMPLETED',tag,flush=True);return old
    dump(folder/'settings.json',c)
    project=native/'model.fsp';t0=time.perf_counter()
    prior=parse_log(folder/'engine_stdout.log')
    if project.exists() and prior.get('completed_successfully'):
        print('RESUME_EXPORT',tag,flush=True);s.load(str(project));wall=prior.get('wall_s')
    else:
        script=build(s,c);(folder/'build_this_case.lsf').write_text(script,encoding='utf-8')
        (folder/'run_this_case.lsf').write_text('case_dir=filedirectory(currentscriptname);\n'+script+'save(case_dir+"/reproduction.fsp");\nrun("FDTD");\nsave;\n',encoding='utf-8')
        s.save(str(project));print('START_SOLVE',tag,json.dumps(clean(c)),flush=True)
        t=time.perf_counter();result=run_cpu_engine(project,folder,threads=8);wall=time.perf_counter()-t
        dump(folder/'engine_result.json',result)
        if result.get('returncode',result.get('return_code',0))!=0:raise RuntimeError('Engine failure '+tag)
        s.load(str(project))
    term=parse_log(folder/'engine_stdout.log');sanitize_log(folder)
    if not term.get('completed_successfully'):raise RuntimeError('No successful engine record '+tag)
    exp=export(s,c,folder)
    out=dict(status='completed',settings=c,termination=term,solve_wall_s=wall,total_wall_s=time.perf_counter()-t0,exports=exp,completed_utc=datetime.now(timezone.utc).isoformat(),model_sha256=hashlib.sha256((folder/'model.fsp').read_bytes()).hexdigest(),raw_field_policy='TOTAL fields, no flat subtraction; all six raw complex components retained')
    dump(folder/'DONE.json',out);print('DONE',tag,json.dumps(clean(out)),flush=True);return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cases',nargs='+',default=['air','flatSi','ah53_cap000']);ap.add_argument('--dx-nm',type=float,default=20);ap.add_argument('--dy-nm',type=float,default=10);ap.add_argument('--pml',type=int,default=32);ap.add_argument('--time-fs',type=float,default=600);ap.add_argument('--shutoff',type=float,default=1e-9);ap.add_argument('--suffix',default='');ap.add_argument('--probe-only',action='store_true');ap.add_argument('--air-reference',help='Completed same-grid empty-domain case explicitly reused, e.g. air_fine');args=ap.parse_args()
    if not NATIVE.exists() or not NATIVE.samefile(ROOT):raise RuntimeError('ASCII root does not resolve to this project')
    with hidden_fdtd(8) as s:
        probe_materials(s)
        if args.probe_only:return
        for case in args.cases:
            c=cfg(case,args.dx_nm,args.dy_nm,args.pml,args.time_fs,args.shutoff);tag=case+args.suffix
            if args.air_reference and case!='air':
                ref=ROOT/'cases'/args.air_reference
                refdone=json.loads((ref/'DONE.json').read_text(encoding='utf-8'))
                check=('dx_m','dy_m','pml_layers','x_span_m','monitor_span_m','y_min_m','y_max_m','source_y_m','source_span_m','waist_radius_m','source_band_m','monitor_y_m')
                assert all(refdone['settings'][k]==clean(c[k]) for k in check),'Reference grid/source mismatch'
                assert refdone['termination']['simulated_time_s']<c['t_max_s'],'Reference did not stop before the requested mark duration'
                assert refdone['termination'].get('final_auto_shutoff',float('inf'))<=c['auto_shutoff'],'Reference did not meet the requested mark shutoff threshold'
                folder=ROOT/'cases'/tag;folder.mkdir(parents=True,exist_ok=True)
                dump(folder/'AIR_REFERENCE.json',dict(air_reference=args.air_reference,path_base='fdtd/cases',reference_reused=True,reason='Identical source/domain/grid/PML; empty-domain pulse already met stricter-or-equal requested shutoff before mark end time; no duplicate air solve',reference_actual_time_s=refdone['termination']['simulated_time_s'],reference_final_auto_shutoff=refdone['termination'].get('final_auto_shutoff')))
            run_one(s,c,tag)

if __name__=='__main__':
    try:main()
    except Exception:
        traceback.print_exc();sys.exit(1)
