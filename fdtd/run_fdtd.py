"""Finite, PML-bounded Si mark: genuine Lumerical FDTD with archived fields.

Native 2D XY coordinates; Ez TE; incident plane wave travels toward -y.
All dimensions in SI. Run with the installed Lumerical Python interpreter.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, json, time, traceback, shutil, re, os
import numpy as np
from runtime_support import hidden_fdtd, run_cpu_engine

ROOT=Path(__file__).resolve().parent
ASCII_ROOT=Path(os.environ.get('SI532_NATIVE_ROOT',r'D:\Si532_spectrum_native_20261031'))/'fdtd'
C0=299792458.0

def jswrite(path,data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else str(x)),encoding='utf-8')

def params(n,dx_nm=20,dy_nm=10,t_fs=600,shutoff=1e-7,pml=24):
    return dict(N=int(n),pitch_m=2e-6,duty=0.5,etch_depth_m=100e-9,wavelength_m=532e-9,
      buffer_m=4e-6,dx_m=dx_nm*1e-9,dy_m=dy_nm*1e-9,t_max_s=t_fs*1e-15,
      auto_shutoff=shutoff,pml_layers=pml,y_min_m=-0.7e-6,y_max_m=1.8e-6,
      source_y_m=0.7e-6,R_y_m=1e-6,R2_y_m=1.2e-6,I_y_m=0.,
      material='Si (Silicon) - Palik',source_wavelength_start_m=500e-9,
      source_wavelength_stop_m=565e-9,polarization='Ez (TE)',boundary='all PML',
      interpolation='specified position',mesh_refinement='conformal variant 0')

def quote(v):
    if isinstance(v,str):return json.dumps(v)
    return format(float(np.asarray(v).ravel()[0]),'.17g')

def parse_log(path):
    path=Path(path)
    if not path.exists():path=path.with_name(path.stem+'_public.log')
    text=path.read_text(encoding='utf-8',errors='replace');out={}
    for key,pat in [('iterations',r'Completed (\d+) iterations'),
      ('simulated_time_s',r'or ([\deE.+-]+)s of Simulation Time'),
      ('peak_memory_GiB',r'Peak CPU memory used in the simulation \(GiB\): ([\deE.+-]+)'),
      ('engine_overall_wall_s',r'Overall wall time measurements in seconds: ([\deE.+-]+)'),
      ('fdtd_loop_wall_s',r'time to run FDTD simulation: ([\deE.+-]+)'),
      ('Nx',r'Total X spatial steps (\d+)'),('Ny',r'Total Y spatial steps (\d+)')]:
        hits=re.findall(pat,text)
        if hits:out[key]=float(hits[-1]) if key not in ['iterations','Nx','Ny'] else int(hits[-1])
    hits=re.findall(r'Auto Shutoff: ([\deE.+-]+)',text)
    if hits:out['final_auto_shutoff']=float(hits[-1])
    out['terminated_early']='Early termination of simulation' in text
    out['completed_successfully']='Simulation completed successfully' in text
    return out

def build(s,cfg,air):
    s.switchtolayout();s.deleteall();lines=['# SI units. Finite mark, all boundaries PML.','newproject;']
    def add(kind,properties):
        getattr(s,kind)();lines.append(kind+';')
        for k,v in properties.items():s.set(k,v);lines.append('set('+quote(k)+','+quote(v)+');')
    width=cfg['N']*cfg['pitch_m']+2*cfg['buffer_m']
    add('addfdtd',{'dimension':1,'x':0.,'x span':width,'y min':cfg['y_min_m'],
      'y max':cfg['y_max_m'],'z':0.,'mesh type':'uniform','dx':cfg['dx_m'],'dy':cfg['dy_m'],
      'mesh refinement':cfg['mesh_refinement'],'simulation time':cfg['t_max_s'],
      'auto shutoff min':cfg['auto_shutoff'],'pml profile':1,'pml layers':cfg['pml_layers'],
      'x min bc':'PML','x max bc':'PML','y min bc':'PML','y max bc':'PML'})
    if not air:
        add('addrect',{'name':'silicon_substrate','material':cfg['material'],'x':0.,'x span':width+4e-6,
            'y min':cfg['y_min_m']-1e-6,'y max':-cfg['etch_depth_m'],'z':0.,'z span':1e-6})
        for j in range(cfg['N']):
            center=(j-(cfg['N']-1)/2)*cfg['pitch_m']
            add('addrect',{'name':'ridge_%04d'%j,'material':cfg['material'],'x':center,
              'x span':cfg['duty']*cfg['pitch_m'],'y min':-cfg['etch_depth_m'],'y max':0.,'z':0.,'z span':1e-6})
    add('addplane',{'name':'source','injection axis':'y','direction':'Backward',
      'plane wave type':'Bloch/periodic','x':0.,'x span':width+2*cfg.get('source_extension_m',0.),'y':cfg['source_y_m'],'z':0.,
      'polarization angle':90.,'angle theta':0.,'wavelength start':cfg['source_wavelength_start_m'],
      'wavelength stop':cfg['source_wavelength_stop_m']})
    for name in ['R','R2']+(['I'] if air else []):
        add('adddftmonitor',{'name':name,'monitor type':'Linear X','x':0.,'x span':width,
          'y':cfg[name+'_y_m'],'z':0.,'override global monitor settings':1,
          'sample spacing':'custom','custom frequency samples':np.array([C0/cfg['wavelength_m']]),
          'spatial interpolation':cfg['interpolation']})
    return '\n'.join(lines)+'\n'

def export_monitor(s,name,path):
    fields={}
    for f in ['x','y','z','f','Ex','Ey','Ez','Hx','Hy','Hz','delta_x','delta_y','delta_z']:
        try:fields[f]=np.asarray(s.getdata(name,f,1))
        except Exception:
            if not f.startswith('delta_'):raise
    fields['x_m']=fields['x'];fields['y_m']=fields['y'];fields['frequency_hz']=fields['f']
    fields['transmission']=np.asarray(s.transmission(name))
    np.savez_compressed(path,**fields)
    return {k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in fields.items()}

def run_one(s,cfg,air,tag=None):
    if not ASCII_ROOT.exists() or ASCII_ROOT.resolve()!=ROOT.resolve():
        raise RuntimeError('SI532_NATIVE_ROOT must point to this project root through an ASCII path or junction.')
    base=tag or 'N%03d_dx%g_dy%g_t%g'%(cfg['N'],cfg['dx_m']/1e-9,cfg['dy_m']/1e-9,cfg['t_max_s']/1e-15)
    directory=ROOT/'cases'/base/('air' if air else 'mark');directory.mkdir(parents=True,exist_ok=True)
    asc=ASCII_ROOT/directory.relative_to(ROOT)
    done=directory/'DONE.json'
    if done.exists():
        d=json.loads(done.read_text(encoding='utf-8'))
        if d.get('status')=='completed' and d['settings']==cfg:return d
    domain_width=cfg['N']*cfg['pitch_m']+2*cfg['buffer_m']
    jswrite(directory/'input.json',{'settings':cfg,'air':air,'domain_width_m':domain_width,
      'source_width_m':domain_width+2*cfg.get('source_extension_m',0.)})
    print('BUILD',base,'air' if air else 'mark',flush=True)
    script=build(s,cfg,air);(directory/'model.lsf').write_text(script+'save("model.fsp");\nrun;\nsave("model.fsp");\n',encoding='utf-8')
    observed={}
    for name in ['x span','y min','y max','dx','dy','pml layers','pml profile','x min bc','x max bc','y min bc','y max bc','simulation time','auto shutoff min','dt stability factor']:
        try:observed[name]=s.getnamed('FDTD',name)
        except Exception as ex:observed[name]={'read_error':str(ex)}
    observed['source x span']=s.getnamed('source','x span')
    jswrite(directory/'observed_settings.json',observed)
    material={}
    for method in ['getindex','getfdtdindex']:
        try:
            if method=='getfdtdindex':v=s.getfdtdindex(cfg['material'],C0/cfg['wavelength_m'],C0/cfg['source_wavelength_stop_m'],C0/cfg['source_wavelength_start_m'])
            else:v=s.getindex(cfg['material'],C0/cfg['wavelength_m'])
            v=complex(np.asarray(v).ravel()[0]);material[method]={'n_real':v.real,'n_imag':v.imag,'epsilon_real':(v*v).real,'epsilon_imag':(v*v).imag}
        except Exception as ex:material[method]={'error':str(ex)}
    jswrite(directory/'material.json',material)
    s.save(str(asc/'model.fsp'))
    print('RUN',base,'air' if air else 'mark',flush=True)
    result=run_cpu_engine(asc/'model.fsp',asc,threads=8,timeout=3600)
    if result['returncode']!=0:raise RuntimeError(str(result))
    s.load(str(asc/'model.fsp'))
    monitors={name:export_monitor(s,name,directory/('fields_'+name+'.npz')) for name in ['R','R2']+(['I'] if air else [])}
    stats=parse_log(directory/'engine_stdout.log')
    for mon in ['FDTD']:
        try:stats['FDTD_results']=str(s.getresult(mon))
        except Exception as ex:stats['FDTD_results_error']=str(ex)
        for attr in ['status','simulation time','auto shutoff min']:
            try:stats[attr]=np.asarray(s.getnamed(mon,attr)).tolist()
            except Exception:pass
    # The genuine engine stdout records cycle counts and autoshutoff progress.
    out=dict(status='completed',settings=cfg,air=air,solve_s=result['elapsed_seconds'],
      engine=result,monitors=monitors,material=material,solver_statistics=stats,
      completed_utc=datetime.now(timezone.utc).isoformat(),
      phase_convention='Lumerical phasors exp(-i omega t), outgoing +y exp(+i ky y)')
    jswrite(done,out);print('DONE',base,'air' if air else 'mark','%.2f s'%out['solve_s'],flush=True)
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--N',nargs='+',type=int,default=[10]);p.add_argument('--dx',type=float,default=20);p.add_argument('--dy',type=float,default=10);p.add_argument('--time-fs',type=float,default=600);p.add_argument('--shutoff',type=float,default=1e-7);p.add_argument('--pml',type=int,default=24);p.add_argument('--tag');p.add_argument('--mark-only',action='store_true');p.add_argument('--air-only',action='store_true');p.add_argument('--source-extension-um',type=float,default=0.);args=p.parse_args()
    with hidden_fdtd(threads=8) as s:
        for n in args.N:
            cfg=params(n,args.dx,args.dy,args.time_fs,args.shutoff,args.pml)
            if args.source_extension_um:cfg['source_extension_m']=args.source_extension_um*1e-6
            for air in ([True] if args.air_only else [False] if args.mark_only else [True,False]):run_one(s,cfg,air,args.tag)

if __name__=='__main__':main()
