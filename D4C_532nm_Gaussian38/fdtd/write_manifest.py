"""Record actual completed engine solves, raw-field integrity and references."""
from pathlib import Path
import csv,hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parent

def main():
    rows=[]
    for path in sorted((ROOT/'cases').glob('*/DONE.json')):
        folder=path.parent;d=json.loads(path.read_text(encoding='utf-8'));c=d['settings'];term=d['termination'];obs=json.loads((folder/'observed_settings.json').read_text(encoding='utf-8'));v=c.get('geometry',{})
        record=folder/'AIR_REFERENCE.json'
        ref=json.loads(record.read_text(encoding='utf-8')) if record.exists() else {}
        source=ref.get('air_reference','air'+folder.name[len(c['case']):])
        r=dict(case=folder.name,physical_kind=c['case'],geometry=v.get('geometry','reference'),cap_nm=v.get('cap_nm'),wavelength_nm=532.,dx_nm=c['dx_m']*1e9,dy_nm=c['dy_m']*1e9,pml_requested=c['pml_layers'],pml_actual_max=int(np.max(np.asarray(obs['pml layers']))),actual_all_PML=all(obs[p]=='PML' for p in ('x min bc','x max bc','y min bc','y max bc')),requested_time_fs=c['t_max_s']*1e15,actual_time_fs=term.get('simulated_time_s',0)*1e15,auto_shutoff_requested=c['auto_shutoff'],auto_shutoff_final=term.get('final_auto_shutoff'),met_requested_shutoff=term.get('final_auto_shutoff',float('inf'))<=c['auto_shutoff'],iterations=term.get('iterations'),engine_returncode=json.loads((folder/'engine_result.json').read_text(encoding='utf-8'))['returncode'],engine_success=term['completed_successfully'],solve_wall_s=d['solve_wall_s'],peak_memory_GiB=term.get('peak_memory_GiB'),dt_s=d['exports']['dt_s'],sourcepower_W_per_m=d['exports']['sourcepower_W_per_m'],reference_air_case=source if c['case']!='air' else folder.name,air_reference_reused=bool(ref),model_sha256=d['model_sha256'])
        valid=True;hashes={}
        for mon in ('R','R2','INC','T'):
            p=folder/f'nearfield_{mon}.npz';hashes[mon]=hashlib.sha256(p.read_bytes()).hexdigest()
            with np.load(p) as z:
                valid &= abs(float(z['wavelength_m'])-532e-9)<1e-15
                valid &= all(k in z and np.iscomplexobj(z[k]) and np.all(np.isfinite(z[k])) for k in ('Ex','Ey','Ez','Hx','Hy','Hz'))
        r['all_raw_complex_fields_valid']=bool(valid);rows.append(r)
        r['ridge_height_nm']=v.get('ridge_height_nm')
        r['ridge_width_um']=v.get('ridge_width_um')
        r['periods']=v.get('periods')
        (folder/'raw_field_hashes.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
    cancelled=[dict(case=p.parent.name,**json.loads(p.read_text(encoding='utf-8'))) for p in sorted((ROOT/'cases').glob('*/CANCELLED.json'))]
    out=dict(status='passed' if all(r['engine_success'] and r['engine_returncode']==0 and r['actual_all_PML'] and r['all_raw_complex_fields_valid'] for r in rows) else 'failed',actual_completed_case_count=len(rows),cancelled_case_count=len(cancelled),cancelled_cases=cancelled,one_engine_at_a_time=True,solver='genuine installed Lumerical FDTD v252, one CPU process, 8 threads',cases=rows)
    (ROOT/'solver_manifest.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    with (ROOT/'solver_manifest.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else []);w.writeheader();w.writerows(rows)
    print(json.dumps(dict(status=out['status'],actual_completed_case_count=len(rows))))

if __name__=='__main__':main()
