"""Aggregate genuine solver diagnostics, with independent propagation checks."""
from pathlib import Path
import json,csv
import numpy as np
from run_fdtd import ROOT,parse_log,jswrite,C0
from quick_extract import extract,air_directory
Z0=376.730313668

def summarize():
    rows=[]
    for case in sorted((ROOT/'cases').glob('*')):
        ad=air_directory(case)
        if not (ad/'DONE.json').exists() or not (case/'mark/DONE.json').exists():continue
        cfg=json.loads((case/'mark/input.json').read_text(encoding='utf-8'))['settings']
        observed_file=case/'mark/observed_settings.json'
        observed=json.loads(observed_file.read_text(encoding='utf-8')) if observed_file.exists() else {}
        row=extract(case);row.update(N=cfg['N'],dx_nm=cfg['dx_m']/1e-9,dy_nm=cfg['dy_m']/1e-9,
          pml_layers=observed.get('pml layers',cfg['pml_layers']),requested_pml_layers=cfg['pml_layers'],max_time_fs=cfg['t_max_s']/1e-15)
        row['air_reference_reused']=(ad!=case/'air');row['air_reference']=str(ad.relative_to(ROOT/'cases'))
        for air in ['air','mark']:
            directory=ad if air=='air' else case/'mark'
            d=json.loads((directory/'DONE.json').read_text(encoding='utf-8'));stats=parse_log(directory/'engine_stdout.log')
            d['solver_statistics'].update(stats);jswrite(directory/'DONE.json',d)
            for k in ['simulated_time_s','final_auto_shutoff','terminated_early','fdtd_loop_wall_s','peak_memory_GiB','Nx','Ny']:row[air+'_'+k]=stats.get(k)
            row[air+'_solve_s']=d['solve_s']
        row['new_air_solve_s']=0. if row['air_reference_reused'] else row['air_solve_s']
        i=np.load(ad/'fields_I.npz');incfield=i['Ez'].ravel()
        a=np.load(ad/'fields_R.npz');b=np.load(case/'mark/fields_R.npz');a2=np.load(ad/'fields_R2.npz');b2=np.load(case/'mark/fields_R2.npz')
        x=b['x'].ravel();L=cfg['N']*cfg['pitch_m'];q=(x>=-L/2-1e-12)&(x<L/2-1e-12)
        k0=2*np.pi/cfg['wavelength_m'];orders=np.arange(-3,4);kx=orders*2*np.pi/cfg['pitch_m'];ky=np.sqrt(k0*k0-kx*kx)
        phase=np.exp(-1j*kx[:,None]*x[q]);inc=incfield[q].mean()
        e=phase@((b['Ez']-a['Ez']).ravel()[q])/q.sum()/inc
        h=phase@((b['Hx']-a['Hx']).ravel()[q])/q.sum()/inc
        e2=phase@((b2['Ez']-a2['Ez']).ravel()[q])/q.sum()/inc
        r=e*np.exp(-1j*ky*cfg['R_y_m']);r2=e2*np.exp(-1j*ky*cfg['R2_y_m'])
        row['R_to_R2_deprop_max_complex_difference']=float(np.max(abs(r2-r)))
        ep=(e+Z0*k0/ky*h)/2;em=(e-Z0*k0/ky*h)/2
        row['R_max_inward_over_outward_amplitude']=float(np.max(abs(em)/abs(ep)))
        rows.append(row)
    jswrite(ROOT/'solver_manifest.json',rows)
    flat=[{k:v for k,v in r.items() if not isinstance(v,list)} for r in rows]
    if flat:
        with (ROOT/'solver_manifest.csv').open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    indexed={r['case']:r for r in rows};controls=[]
    for kind,baseline,changed in [
      ('maximum simulation time, with automatic early stop allowed','N200_dx20_dy10_t600','N200_dx20_dy10_t1200'),
      ('forced longer integration, stricter autoshutoff; existing full-duration air reused','N200_dx20_dy10_t600','N200_dx20_dy10_t1200_strictmark'),
      ('mesh at fixed physical PML thickness','N050_dx20_dy10_t600','N050_dx10_dy5_t600'),
      ('PML thickness only','N020_dx10_dy5_t600','N020_dx10_dy5_t600_pml48'),
      ('mesh at fixed physical PML thickness','N020_dx20_dy10_t600','N020_dx10_dy5_t600_pml48'),
      ('mesh refinement; finest effective PML capped at 64, thickness not fully matched','N020_dx10_dy5_t600_pml48','N020_dx5_dy2.5_t600')]:
        if baseline not in indexed or changed not in indexed:continue
        aa=indexed[baseline];bb=indexed[changed];ra=np.array(aa['r_re'])+1j*np.array(aa['r_im']);rb=np.array(bb['r_re'])+1j*np.array(bb['r_im']);delta=rb-ra
        controls.append(dict(control=kind,baseline=baseline,changed=changed,orders=list(range(-3,4)),
          baseline_mark_actual_time_s=aa['mark_simulated_time_s'],changed_mark_actual_time_s=bb['mark_simulated_time_s'],
          complex_delta_real=delta.real.tolist(),complex_delta_imag=delta.imag.tolist(),absolute_complex_delta=abs(delta).tolist(),
          max_absolute_complex_delta=float(max(abs(delta))),relative_complex_L2_delta=float(np.linalg.norm(delta)/np.linalg.norm(ra))))
    jswrite(ROOT/'numerical_controls.json',controls)
    for r in rows:print(r['case'],'max propagation change',r['R_to_R2_deprop_max_complex_difference'],'incoming/outgoing',r['R_max_inward_over_outward_amplitude'])
    return rows

if __name__=='__main__':summarize()
