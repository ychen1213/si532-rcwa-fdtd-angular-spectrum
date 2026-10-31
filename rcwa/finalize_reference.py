"""Produce an audited RCWA coefficient adapter with explicit phase origins.

Original solver datasets and reflection_orders.csv are never overwritten.
"""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
from run_rcwa import ROOT,PERIOD,WAVELENGTH,REFERENCE,write_json

RUNS=['harmonics_020_dx5nm','harmonics_040_dx5nm','harmonics_080_dx5nm',
      'harmonics_080_dx2.5nm','harmonics_080_dx1.25nm','harmonics_160_dx1.25nm']

def load(name):
    out=ROOT/name
    mesh=json.loads((out/'actual_mesh.json').read_text())
    summary=json.loads((out/'summary.json').read_text())
    with np.load(out/'reflection_propagating.npz') as a:
        m=a['orders'].copy(); r=a['r'].copy(); power=a['R'].copy()
    actual_cells=int(mesh['actual mesh cells x']);dx=PERIOD/actual_cells
    correction=np.exp(-1j*(2*np.pi/PERIOD*m)*dx/2)
    centered=r*correction
    with np.load(out/'sampled_index.npz') as a:
        x=a['x'];n=a['index_z'][:,0,2,0]
        sampled_ridge_center=float(np.mean(x[np.abs(n)>2]))
    audit=[]
    for order in [1,2,3]:
        pos=int(np.flatnonzero(m==order)[0]);neg=int(np.flatnonzero(m==-order)[0])
        raw_ratio=r[pos]/r[neg]
        expected=np.exp(1j*(2*np.pi/PERIOD*order)*dx)
        audit.append({'run':name,'m':order,'actual_cells_x':actual_cells,'actual_dx_nm':dx*1e9,
            'sampled_ridge_center_nm':sampled_ridge_center*1e9,
            'raw_pair_phase_difference_deg':np.degrees(np.angle(raw_ratio)),
            'expected_half_cell_pair_phase_deg':np.degrees(np.angle(expected)),
            'raw_pair_phase_law_error':abs(raw_ratio-expected),
            'corrected_pair_complex_error':abs(centered[pos]-centered[neg])})
    return {'m':m,'r':r,'centered':centered,'R':power,'dx':dx,'mesh':actual_cells,
            'correction':correction,'summary':summary,'audit':audit}

def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    runs={name:load(name) for name in RUNS}
    write_csv(ROOT/'phase_origin_audit.csv',[item for run in runs.values() for item in run['audit']])
    convergence=[]
    for a,b in zip(RUNS[:-1],RUNS[1:]):
        low=runs[a];high=runs[b]
        convergence.append({'from':a,'to':b,'from_actual_dx_nm':low['dx']*1e9,'to_actual_dx_nm':high['dx']*1e9,
            'max_complex_r_change_before_lateral_correction':float(max(abs(high['r']-low['r']))),
            'max_complex_r_change_centered':float(max(abs(high['centered']-low['centered']))),
            'max_R_change':float(max(abs(high['R']-low['R']))),
            'from_total_R':low['summary']['reflection_total'],'to_total_R':high['summary']['reflection_total']})
    write_csv(ROOT/'convergence_corrected.csv',convergence)
    write_json(ROOT/'convergence_corrected.json',convergence)
    final_name=RUNS[-1];final=runs[final_name];m=final['m']
    # field_amplitudes uses the RCWA cell's x_min as lateral phase origin.
    # The independent real-space field-monitor projection establishes this.
    xmin_to_global=np.exp(+1j*np.pi*m)
    ridge_center=final['centered']*xmin_to_global
    fdtd_factor=np.exp(-1j*np.pi*m)
    fdtd=ridge_center*fdtd_factor
    with np.load(ROOT/final_name/'field_amplitudes.npz') as a:
        ix=np.array([int(np.flatnonzero(a['n'].ravel()==order)[0]) for order in m])
        native_ey=a['eky_bs'][0,0,ix,0].copy()
    k0=2*np.pi/WAVELENGTH;kx=2*np.pi/PERIOD*m;kz=np.sqrt((k0*k0-kx*kx).astype(complex))
    z_factor=np.exp(1j*(k0+kz)*REFERENCE)
    np.savez_compressed(ROOT/'reference_for_analysis.npz',
        m=m,r_native_origin=fdtd,r_center_origin=ridge_center,
        r_rcwa_cell_origin=final['centered'],
        r_before_lateral_correction=final['r'],native_eky_bs=native_ey,
        native_s_input_sign=np.array(-1),reference_z_factor=z_factor,
        lateral_half_mesh_factor=final['correction'],rcwa_xmin_to_physical_origin_factor=xmin_to_global,
        fdtd_groove_origin_factor=fdtd_factor,
        R=final['R'],kx=kx,kz=kz,wavelength_m=np.array(WAVELENGTH),period_m=np.array(PERIOD),
        actual_rcwa_dx_m=np.array(final['dx']))
    table=[]
    for i,order in enumerate(m):
        table.append({'m':int(order),'angle_deg':float(np.degrees(np.arcsin(order*WAVELENGTH/PERIOD))),
            'r_center_real':ridge_center[i].real,'r_center_imag':ridge_center[i].imag,
            'r_fdtd_origin_real':fdtd[i].real,'r_fdtd_origin_imag':fdtd[i].imag,
            'abs_r':abs(fdtd[i]),'phase_fdtd_origin_deg':np.degrees(np.angle(fdtd[i])),
            'R':final['R'][i]})
    write_csv(ROOT/'final_reference.csv',table)
    metadata={'solver':'Ansys Lumerical RCWA','final_run':final_name,'total_harmonics':321,'actual_dx_nm':final['dx']*1e9,
        'scalar_field':'canonical Ey; FDTD native Ez','time_convention':'exp(-i omega t)',
        'phase_reference_z_nm':0,'incident_Ey_at_z0':1,
        'center_origin':'ridge center x=0; ridge occupies [-0.5,+0.5] um',
        'native_origin_name_note':'r_native_origin means FDTD native x origin, not raw RCWA native coefficients',
        'fdtd_native_origin':'even N finite mark centered as x_j=(j-(N-1)/2)*period; groove center x=0',
        'native_rcwa_lateral_phase_origin':'RCWA x_min=-1um; independently verified by physical-x field monitor Fourier projection',
        'rcwa_xmin_to_global_factor':'exp(-i*kxm*x_min)=exp(+i*m*pi)=(-1)^m',
        'fdtd_translation_m':PERIOD/2,'fdtd_origin_factor':'exp(-i*m*pi)=(-1)^m',
        'raw_to_center_formula':'r_center=-native_eky_bs*exp(i*(k0+kzm)*(-200nm))*exp(-i*kxm*dx_actual/2)*exp(+i*m*pi)',
        'raw_to_fdtd_formula':'r_fdtd=-native_eky_bs*exp(i*(k0+kzm)*(-200nm))*exp(-i*kxm*dx_actual/2); the two half-period factors cancel',
        'mesh_phase_evidence':'sampled index is mirror symmetric at x=0; raw r_m/r_-m=exp(i*m*G*dx_actual) for m=1,2,3; corrected pairs equal to numerical precision',
        'lateral_calibration_scope':'Measured v252 2D TE RCWA sampling convention, validated across 400,640,800,1600 mesh cells. Not a claim that the geometry is shifted.',
        'flat_fresnel_complex_error':json.loads((ROOT/'flat_calibration'/'summary.json').read_text())['fresnel_complex_error'],
        'max_phase_law_error':max(item['raw_pair_phase_law_error'] for run in runs.values() for item in run['audit']),
        'max_mirror_error_after_correction':max(item['corrected_pair_complex_error'] for run in runs.values() for item in run['audit']),
        'convergence_last_comparison':convergence[-1], 'all_original_datasets_preserved':True}
    write_json(ROOT/'reference_metadata.json',metadata)
    print('FINAL_TABLE',json.dumps(table),flush=True)
    print('CONVERGENCE',json.dumps(convergence),flush=True)
    print('AUDIT',json.dumps(metadata),flush=True)

if __name__=='__main__':main()
