"""Rebuild the complete RCWA reference, reusing verified saved artifacts.

From the repository root:
    <Lumerical Python executable> rcwa/reproduce_reference.py

No solver is launched when all required raw reference data already exist.
Use --dry-run to display required stages without executing them.
"""
from __future__ import annotations
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
RUNS=['flat_calibration','harmonics_020_dx5nm','harmonics_040_dx5nm',
      'harmonics_080_dx5nm','harmonics_080_dx2.5nm',
      'harmonics_080_dx1.25nm','harmonics_160_dx1.25nm']
CORE_FILES=['summary.json','model.fsp','native_datasets.mat',
            'reflection_propagating.npz','field_amplitudes.npz']

def execute(script,dry_run=False):
    print(('PLAN ' if dry_run else 'RUN  ')+script,flush=True)
    if not dry_run:
        subprocess.run([sys.executable,str(ROOT/script)],cwd=ROOT,check=True)

def verify():
    with np.load(ROOT/'reference_for_analysis.npz') as a:
        assert np.array_equal(a['m'],np.arange(-3,4))
        reconstructed=(-a['native_eky_bs']*a['reference_z_factor']*
            a['lateral_half_mesh_factor']*a['rcwa_xmin_to_physical_origin_factor']*
            a['fdtd_groove_origin_factor'])
        assert max(abs(reconstructed-a['r_native_origin']))<1e-13
        assert max(abs(a['r_native_origin']-a['r_native_origin'][::-1]))<1e-12
        power=np.real(a['kz'])/(2*np.pi/float(a['wavelength_m']))*abs(a['r_native_origin'])**2
        assert max(abs(power-a['R']))<1e-13
    origin=json.loads((ROOT/'spatial_origin_projection.json').read_text())
    assert origin['max_complex_error']<1e-10
    metadata=json.loads((ROOT/'reference_metadata.json').read_text())
    assert metadata['final_run']=='harmonics_160_dx1.25nm'
    assert metadata['total_harmonics']==321
    print('VERIFIED: all six grating cases, flat calibration, phase origins, and final coefficient/power adapter.',flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true',help='Show required stages without running them.')
    args=parser.parse_args()
    if not (ROOT/'material.json').exists():
        execute('probe_material.py',args.dry_run)
    else:print('REUSE material.json',flush=True)
    missing=[name for name in RUNS if not all((ROOT/name/file).exists() for file in CORE_FILES)]
    if missing:
        print('Missing complete solver cases: '+', '.join(missing),flush=True)
        execute('run_rcwa.py',args.dry_run)
    else:print('REUSE flat calibration and all six grating solver cases',flush=True)
    if any(not (ROOT/name/'actual_mesh.json').exists() or
           not (ROOT/name/'sampled_index.npz').exists() for name in RUNS[1:]):
        execute('inspect_mesh.py',args.dry_run)
    else:print('REUSE actual mesh and raw index exports',flush=True)
    execute('finalize_reference.py',args.dry_run)
    execute('audit_spatial_origin.py',args.dry_run)
    if not args.dry_run:verify()

if __name__=='__main__':main()
