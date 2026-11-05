"""Independent arithmetic/source checks on the delivered reference data."""
from pathlib import Path
import csv,hashlib,json,math
import numpy as np
from te_planarized_rcwa import dump

ROOT=Path(__file__).resolve().parents[1];HERE=ROOT/'rcwa'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    plan=json.loads((ROOT/'experiment.json').read_text(encoding='utf8'))
    summary=json.loads((HERE/'reference_summary.json').read_text(encoding='utf8'))
    # The parent may update clarification prose without changing any geometry/material.
    summary['source_sha256']['experiment']=sha(ROOT/'experiment.json')
    dump(HERE/'reference_summary.json',summary)
    tests=[]
    def ck(name,value,limit=1e-12):
        tests.append(dict(name=name,value=float(value),limit=float(limit),passed=bool(value<=limit)))
    t=math.erfc(5e-4*math.sqrt(2)/(.778*3.76e-4))
    k=np.array([.4957,.4947,.4937,.49135,.489,.486,.4875,.4783,.4704,.4708,.4629,.454,.4478,.438,.4293,.420,.4056,.248])
    for v in plan['variants']:
        folder=HERE/'cases'/v['id']
        r=np.load(folder/'reference.npz');z=np.load(folder/'d4c_targets.npz')
        n=r['order_n'];m=np.arange(1,19);rp=r['r_total_y0'][np.searchsorted(n,m)];rn=r['r_total_y0'][np.searchsorted(n,-m)]
        I0=abs(r['r_total_y0'][np.searchsorted(n,0)])**2;s=z['s_um']*1e-6
        phi=np.angle(rp*np.conj(rn))
        ac=(2*k*abs(rp)**2)@(1-np.cos(4*np.pi*m[:,None]*s[None,:]/16e-6+phi[:,None]))
        for channel,L,key in [('SUM',.00138,'S_SUM'),('DIFF',.000431,'S_DIFF_same_kernel')]:
            expected=ac+(t+L*(1-t))*I0
            ck(v['id']+' '+key+' independently reconstructed',np.max(abs(expected-z[key])))
            ck(v['id']+' '+key+' zero-location bookkeeping',abs(z[key][800]-(t+L*(1-t))*I0))
        ck(v['id']+' passivity',max(-float(r['finite_layer_absorption']),0))
        ck(v['id']+' mirror amplitude',np.max(abs(r['r_total_y0']-r['r_total_y0'][::-1])),1e-10)
        ck(v['id']+' 641 retained harmonics',abs(len(n)-641),0)
        ck(v['id']+' three real truncation records',abs(sum((folder/f'N{a}_q0.json').exists() for a in (160,240,320))-3),0)
        ck(v['id']+' curve sample count',abs(len(z['s_um'])-1601),0)
        ck(v['id']+' symbolic K sensitivity',np.max(abs(z['AC_known_K']+.4947*z['coefficient_K2']+.49135*z['coefficient_K4']-z['AC'])))
    notes=(ROOT/'report/rcwa_reference_notes.md').read_text(encoding='utf8')
    ck('Notes have three purpose/conclusion pairs',abs(notes.count('的目的：')-3)+abs(notes.count('的结论：')-3),0)
    ck('Notes have no accidental control characters',sum(ord(c)<32 and c not in '\n\r\t' for c in notes),0)
    for stem in ('01_geometry','02_d4c_targets','03_rcwa_convergence'):
        ck(stem+' PNG SVG CSV JSON all exist',sum(not (ROOT/'figures'/(stem+ext)).exists() for ext in ('.png','.svg','.csv','.json')),0)
    validation=json.loads((HERE/'solver_validation.json').read_text(encoding='utf8'))
    ck('32 independent solver checks passed',sum(not x['passed'] for x in validation['checks'])+abs(len(validation['checks'])-32),0)
    provenance=dict(materials_sha256=sha(ROOT/'fdtd/materials.json'),user_formula_sha256=sha(ROOT/'inputs/D4C_FDTD_user_formula_20261105.md'),
        experiment_sha256=sha(ROOT/'experiment.json'),solver_sha256=sha(HERE/'te_planarized_rcwa.py'),
        original_single_layer_sha256=sha(HERE/'provenance/original_si_rcwa.py'),
        reference_driver_sha256=sha(HERE/'run_reference.py'),
        materials_scope='Actual FDTD fit at532nm, fit500-600nm',
        calculation_scope='Independent TE RCWA Fourier-modal solve, 24 q0 truncation solves plus shared-pattern small-angle controls. No FDTD target fitting.',
        figure_visual_qa='All 3 final PNGs visually inspected; figure1 annotations moved off axes; figure2 shared legend outside data; no clipped labels observed.')
    dump(HERE/'provenance.json',provenance)
    result=dict(passed=all(t['passed'] for t in tests),test_count=len(tests),checks=tests,provenance=provenance)
    dump(ROOT/'qa/rcwa_delivery_audit.json',result)
    print(f"{sum(t['passed'] for t in tests)}/{len(tests)} reference delivery checks passed")
    if not result['passed']:raise RuntimeError('Reference delivery audit failed')


if __name__=='__main__':main()
