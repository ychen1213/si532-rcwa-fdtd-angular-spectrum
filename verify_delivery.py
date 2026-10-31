"""Check experiment provenance and report/data consistency before publication."""
from pathlib import Path
import json,re,hashlib,csv
import numpy as np

ROOT=Path(__file__).resolve().parent

def main():
    checks=[]
    def check(name,ok,detail=None):
        checks.append(dict(name=name,passed=bool(ok),detail=detail))
    ref=np.load(ROOT/'rcwa'/'reference_for_analysis.npz')
    m=ref['m'];r=ref['r_native_origin']
    check('RCWA has exactly seven air propagating orders',np.array_equal(m,np.arange(-3,4)))
    check('RCWA mirror symmetry at common groove origin',max(abs(r-r[::-1]))<1e-10,float(max(abs(r-r[::-1]))))
    power=np.sqrt(1-(m*.532/2)**2)*abs(r)**2
    check('Passive RCWA reflection',0<sum(power)<1,float(sum(power)))
    metadata=json.loads((ROOT/'rcwa'/'reference_metadata.json').read_text())
    check('Independent complex Fresnel calibration',metadata['flat_fresnel_complex_error']<1e-10,metadata['flat_fresnel_complex_error'])
    origin=json.loads((ROOT/'rcwa'/'spatial_origin_projection.json').read_text())
    check('Independent physical-x RCWA monitor projection',origin['max_complex_error']<1e-10,origin['max_complex_error'])
    summaries=json.loads((ROOT/'analysis'/'case_summary.json').read_text())
    required_controls=['N020_dx10_dy5_t600','N020_dx10_dy5_t600_pml48',
      'N020_dx5_dy2.5_t600','N050_dx10_dy5_t600','N200_dx20_dy10_t1200',
      'N200_dx20_dy10_t1200_strictmark']
    names={t['case'] for t in summaries}
    check('All reported mesh, PML and strict time controls are analyzed',all(t in names for t in required_controls))
    strict=ROOT/'fdtd'/'cases'/'N200_dx20_dy10_t1200_strictmark'
    check('Strict time control explicitly reuses existing full-length air',
      (strict/'AIR_REFERENCE.json').exists() and not (strict/'air'/'DONE.json').exists())
    if (strict/'mark'/'DONE.json').exists():
        cfg=json.loads((strict/'mark'/'DONE.json').read_text())['settings']
        check('Strict mark uses 1200fs cap and 1e-12 shutoff',abs(cfg['t_max_s']-1.2e-12)<1e-20 and cfg['auto_shutoff']==1e-12)
    else:check('Strict mark actually completed',False)
    run_manifest=list(csv.DictReader(open(ROOT/'fdtd'/'solver_manifest.csv',encoding='utf-8-sig')))
    strict_rows=[t for t in run_manifest if t['case']=='N200_dx20_dy10_t1200_strictmark']
    check('Strict mark actually simulated at least 1200fs',bool(strict_rows) and float(strict_rows[0]['mark_simulated_time_s'])>=1.199e-12)
    finest=ROOT/'fdtd'/'cases'/'N020_dx5_dy2.5_t600'/'mark'/'observed_settings.json'
    check('Finest mesh reports effective rather than requested PML',finest.exists() and json.loads(finest.read_text())['pml layers']==64)
    actual=json.loads((ROOT/'fdtd'/'saved_project_audit.json').read_text())
    for tag in ['N010_dx20_dy10_t600','N400_dx20_dy10_t600']:
        a=next(t for t in actual if t['case']==tag)
        check(tag+' actual saved four boundaries are PML',all(a['FDTD'][k]=='PML' for k in ['x min bc','x max bc','y min bc','y max bc']))
        check(tag+' actual saved dimension is 2D',a['FDTD']['dimension']=='2D')
    for n in [10,20,50,100,200,400]:
        tag=f'N{n:03d}_dx20_dy10_t600'
        pair=ROOT/'fdtd'/'cases'/tag
        check(tag+' actual air and mark results',all((pair/t/'DONE.json').exists() for t in ['air','mark']))
        for kind in ['air','mark']:
            done=json.loads((pair/kind/'DONE.json').read_text())
            check(tag+'/'+kind+' engine completed successfully',done['engine']['returncode']==0 and done['status']=='completed')
            cfg=done['settings']
            check(tag+'/'+kind+' finite PML plane-wave settings',cfg['boundary']=='all PML' and cfg['N']==n and cfg['duty']==.5)
            z=np.load(pair/kind/'fields_R.npz')
            check(tag+'/'+kind+' actual complex R field at 532nm',np.iscomplexobj(z['Ez']) and np.all(np.isfinite(z['Ez'])) and abs(float(z['frequency_hz'].squeeze())-299792458/532e-9)<1)
    table=list(csv.DictReader(open(ROOT/'analysis'/'all_extracted_orders.csv',encoding='utf-8')))
    check('Exactly seven extracted rows per complete case',len(table)==7*len(summaries),len(table))
    check('Every extracted amplitude retains real and imaginary components',all(all(np.isfinite(float(t[k])) for k in ['fdtd_re','fdtd_im','rcwa_re','rcwa_im']) for t in table))
    md=(ROOT/'report'/'REPORT_zh.md').read_text(encoding='utf-8')
    best=min(summaries,key=lambda s:s['relative_complex_L2'])['case']
    readme=(ROOT/'README.md').read_text(encoding='utf-8')
    check('Report and README use the current actual best case',best in md and best in readme,best)
    check('No accidental non-whitespace control characters in report',not any(ord(c)<32 and c not in '\n\r\t' for c in md))
    images=re.findall(r'!\[[^\]]*\]\(([^)]+)\)',md)
    check('All report figures exist',len(images)==12 and all((ROOT/'report'/p).exists() for p in images),len(images))
    check('Every figure has purpose and conclusion',md.count('**图的目的：**')==12 and md.count('**图的结论：**')==12)
    check('Report explicitly distinguishes scalar sinc energy from weighted power','标量谱能量' in md and 'Poynting' in md)
    check('Report generated from all six advertised sizes',all(any(t['N']==n and t['case']==f'N{n:03d}_dx20_dy10_t600' for t in summaries) for n in [10,20,50,100,200,400]))
    check('PDF and offline HTML exist',all((ROOT/'report'/p).exists() for p in ['REPORT_zh.pdf','REPORT_zh.html']))
    html=(ROOT/'report'/'REPORT_zh.html').read_text(encoding='utf-8')
    check('Offline HTML embeds all figures',html.count('data:image/png;base64,')>=12,html.count('data:image/png;base64,'))
    check('Offline HTML includes MathML formulas','<math' in html)
    result={'checks':checks,'passed':sum(t['passed'] for t in checks),'total':len(checks),'all_passed':all(t['passed'] for t in checks)}
    (ROOT/'analysis'/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
    if not result['all_passed']:
        print(json.dumps([t for t in checks if not t['passed']],ensure_ascii=False,indent=2))
        raise SystemExit(1)

if __name__=='__main__':main()
