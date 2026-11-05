"""Final provenance/completeness and cross-file scientific consistency audit.

This is a delivery audit, not a second scan solver. It checks real completed
engine records, saved settings/raw fields, locked train/heldout separation,
and invariant zero-accounting differences between saved candidate curves.
--allow-partial records missing artifacts as PENDING and never calls them passed.
"""
from __future__ import annotations
from pathlib import Path
import argparse,csv,hashlib,json,math,re,sys
from datetime import datetime,timezone
from urllib.parse import unquote,urlparse
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
COMPS=('Ex','Ey','Ez','Hx','Hy','Hz')
T0=math.erfc(5e-4*math.sqrt(2)/(.778*3.76e-4))
LEAK={'SUM':.00138,'DIFF':.000431}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def readjson(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def readcsv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


class Audit:
    def __init__(self,allow_partial=False):self.items=[];self.stats={};self.hash_cache={};self.allow_partial=allow_partial
    def check(self,name,ok,detail=None):
        self.items.append(dict(check=name,status='PASS' if bool(ok) else 'FAIL',detail=detail))
    def pending(self,name,detail=None):self.items.append(dict(check=name,status='PENDING',detail=detail))
    def need(self,path,name=None):
        path=Path(path)
        if not path.is_file():self.pending(name or 'Missing '+str(path.relative_to(ROOT)),str(path));return False
        return True
    def close(self,name,a,b,atol=1e-12,rtol=1e-10):
        a=np.asarray(a);b=np.asarray(b)
        ok=a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and np.allclose(a,b,atol=atol,rtol=rtol)
        err=float(np.max(np.abs(a-b))) if a.shape==b.shape and a.size else None
        self.check(name,ok,dict(max_abs_error=err,atol=atol,rtol=rtol))
    def hash(self,p):
        p=Path(p)
        if p not in self.hash_cache:self.hash_cache[p]=sha(p)
        return self.hash_cache[p]
    def group(self,name,fn):
        try:fn()
        except Exception as exc:self.check(name+' can be inspected',False,type(exc).__name__+': '+str(exc))


def inspect_case(a,name,variant=None):
    folder=ROOT/'fdtd/cases'/name
    if not a.need(folder/'DONE.json',name+' genuine completed solve'):return None
    d=readjson(folder/'DONE.json');s=d.get('settings',{});t=d.get('termination',{})
    a.check(name+' DONE engine success',d.get('status')=='completed' and t.get('completed_successfully') is True)
    for file in ('settings.json','observed_settings.json','engine_result.json','model.fsp','build_this_case.lsf','engine_stdout_public.log'):
        a.need(folder/file,name+' '+file)
    if (folder/'engine_result.json').exists():
        e=readjson(folder/'engine_result.json');a.check(name+' actual engine return code',e.get('returncode',e.get('return_code'))==0)
    if (folder/'settings.json').exists():a.check(name+' settings agree with DONE',readjson(folder/'settings.json')==s)
    if (folder/'model.fsp').exists():a.check(name+' saved model matches DONE hash',a.hash(folder/'model.fsp')==d.get('model_sha256'))
    a.close(name+' target wavelength',s.get('wavelength_m',0),532e-9,atol=1e-15)
    a.close(name+' fine dx/dy',[s.get('dx_m',0),s.get('dy_m',0)],[10e-9,5e-9],atol=1e-15)
    a.close(name+' Gaussian waist',s.get('waist_radius_m',0),19e-6,atol=1e-14)
    actual=float(t.get('simulated_time_s',0));maximum=float(s.get('t_max_s',0))
    a.check(name+' positive actual integration duration',0<actual<=maximum*1.001,dict(actual_s=actual,requested_s=maximum))
    if t.get('terminated_early'):
        a.check(name+' early termination met stated shutoff',float(t.get('final_auto_shutoff',math.inf))<=float(s.get('auto_shutoff',0))*1.1)
    if (folder/'observed_settings.json').exists():
        obs=readjson(folder/'observed_settings.json')
        a.check(name+' saved four actual boundaries are PML',all(obs.get(k)=='PML' for k in ('x min bc','x max bc','y min bc','y max bc')))
        a.check(name+' saved actual 2D uniform mesh',obs.get('dimension')=='2D' and obs.get('mesh type')=='uniform')
        a.check(name+' actual PML layer count',int(obs.get('pml layers',0))==int(s.get('pml_layers',-1))==64,
            dict(actual=obs.get('pml layers'),requested=s.get('pml_layers')))
        a.close(name+' source settings',[obs.get('source',{}).get(k,0) for k in ('waist radius w0','amplitude','phase','polarization angle')],[19e-6,1.,0.,90.],atol=1e-13)
        if variant:
            a.check(name+' same intended variant',s.get('geometry')==variant)
            geom=obs.get('geometry',{})
            a.check(name+' actual ridge count',geom.get('ridge_count')==variant['periods']*len(variant['ridge_centers_local_um']))
            a.check(name+' actual oxide presence',geom.get('oxide_present')==(variant['cap_nm']>0))
            if variant['cap_nm']>0:
                ox=geom.get('SiO2_planarization',{})
                a.close(name+' actual planarized oxide height',[ox.get('y min',-1),ox.get('y max',-1)],[0.,(variant['ridge_height_nm']+variant['cap_nm'])*1e-9],atol=1e-14)
            ridge=geom.get('Si_ridge_00_0',{})
            a.close(name+' actual ridge height/width',[ridge.get('y max',0)-ridge.get('y min',0),ridge.get('x span',0)],
                [variant['ridge_height_nm']*1e-9,variant['ridge_width_um']*1e-6],atol=1e-14)
        mesh=obs.get('mesh_actual',{})
        for axis,target in (('x_m',10e-9),('y_m',5e-9)):
            z=mesh.get(axis,{})
            a.close(name+' actual '+axis+' spacing',[z.get('step_min_m',0),z.get('step_max_m',0)],[target,target],atol=1e-14)
    for mon in ('R','R2','INC','T'):
        p=folder/f'nearfield_{mon}.npz'
        if not a.need(p,name+' raw '+mon+' field'):continue
        with np.load(p) as z:
            a.check(name+' '+mon+' all six complex raw components',all(k in z.files and np.iscomplexobj(z[k]) for k in COMPS))
            x=z['x_m'].squeeze()
            a.check(name+' '+mon+' finite monotonic raw samples',len(x)>100 and np.all(np.diff(x)>0) and all(z[k].size==x.size and np.isfinite(z[k]).all() for k in COMPS))
            a.close(name+' '+mon+' 532nm recorded',z['wavelength_m'],532e-9,atol=1e-15)
            a.check(name+' '+mon+' Ez is nonzero',np.linalg.norm(z['Ez'])>0)
    airfile=folder/'AIR_REFERENCE.json'
    airname=readjson(airfile).get('air_reference') if airfile.exists() else 'air_fine'
    ref=ROOT/'fdtd/cases'/str(airname)/'DONE.json'
    if name!='air_fine' and a.need(ref,name+' completed same-source air reference'):
        rs=readjson(ref)['settings']
        keys=('dx_m','dy_m','pml_layers','x_span_m','monitor_span_m','y_min_m','y_max_m','source_y_m','source_span_m','waist_radius_m','source_band_m','monitor_y_m')
        a.check(name+' air reference has same grid and source',all(s.get(k)==rs.get(k) for k in keys))
        if airfile.exists():a.check(name+' explicit air reuse is not a duplicate fake solve',readjson(airfile).get('reference_reused') is True and airname!=name)
    return d


def inspect_spectrum(a,name,suffix='_fine'):
    p=ROOT/'outputs/spectra'/(name+suffix+'.npz');mp=p.with_suffix('.json')
    if not a.need(p,name+' processed fine complex spectrum') or not a.need(mp,name+' spectrum provenance'):return
    meta=readjson(mp)
    with np.load(p) as z:
        required=('kx_rad_m','A_ref_V','A_ref_R2_V','A_ref_EH_V','Ainc_ref_V','Ainc_y0_V','A_zero_degree0_V','A_zero_degree2_V','Pin_W_m')
        a.check(name+' all spectrum/decomposition arrays',all(k in z.files for k in required))
        k=z['kx_rad_m'];a.check(name+' symmetric angular grid',np.all(np.diff(k)>0) and np.allclose(k,-k[::-1],atol=1e-5,rtol=1e-10))
        a.check(name+' finite complex spectra',all(z[key].shape==k.shape and np.isfinite(z[key]).all() and np.iscomplexobj(z[key]) for key in required if key.startswith('A')))
        a.check(name+' fixed measured positive Pin',float(z['Pin_W_m'])>0 and np.isclose(float(z['Pin_W_m']),meta['Pin_W_m']))
        a.close(name+' reference plane',z['reference_y_m'],.85e-6,atol=1e-15)
    for label,item in meta.get('monitors',{}).items():
        src=(ROOT/item['file']).resolve()
        a.check(name+' '+label+' monitor path stays in project',src.is_relative_to(ROOT.resolve()))
        if a.need(src,name+' '+label+' actual raw source'):
            a.check(name+' '+label+' raw SHA256 matches',a.hash(src)==item.get('sha256'))
    for deg in (0,2):
        pm=meta.get(f'projection_degree{deg}',{})
        a.check(name+f' degree{deg} source-only projection provenance','empty-domain FDTD' in pm.get('source','') and 'no RCWA' in pm.get('source',''))
        a.check(name+f' degree{deg} residual recorded without concealing accuracy',np.isfinite(pm.get('relative_field_L2',math.nan)) and pm.get('columns',0)>=pm.get('rank',math.inf)>0)


def inspect_rcwa(a,variants):
    p=ROOT/'qa/rcwa_delivery_audit.json'
    if a.need(p,'Independent RCWA delivery audit'):
        q=readjson(p);a.check('Existing RCWA delivery audit passed all 86 checks',q.get('passed') is True and q.get('test_count')==86 and all(x['passed'] for x in q['checks']))
    p=ROOT/'rcwa/solver_validation.json'
    if a.need(p,'Thin-film/solver physics validation'):
        q=readjson(p);a.check('32 RCWA physics checks passed',q.get('passed') is True and len(q['checks'])==32 and all(x['passed'] for x in q['checks']))
    table={int(r['m']):r for r in readcsv(ROOT/'rcwa/K_table.csv')}
    a.close('User-authorized K2/K4 interpolation',[float(table[2]['K']),float(table[4]['K'])],[.4947,.49135],atol=1e-14)
    a.close('New supplied K13/K14 preserved',[float(table[13]['K']),float(table[14]['K'])],[.4478,.4380],atol=1e-14)
    mat=readjson(ROOT/'fdtd/materials.json')['materials']
    for v in variants:
        folder=ROOT/'rcwa/cases'/v['id'];p=folder/'reference.npz';tp=folder/'d4c_targets.npz'
        if not a.need(p,v['id']+' RCWA reference') or not a.need(tp,v['id']+' complete D4C target'):continue
        with np.load(p) as r,np.load(tp) as z:
            a.check(v['id']+' primary complete 1..18 target',all(x in z.files for x in ('S_SUM','S_DIFF_same_kernel','AC','T_m','r_plus','r_minus')) and len(z['T_m'])==18)
            a.close(v['id']+' confirmed Tm equals field modulus squared',z['T_m'],abs(z['r_plus'])**2)
            a.close(v['id']+' supplied T0',z['T0'],T0)
            for channel,key in (('SUM','S_SUM'),('DIFF','S_DIFF_same_kernel')):
                expected=(T0+LEAK[channel]*(1-T0))*float(z['I0'])
                a.close(v['id']+' D4C '+channel+' zero counted once',z[key]-z['AC'],np.full_like(z[key],expected))
            cfg=json.loads(str(r['config_json']))
            for key,material in (('ridge_n','Si'),('substrate_n','Si'),('cap_n','SiO2')):
                ref=mat[material]['FDTD_fit_index']
                a.close(v['id']+' matched '+key,[cfg[key]['real'],cfg[key]['imag']],[ref['n'],ref['k']],atol=1e-13)
            a.check(v['id']+' full final RCWA vector',len(r['order_n'])==641 and np.isfinite(r['r_total_y0']).all())


def inspect_comparison(a,variants):
    folder=ROOT/'outputs/comparison_fine';expected={v['id'] for v in variants}
    paths=[folder/x for x in ('summary.json','scan_curves.npz','candidate_errors.csv','zero_power_audit.csv','shared_fits.json','shared_fit_errors.csv','shared_stop_sweep.csv')]
    if not all([a.need(p,'Final comparison '+p.name) for p in paths]):return
    summary=readjson(paths[0])
    ready=set(summary.get('available_structures',[]))==expected and summary.get('all_eight_complete') is True
    if not ready:
        detail=dict(available=summary.get('available_structures',[]),expected=sorted(expected))
        if a.allow_partial:a.pending('Final eight-structure comparison and held-out fits not ready',detail)
        else:a.check('Final comparison contains all eight completed structures',False,detail)
        a.stats['comparison_structures']=len(summary.get('available_structures',[]))
        return
    a.check('Final comparison contains all eight completed structures',True)
    a.check('No per-Mark offsets or fitted per-order gains in final comparison',summary.get('no_per_mark_offsets') is True and summary.get('no_per_order_fitted_gain') is True)
    a.check('All final power integrals are angular / infinite image','angular' in summary.get('all_integrals','') and 'infinite image' in summary.get('all_integrals',''))
    candidates=readcsv(paths[2]);zeros=readcsv(paths[3]);fits=readjson(paths[4]);fitrows=readcsv(paths[5])
    train={v['id'] for v in variants if v['cap_nm'] in (0,100)};hold=expected-train
    trainpairs={(i,ch) for i in train for ch in LEAK};holdpairs={(i,ch) for i in hold for ch in LEAK}
    with np.load(paths[1]) as z:
        s=z['s_um'];a.close('Final scan uses fixed -4..4um 0.005um grid',s,np.arange(-800,801)*.005,atol=1e-12)
        a.check('All stored scan vectors have full finite shape',all(z[k].shape==s.shape and np.isfinite(z[k]).all() for k in z.files if k!='s_um'))
        for v in variants:
            name=v['id']
            with np.load(ROOT/'rcwa/cases'/name/'d4c_targets.npz') as target:
                for ch,key in (('SUM','S_SUM'),('DIFF','S_DIFF_same_kernel')):
                    ck=f'{name}__{ch}__D4C'
                    a.check(name+' '+ch+' target vector present',ck in z.files)
                    if ck in z.files:a.close(name+' '+ch+' final curve target equals original RCWA target',z[ck],target[key])
        for row in candidates:
            key=f"{row['id']}__{row['channel']}__{row['model']}"
            a.check(key+' saved actual comparison curve',key in z.files)
        # These differences isolate zero accounting without reproducing the B solver.
        for name in expected:
            lookup={r['model']:r for r in zeros if r['id']==name}
            needed=('weighted_matchedstop_mode2_scalar','weighted_matchedstop_mode2_physical','weighted_matchedstop_mode2_double','weighted_nostop_mode2_scalar')
            if not all(x in lookup for x in needed):a.check(name+' all zero-accounting controls saved',False);continue
            zz=lookup['weighted_matchedstop_mode2_scalar'];raw=float(zz['P0_over_Pin']);passed=float(zz['P0pass_over_Pin'])
            a.check(name+' measured zero transmission has valid raw/pass power',raw>0 and 0<=passed<=raw*(1+1e-10))
            for ch in LEAK:
                prefix=name+'__'+ch+'__';sc=z[prefix+'weighted_matchedstop_mode2_scalar']
                physical=z[prefix+'weighted_matchedstop_mode2_physical'];double=z[prefix+'weighted_matchedstop_mode2_double']
                a.close(name+' '+ch+' repeated model differs by exactly passed P0',double-sc,np.full_like(sc,passed),atol=1e-12)
                delta=(1-LEAK[ch])*(passed-T0*raw)
                a.close(name+' '+ch+' physical-vs-scalar once accounting invariant',physical-sc,np.full_like(sc,delta),atol=1e-12)
            raw=float(lookup['weighted_nostop_mode2_scalar']['P0_over_Pin'])
            a.close(name+' scalar channels differ only by supplied leakage',z[name+'__SUM__weighted_nostop_mode2_scalar']-z[name+'__DIFF__weighted_nostop_mode2_scalar'],
                np.full_like(s,(LEAK['SUM']-LEAK['DIFF'])*(1-T0)*raw),atol=1e-12)
        a.check('Expected structured shared-fit families present',all(key in fits for key in ('decband_K_over_gamma_zero_transmission','decband_K_over_gamma_gain+zero_transmission','dec2_K_over_gamma_zero_transmission','dec2_K_over_gamma_gain+zero_transmission')))
        for key,fit in fits.items():
            pairs={(r['id'],r['channel']) for r in fit.get('records',[])}
            a.check(key+' trains only on both geometries cap0/cap100',pairs==trainpairs)
            a.check(key+' one parameter set, no per-Mark offsets',fit.get('shared_across_marks') is True and fit.get('per_mark_offsets') is False and set(fit['parameters'])<=set(('gain','zero_transmission','leak_sum','leak_diff')))
            a.check(key+' fit is numerically identifiable',fit.get('success') is True and fit.get('jacobian_rank')==len(fit.get('free',[])))
            held=fit.get('heldout') or {}
            a.check(key+' held-out both geometries cap150/cap200',{(r['id'],r['channel']) for r in held.get('records',[])}==holdpairs)
            a.check(key+' held-out parameters are frozen',held.get('parameters')==fit['parameters'] and held.get('free')==[])
            rr=[r for r in fitrows if r['fit']==key]
            a.check(key+' per-case split table has all sixteen records',len(rr)==16 and {(r['id'],r['channel']) for r in rr}==trainpairs|holdpairs and all(r['split']==('train' if r['id'] in train else 'heldout') for r in rr))
            a.check(key+' all saved fitted curves present',all(f"{i}__{ch}__fit_{key}" in z.files for i,ch in trainpairs|holdpairs))
    sweep=readcsv(paths[6]);best=summary.get('best_shared_physical_stop',{})
    if sweep:
        selected=min(sweep,key=lambda r:float(r['train_mean_relative_L2_pct']))
        a.close('One common stop selected only by training score',best.get('h_rad_um',math.nan),float(selected['h_rad_um']),atol=1e-12)
    a.stats['shared_fit_count']=len(fits);a.stats['comparison_structures']=len(summary.get('available_structures',[]))


def inspect_reports(a,explicit=None):
    folder=ROOT/'report'
    if explicit:
        main=Path(explicit);main=main if main.is_absolute() else ROOT/main
    else:
        candidates=[p for p in folder.glob('*.md') if p.name!='rcwa_reference_notes.md']
        main=max(candidates,key=lambda p:p.stat().st_size) if candidates else folder/'REPORT_zh.md'
    if not a.need(main,'Final main report (not only reference notes)'):return
    files=list(folder.glob('*.md'))+list(folder.glob('*.html'))
    if (ROOT/'README.md').exists():files.append(ROOT/'README.md')
    linked=0;figs=set();main_figs=set();main_embedded_images=0
    for file in files:
        text=file.read_text(encoding='utf-8-sig')
        if file.resolve()==main.resolve():main_embedded_images=len(re.findall(r'<img[^>]+src=[\"\']data:image',text,re.I))
        a.check(file.name+' has no accidental control characters',not any(ord(c)<32 and c not in '\n\r\t' for c in text))
        links=re.findall(r'!?\[[^\]\n]*\]\(([^\n]*?)\)',text)+re.findall(r'(?:src|href)=[\"\']([^\"\']+)[\"\']',text)
        for link in links:
            target=link.strip().strip('<>')
            if not target or target.startswith(('#','data:','mailto:')):continue
            if urlparse(target).scheme and not re.match(r'^[A-Za-z]:[\\/]',target):continue
            target=unquote(target.split('#')[0].split('?')[0])
            target=re.sub(r':\d+$','',target)
            # Pandoc links may carry a quoted title following the path.
            target=re.sub(r'\s+[\"\'].*[\"\']$','',target)
            p=Path(target);p=p if p.is_absolute() else file.parent/p
            a.check(file.name+' local link resolves: '+target,p.exists())
            linked+=1
            if p.suffix.lower() in ('.png','.svg','.jpg','.jpeg'):
                figs.add(str(p.resolve()))
                if file.resolve()==main.resolve():main_figs.add(str(p.resolve()))
    a.check('Final main report includes actual result figures',len(main_figs)+main_embedded_images>=3,
        dict(main_unique_local_images=len(main_figs),main_embedded_images=main_embedded_images,all_linked_images=len(figs)))
    a.stats['report']=str(main.relative_to(ROOT));a.stats['local_links_checked']=linked;a.stats['linked_figures']=len(figs)
    pdf=main.with_suffix('.pdf')
    if pdf.exists():a.check('Final PDF is nonempty',pdf.stat().st_size>10000)


def inspect_physical_component(a,variants):
    """Check the later gain-on-nonzero-only branch and its same held-out split."""
    folder=ROOT/'outputs/physical_component_fine'
    if not all([a.need(folder/f,'Physical-component '+f) for f in ('summary.json','shared_fits.json','scan_curves.npz')]):return
    s=readjson(folder/'summary.json');fits=readjson(folder/'shared_fits.json')
    expected={v['id'] for v in variants};train={v['id'] for v in variants if v['cap_nm'] in (0,100)};hold=expected-train
    a.check('Physical-component comparison retains original eight only',set(s.get('available_structures',[]))==expected and s.get('all_eight_complete') is True)
    a.check('Physical-component training/heldout sets remain original',set(s.get('training_structures',[]))==train and set(s.get('heldout_structures',[]))==hold)
    for key in ('actual_nearest_band_NA0.60__nonzero_only__gain_ac','actual_nearest_band_NA0.60__nonzero_only_baseline_aware__gain_ac+zero_transmission'):
        if key not in fits:a.check(key+' main branch exists',False);continue
        f=fits[key];held=f.get('heldout',{})
        a.check(key+' gain is explicitly nonzero-only',f.get('gain_scope')=='nonzero_only' and 'gain_ac' in f.get('parameters',{}))
        a.check(key+' no per-mark offsets',f.get('per_mark_offsets') is False and f.get('shared_across_marks') is True)
        a.check(key+' training records remain cap0/100',{r['id'] for r in f.get('records',[])}==train)
        a.check(key+' heldout records remain cap150/200',{r['id'] for r in held.get('records',[])}==hold)
        a.check(key+' heldout parameters frozen',held.get('parameters')==f.get('parameters') and held.get('free')==[])
    for item in s.get('inputs',[]):
        name=item['id']
        for key,path in [('spectrum_sha256',ROOT/'outputs/spectra'/(name+'_fine.npz')),
            ('projection_metadata_sha256',ROOT/'outputs/spectra'/(name+'_fine.json')),
            ('target_curve_sha256',ROOT/'rcwa/cases'/name/'d4c_targets.npz')]:
            if a.need(path,'Physical main input '+name+' '+key):a.check('Physical main current hash '+name+' '+key,a.hash(path)==item.get(key))


def inspect_extra(a,variants):
    """Ninth structure must be independently held out from every original fit."""
    inputfile=ROOT/'inputs/extra_weak_mark.json'
    if not inputfile.exists():return
    v=readjson(inputfile)['variant'];ID=v['id']
    inspect_case(a,ID+'_fine_t200',v);inspect_spectrum(a,ID,suffix='_fine_t200')
    a.check('Ninth geometry remains separate from original variants',ID not in {x['id'] for x in variants})
    a.close('Ninth geometry h180/cap400',[v['ridge_height_nm'],v['cap_nm']],[180.,400.])
    folder=ROOT/'outputs/extra_weak';paths=[folder/x for x in ('summary.json','scan_curves.npz','scan_curves.csv','errors.csv')]
    if not all([a.need(p,'Extra heldout '+p.name) for p in paths]):return
    summary=readjson(paths[0]);target=ROOT/'rcwa/cases'/ID/'d4c_targets.npz'
    if summary.get('primary_case_suffix')!='_fine_t200':
        if a.allow_partial:a.pending('Weak primary analysis must be refreshed to the 200fs case')
        else:a.check('Weak primary analysis uses 200fs after 100fs failed',False)
        return
    a.check('Weak primary analysis uses 200fs after 100fs failed',True)
    if not a.need(target,'Extra final D4C target'):return
    a.check('Extra mark has no training additions',summary.get('training_additions') is False and summary.get('extra_mark_was_not_used_for_any_fit') is True and summary.get('extra_independent_holdout') is True)
    a.check('Extra output geometry equals its independent input',summary.get('structure',{}).get('variant')==v)
    cal=(ROOT/summary['frozen_calibration_file']).resolve()
    a.check('Extra calibration points to original physical-component fits',cal== (ROOT/'outputs/physical_component_fine/shared_fits.json').resolve())
    if not a.need(cal,'Frozen original calibration file'):return
    fits=readjson(cal)
    a.check('Extra frozen calibration SHA is current',a.hash(cal)==summary.get('calibration_sha256',summary.get('frozen_calibration_sha256')))
    pairs=[('frozen_nonzero_gain_parameters','actual_nearest_band_NA0.60__nonzero_only__gain_ac'),
        ('frozen_baseline_aware_parameters','actual_nearest_band_NA0.60__nonzero_only_baseline_aware__gain_ac+zero_transmission')]
    for field,key in pairs:a.check('Extra '+field+' exactly equals pre-existing fit',summary.get(field)==fits[key]['parameters'])
    oldfile=ROOT/'outputs/comparison_fine/shared_fits.json';old=readjson(oldfile)
    a.check('Extra failed curve-only example also freezes its old parameters',summary.get('frozen_failed_curve_only_parameters')==old['dec2_K_over_gamma_gain+zero_transmission']['parameters'])
    # No identifier of the ninth mark is allowed in training OR original holdout fits.
    for label,collection in [('original',old),('physical-component',fits)]:
        ids=[]
        for f in collection.values():
            ids.extend(r['id'] for r in f.get('records',[]));ids.extend(r['id'] for r in (f.get('heldout') or {}).get('records',[]))
        a.check('Ninth mark absent from every '+label+' fit record',ID not in ids)
    spectral=ROOT/'outputs/spectra'/(ID+'_fine_t200.npz')
    for field,path in [('target_curve_sha256',target),('spectrum_sha256',spectral)]:
        a.check('Extra provenance records current '+field,summary.get(field)==a.hash(path))
    conv=readcsv(ROOT/'rcwa/cases'/ID/'convergence.csv');last=conv[-1]
    a.check('Weak extra reference nonzero convergence below 0.05%',float(last['nonzero_relative_complex_L2'])<.0005)
    a.stats['extra_reference_max_order']=int(last['max_order'])
    raw=float(summary['physical_zero_over_Pin']);passed=float(summary['physical_zero_after_stop_over_Pin'])
    a.check('Extra actual zero raw/passed are sensible',0<=passed<=raw and raw>0)
    with np.load(paths[1]) as z,np.load(target) as t:
        s=z['s_um'];a.close('Extra scan uses same fixed displacement grid',s,t['s_um'])
        a.check('Every extra stored scan vector is full finite data',all(z[k].shape==s.shape and np.isfinite(z[k]).all() for k in z.files))
        a.check('Extra target peak is 0.00x and nonzero bottom',.001<=float(t['S_SUM'].max())<=.009 and float(t['S_SUM'].min())>0)
        for ch,key in (('SUM','S_SUM'),('DIFF','S_DIFF_same_kernel')):
            a.close('Extra '+ch+' original q0 target unchanged',z['D4C_'+ch],t[key])
            fixed=z['fixed_T0_once_'+ch];physical=z['physical_stop_once_'+ch];double=z['repeated_passed_zero_'+ch]
            diff=(1-LEAK[ch])*T0*raw+LEAK[ch]*passed
            a.close('Extra '+ch+' repeated-vs-physical zero accounting',double-physical,np.full_like(s,diff),atol=1e-12)
            for model,field in [('frozen_nonzero_gain','frozen_nonzero_gain_parameters'),('frozen_baseline_aware','frozen_baseline_aware_parameters')]:
                par=summary[field];gain=par['gain_ac'];tau=par['zero_transmission']
                c=lambda x:x+LEAK[ch]*(1-x)
                a.close('Extra '+ch+' '+model+' only scales nonzero contribution',z[model+'_'+ch]-gain*fixed,
                    np.full_like(s,(c(tau)-gain*c(T0))*raw),atol=1e-12)
        a.close('Extra scalar channels differ only by leakage',z['fixed_T0_once_SUM']-z['fixed_T0_once_DIFF'],
            np.full_like(s,(LEAK['SUM']-LEAK['DIFF'])*(1-T0)*raw),atol=1e-12)
        table=readcsv(paths[2]);a.check('Extra CSV has every displacement',len(table)==len(s))
        if len(table)==len(s):
            for key in z.files:
                a.close('Extra CSV matches saved '+key,np.array([float(r[key]) for r in table]),z[key],atol=1e-12)
    a.stats['extra_independent_holdout']=ID
    # Requested weak-signal time controls must complete before final delivery.
    controlpath=folder/'control_curves.npz'
    for suffix in ('_fine','_fine_t400'):
        inspect_case(a,ID+suffix,v);inspect_spectrum(a,ID,suffix=suffix)
    cancelled=ROOT/'fdtd/cases'/(ID+'_finerY')
    if (cancelled/'CANCELLED.json').exists():
        a.check('Cancelled weak finerY run is not claimed completed',not (cancelled/'DONE.json').exists() and all(r.get('suffix')!='_finerY' for r in summary.get('time_mesh_control_rows',[])))
    if not a.need(controlpath,'Weak100/200/400fs saved control curves'):return
    with np.load(controlpath) as z,np.load(target) as t:
        need={suffix+'_'+ch for suffix in ('_fine','_fine_t200','_fine_t400') for ch in LEAK}
        if not need<=set(z.files):
            detail=dict(missing=sorted(need-set(z.files)))
            if a.allow_partial:a.pending('Weak400fs time comparison not yet prepared',detail)
            else:a.check('Weak100/200/400fs control vectors complete',False,detail)
            return
        a.check('Weak100/200/400fs six control vectors complete',summary.get('controls_complete') is True)
        diffs={}
        for ch,key in (('SUM','S_SUM'),('DIFF','S_DIFF_same_kernel')):
            span=float(np.ptp(t[key]));value=100*float(np.max(abs(z['_fine_t200_'+ch]-z['_fine_t400_'+ch])))/span
            diffs[ch]=value
            a.close('Weak '+ch+' reported time error uses own weak span',summary['weak_200_vs_400_max_difference_over_own_span_pct'][ch],value,atol=1e-10)
        a.check('Weak200fs accepted only after <0.2% own-span 400fs control',summary.get('time_control_passed') is True and max(diffs.values())<.2,diffs)
        a.stats['weak_200_vs_400_own_span_pct']=diffs


def inspect_solver_manifest(a):
    """Audit every completed solve, including coarse/time/mesh controls.

    Generic checks use each case's saved mesh/PML rather than imposing the
    primary fine-grid settings on deliberately different numerical controls.
    """
    base=ROOT/'fdtd';paths=[base/'solver_manifest.json',base/'solver_manifest.csv',base/'saved_project_audit.json']
    if not all([a.need(p,'Global solver provenance '+p.name) for p in paths]):return
    manifest=readjson(paths[0]);table=readcsv(paths[1]);native=readjson(paths[2])
    rows={x['case']:x for x in manifest['cases']};csvrows={x['case']:x for x in table};saved={x['case']:x for x in native['cases']}
    actual={p.parent.name for p in (base/'cases').glob('*/DONE.json')}
    a.check('All 21 genuine DONE cases enumerated exactly once',len(actual)==21 and len(rows)==21 and len(table)==21 and actual==set(rows)==set(csvrows)==set(saved))
    a.check('Solver manifest completed and cancelled counts match disk',manifest.get('status')=='passed' and manifest.get('actual_completed_case_count')==len(actual) and manifest.get('cancelled_case_count')==1)
    a.check('Saved native project readback covers all21 without new solves',native.get('status')=='passed' and native.get('read_only') is True and native.get('new_solver_calls')==0 and native.get('case_count')==len(actual))
    for name in sorted(actual):
        folder=base/'cases'/name;d=readjson(folder/'DONE.json');s=d['settings'];term=d['termination'];row=rows[name];cs=csvrows[name];ns=saved[name]
        obs=readjson(folder/'observed_settings.json');engine=readjson(folder/'engine_result.json')
        a.check('All21 '+name+' real successful engine',d.get('status')=='completed' and term.get('completed_successfully') is True and engine.get('returncode',engine.get('return_code'))==0 and row.get('engine_returncode')==0 and row.get('engine_success') is True)
        a.check('All21 '+name+' model hash agrees with manifest',a.hash(folder/'model.fsp')==d.get('model_sha256')==row.get('model_sha256')==cs.get('model_sha256'))
        a.check('All21 '+name+' actual saved four PML boundaries',all(obs.get(k)==ns.get(k)=='PML' for k in ('x min bc','x max bc','y min bc','y max bc')) and row.get('actual_all_PML') is True and cs.get('actual_all_PML')=='True')
        a.close('All21 '+name+' PML requested/observed/native/CSV',
            [s['pml_layers'],obs['pml layers'],ns['pml layers'],float(cs['pml_actual_max']),row['pml_actual_max']],np.full(5,row['pml_requested']),atol=1e-12)
        a.close('All21 '+name+' actual recorded duration',
            [row['actual_time_fs'],float(cs['actual_time_fs'])],np.full(2,term['simulated_time_s']*1e15),atol=1e-7)
        a.close('All21 '+name+' actual numerical mesh',
            [obs['dx'],obs['dy'],ns['dx'],ns['dy']], [s['dx_m'],s['dy_m'],s['dx_m'],s['dy_m']],atol=1e-15)
        for mon in ('R','R2','INC','T'):
            path=folder/f'nearfield_{mon}.npz'
            if not a.need(path,'All21 '+name+' raw '+mon):continue
            with np.load(path) as z:
                x=np.asarray(z['x_m']).squeeze()
                a.check('All21 '+name+' '+mon+' six finite complex fields',len(x)>100 and np.all(np.diff(x)>0) and all(k in z.files and np.iscomplexobj(z[k]) and z[k].size==x.size and np.isfinite(z[k]).all() for k in COMPS) and row.get('all_raw_complex_fields_valid') is True)
                a.close('All21 '+name+' '+mon+' true target wavelength',z['wavelength_m'],532e-9,atol=1e-15)
    cancelled=manifest.get('cancelled_cases',[])
    a.check('Cancelled weak finerY transparently excluded from21',len(cancelled)==1 and cancelled[0].get('case')=='ah53_h180_cap400_finerY' and cancelled[0].get('completed') is False and cancelled[0].get('DONE_written') is False and cancelled[0]['case'] not in actual and (base/'cases'/cancelled[0]['case']/'CANCELLED.json').exists())
    a.stats['completed_solver_case_count']=len(actual);a.stats['cancelled_solver_case_count']=len(cancelled)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--allow-partial',action='store_true');ap.add_argument('--report');args=ap.parse_args()
    a=Audit(args.allow_partial);plan=readjson(ROOT/'experiment.json');variants=plan['variants']
    a.check('Exactly eight intended structures',len(variants)==8 and len({v['id'] for v in variants})==8)
    a.check('User confirmed Tm and K2/K4 in current contract','confirmed' in plan['working_formula_contract']['T_m'] and 'authorized' in plan['working_formula_contract']['K2_K4'].lower())
    a.group('Fine empty-domain reference',lambda:inspect_case(a,'air_fine'))
    complete=[]
    for v in variants:
        a.group(v['id']+' FDTD provenance',lambda v=v:inspect_case(a,v['id']+'_fine',v))
        if (ROOT/'fdtd/cases'/(v['id']+'_fine')/'DONE.json').exists():complete.append(v['id'])
        a.group(v['id']+' processed fields',lambda v=v:inspect_spectrum(a,v['id']))
    a.stats['fine_DONE_cases_found']=complete
    a.group('RCWA reference',lambda:inspect_rcwa(a,variants))
    a.group('Final comparison',lambda:inspect_comparison(a,variants))
    a.group('Physical-component main comparison',lambda:inspect_physical_component(a,variants))
    a.group('Extra weak heldout',lambda:inspect_extra(a,variants))
    a.group('Every completed solver case and cancellation',lambda:inspect_solver_manifest(a))
    a.group('Final report links',lambda:inspect_reports(a,args.report))
    counts={k:sum(x['status']==k for x in a.items) for k in ('PASS','FAIL','PENDING')}
    result=dict(status='complete' if not counts['FAIL'] and not counts['PENDING'] else 'failed' if counts['FAIL'] else 'partial',
        passed=not counts['FAIL'] and not counts['PENDING'],allow_partial=args.allow_partial,
        inspected_utc=datetime.now(timezone.utc).isoformat(),counts=counts,stats=a.stats,checks=a.items,
        scope='Delivery completeness, real solve provenance and scientific cross-file consistency. No missing simulation is counted as passed; this audit is not absolute Maxwell accuracy certification.')
    out=ROOT/'qa'/('delivery_audit_partial.json' if args.allow_partial else 'delivery_audit.json')
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:result[k] for k in ('status','passed','counts','stats')},ensure_ascii=False,indent=2))
    for x in a.items:
        if x['status']=='FAIL':print('FAIL:',x['check'],x['detail'])
    if counts['FAIL'] or (counts['PENDING'] and not args.allow_partial):return 1
    return 0


if __name__=='__main__':sys.exit(main())
