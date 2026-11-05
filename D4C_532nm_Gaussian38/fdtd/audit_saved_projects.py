"""Read actual completed native projects, including BC, mesh and geometry.

This opens a fresh hidden CAD session after solves; no simulation is run.
"""
from pathlib import Path
import json
import numpy as np
from runtime_support import hidden_fdtd
from run_fdtd import ROOT,NATIVE,dump,clean,C0,SI,OX

def main():
    rows=[]
    with hidden_fdtd(8) as s:
        for path in sorted((ROOT/'cases').glob('*/DONE.json')):
            folder=path.parent;c=json.loads((folder/'settings.json').read_text(encoding='utf-8'));s.load(str(NATIVE/'cases'/folder.name/'model.fsp'))
            actual={k:s.getnamed('FDTD',k) for k in ('dimension','x min bc','x max bc','y min bc','y max bc','pml layers','mesh type','dx','dy','dt','simulation time','auto shutoff min')}
            actual['PML_parameters']={}
            for prop in ('pml profile','pml kappa','pml sigma','pml alpha','pml polynomial','pml alpha polynomial','pml min layers','pml max layers'):
                try:actual['PML_parameters'][prop]=s.getnamed('FDTD',prop)
                except Exception:pass
            actual['case']=folder.name;actual['source']={k:s.getnamed('source',k) for k in ('x','y','x span','waist radius w0','distance from waist','amplitude','phase','polarization angle')}
            actual['numerical_permittivity']={key:complex(np.asarray(s.getnumericalpermittivity(material,C0/532e-9,C0/600e-9,C0/500e-9,float(np.asarray(actual['dt']).squeeze()))).squeeze()) for key,material in [('Si',SI),('SiO2',OX)]}
            actual['geometry']={}
            actual['oxide_object_count']=int(np.asarray(s.getnamednumber('SiO2_planarization')).squeeze())
            names=[]
            if c['case']!='air':names.append('Si_substrate')
            if c.get('geometry',{}).get('cap_nm',0)>0:names.append('SiO2_planarization')
            v=c.get('geometry',{});centers=v.get('ridge_centers_local_um',[])
            names += [f'Si_ridge_{i:02d}_{j}' for i in range(v.get('periods',0)) for j in range(len(centers))]
            for name in names:actual['geometry'][name]={k:s.getnamed(name,k) for k in ('material','x','x span','y min','y max','mesh order')}
            actual['passed_all_PML']=all(actual[k]=='PML' for k in ('x min bc','x max bc','y min bc','y max bc'))
            actual['requested_pml_layers']=c['pml_layers'];actual['actual_pml_max']=int(np.max(np.asarray(actual['pml layers'])))
            actual['cap000_bare_si']=(v.get('cap_nm')!=0 or actual['oxide_object_count']==0)
            actual['cap_positive_top_correct']=(v.get('cap_nm',0)==0 or abs(float(actual['geometry']['SiO2_planarization']['y max'])-(v['ridge_height_nm']+v['cap_nm'])*1e-9)<1e-15)
            actual['passed']=actual['passed_all_PML'] and actual['cap000_bare_si'] and actual['cap_positive_top_correct']
            rows.append(clean(actual));print('AUDITED',folder.name,actual['passed'],flush=True)
    out=dict(status='passed' if all(r['passed'] for r in rows) else 'failed',read_only=True,new_solver_calls=0,case_count=len(rows),cases=rows)
    dump(ROOT/'saved_project_audit.json',out)

if __name__=='__main__':main()
