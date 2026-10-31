"""Publish numerical solver evidence without local machine/licensing details.

Original solver logs are preserved unchanged beside the public derivatives.
"""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parent

KEEP=re.compile(r'^\s*(?:Ansys Lumerical|CPU threads per process|Running [23]D Simulation|Simulation size in gridpoints|Starting meshing|Completed meshing|Meshing complete|Refined \d|Finished \d+ iterations|Total [XYZ] spatial steps|FDTD loop runs at|[\d.eE+-]+ seconds total loop time|[\d.]+% complete|Early termination of simulation|Completed \d+ iterations|Simulation completed successfully|Peak CPU memory|Overall wall time|time to mesh|time to run|time to checkpoint|time to finalize|miscellaneous time|FDTD solver speed|total FDTD solver speed)')
PRIVATE=re.compile(r'license|hostname|Running on host|license host|config file|Running fsp file|[A-Za-z]:[\\/]|/Users/|/home/',re.I)

def sanitize():
    outputs=[]
    dest=ROOT/'public_log_manifest.json'
    previous=json.loads(dest.read_text(encoding='utf-8')) if dest.exists() else []
    for src in list(ROOT.glob('cases/*/*/engine_stdout.log'))+list(ROOT.glob('cases/*/*/model_p0.log')):
        dst=src.with_name(src.stem+'_public.log')
        lines=[line for line in src.read_text(encoding='utf-8',errors='replace').splitlines() if KEEP.match(line) and not PRIVATE.search(line)]
        dst.write_text('# Public numerical excerpt; original local solver log is preserved unchanged.\n'+'\n'.join(lines)+'\n',encoding='utf-8')
        outputs.append({'source':str(src.relative_to(ROOT)),'public':str(dst.relative_to(ROOT)),'retained_lines':len(lines)})
    merged={r['public']:r for r in previous};merged.update({r['public']:r for r in outputs})
    outputs=list(merged.values());dest.write_text(json.dumps(outputs,indent=2),encoding='utf-8')
    return outputs

def audit_saved_projects(s,cases=None):
    from run_fdtd import ASCII_ROOT
    cases=cases or ['N010_dx20_dy10_t600','N400_dx20_dy10_t600','N020_dx10_dy5_t600']
    rows=[]
    for tag in cases:
        p=ROOT/'cases'/tag/'mark'/'model.fsp'
        if not p.exists():continue
        s.load(str(ASCII_ROOT/p.relative_to(ROOT)));record={'case':tag,'read_method':'saved solved model.fsp; read-only property queries; no solve'}
        for obj,properties in [('FDTD',['dimension','x min bc','x max bc','y min bc','y max bc','x span','y min','y max','dx','dy','pml layers','pml profile','simulation time','auto shutoff min','dt stability factor']),('source',['plane wave type','injection axis','direction','x','x span','y','z','polarization angle','angle theta','wavelength start','wavelength stop'])]:
            record[obj]={}
            for prop in properties:
                try:
                    v=s.getnamed(obj,prop)
                    record[obj][prop]=v.tolist() if hasattr(v,'tolist') else v
                except Exception as ex:record[obj][prop]={'read_error':str(ex)}
        record['solver_version']=s.version()
        target=p.parent/'saved_project_audit.json';target.write_text(json.dumps(record,indent=2),encoding='utf-8');rows.append(record)
    merged={}
    dest=ROOT/'saved_project_audit.json'
    if dest.exists():merged={r['case']:r for r in json.loads(dest.read_text(encoding='utf-8'))}
    merged.update({r['case']:r for r in rows});dest.write_text(json.dumps(list(merged.values()),indent=2),encoding='utf-8')
    return rows

if __name__=='__main__':
    out=sanitize();print('Public numerical log excerpts:',len(out))
