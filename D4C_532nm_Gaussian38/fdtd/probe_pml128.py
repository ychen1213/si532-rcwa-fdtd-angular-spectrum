"""Probe configurable PML max layers in a fresh unsaved project, no solve."""
from pathlib import Path
import json
import numpy as np
from runtime_support import hidden_fdtd
from run_fdtd import ROOT,dump,clean

def main():
    trials=[];success=None
    with hidden_fdtd(8) as s:
        s.addfdtd();s.set('dimension','2D')
        for profile in (1,4):
            s.set('pml profile',profile)
            for prop in ('pml max layers','pml maximum layers','maximum pml layers'):
                row=dict(profile=profile,property=prop)
                try:
                    row['before']=s.getnamed('FDTD',prop)
                    s.setnamed('FDTD',prop,128)
                    s.setnamed('FDTD','pml layers',128)
                    row['max_after']=s.getnamed('FDTD',prop);row['layers_actual']=s.getnamed('FDTD','pml layers')
                    row['passed']=int(np.max(np.asarray(row['layers_actual'])))==128
                    if row['passed']:success=dict(profile=profile,max_layer_property=prop,max_layers=128)
                except Exception as e:row['error']=str(e)
                trials.append(clean(row))
                if success:break
            if success:break
    result=dict(no_simulation=True,no_saved_project_modified=True,success=success,trials=trials)
    dump(ROOT/'pml_128_probe.json',result);print(json.dumps(result,indent=2))
    if not success:raise RuntimeError('Could not establish actual128-layer PML')

if __name__=='__main__':main()
