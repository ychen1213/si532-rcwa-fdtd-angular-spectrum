"""Small independent diagnostics; full report analysis lives at project root."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
k0=2*np.pi/532e-9

def air_directory(case):
    link=case/'AIR_REFERENCE.json'
    if link.exists():return ROOT/'cases'/json.loads(link.read_text(encoding='utf-8'))['air_reference']
    return case/'air'

def extract(case):
    ad=air_directory(case)
    a=np.load(ad/'fields_R.npz');b=np.load(case/'mark/fields_R.npz');i=np.load(ad/'fields_I.npz')
    cfg=json.loads((case/'mark/input.json').read_text(encoding='utf-8'))['settings']
    x=b['x_m'].ravel(); L=cfg['N']*cfg['pitch_m'];q=(x>=-L/2-1e-12)&(x<L/2-1e-12)
    e=(b['Ez']-a['Ez']).ravel();inc=i['Ez'].ravel()[q].mean()
    r=np.array([(e[q]*np.exp(-1j*2*np.pi*m*x[q]/2e-6)).mean()/inc*np.exp(-1j*np.sqrt(k0*k0-(2*np.pi*m/2e-6)**2)*cfg['R_y_m']) for m in range(-3,4)])
    return dict(case=case.name,r_re=r.real.tolist(),r_im=r.imag.tolist(),incident_re=inc.real,incident_im=inc.imag,
       air_crop_nonuniformity=float(np.max(abs(i['Ez'].ravel()[q]-inc))/abs(inc)),
       samples=int(q.sum()),mark_solve_s=json.loads((case/'mark/DONE.json').read_text(encoding='utf-8'))['solve_s'])

if __name__=='__main__':
    rows=[]
    for c in sorted((ROOT/'cases').iterdir()):
        if (c/'mark/DONE.json').exists() and (air_directory(c)/'DONE.json').exists():rows.append(extract(c))
    print(json.dumps(rows,indent=2));(ROOT/'quick_summary.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
