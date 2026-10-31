"""Obtain database and fitted indices for matching RCWA/FDTD at 532 nm."""
from pathlib import Path
import json
import numpy as np
from run_rcwa import session
root=Path(__file__).resolve().parent
f=299792458/532e-9
with session() as s:
    material='Si (Silicon) - Palik'
    n=complex(np.asarray(s.getindex(material,f)).ravel()[0])
    nf=complex(np.asarray(s.getfdtdindex(material,f,299792458/565e-9,299792458/500e-9)).ravel()[0])
    result={'lambda_nm':532,'material':material,'getindex':{'real':n.real,'imag':n.imag},'getfdtdindex':{'real':nf.real,'imag':nf.imag},'fit_band_nm':[500,565]}
    print(json.dumps(result),flush=True)
    (root/'material.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    fixed=s.addmaterial('(n,k) Material')
    s.setmaterial(fixed,'name','Si532_fitted_fixed')
    s.setmaterial('Si532_fitted_fixed','Refractive Index',nf.real)
    s.setmaterial('Si532_fitted_fixed','Imaginary Refractive Index',nf.imag)
    checked=complex(np.asarray(s.getindex('Si532_fitted_fixed',f)).ravel()[0])
    assert abs(checked-nf)<1e-14
    print('Fixed (n,k) material verified:',checked,flush=True)
