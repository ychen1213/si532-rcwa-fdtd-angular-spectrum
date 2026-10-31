"""Archive actual RCWA sampling and inspect the sampled, centered ridge."""
from pathlib import Path
import json
import numpy as np
from run_rcwa import ROOT,session,write_json
with session() as s:
    for out in sorted(ROOT.glob('harmonics*')):
        s.load(str(out/'model.fsp'))
        metadata={key:s.getnamed('RCWA',key) for key in ['actual mesh cells x','interface actual positions']}
        s.eval('ii=getresult("RCWA","index"); ix=ii.x; iz=ii.z; iattr=getattribute(ii);')
        metadata['index_attributes']=s.getv('iattr')
        print(out.name,metadata,flush=True)
        x=np.asarray(s.getv('ix')).ravel();z=np.asarray(s.getv('iz')).ravel()
        s.eval('inz=ii.index_z;')
        index_z=np.asarray(s.getv('inz'))
        np.savez_compressed(out/'sampled_index.npz',x=x,z=z,index_z=index_z)
        write_json(out/'actual_mesh.json',metadata)
        print('x',x.shape,x[:4]*1e9,x[-4:]*1e9,'z',z*1e9,'index',index_z.shape,flush=True)
