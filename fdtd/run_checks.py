"""Small, explicit numerical controls after the finite-length production series."""
from datetime import datetime,timezone
import traceback
from run_fdtd import hidden_fdtd,params,run_one,jswrite,ROOT
from publish_diagnostics import audit_saved_projects,sanitize
from run_strict_time import strict_time

def main():
    jobs=[
      (params(200,t_fs=1200),'N200_dx20_dy10_t1200'),
      (params(50,dx_nm=10,dy_nm=5,pml=48),'N050_dx10_dy5_t600'),
      (params(20,dx_nm=10,dy_nm=5,pml=48),'N020_dx10_dy5_t600_pml48'),
      (params(20,dx_nm=5,dy_nm=2.5,pml=96),'N020_dx5_dy2.5_t600')]
    with hidden_fdtd(threads=8) as s:
        audit_saved_projects(s)
        for cfg,tag in jobs:
            for air in [True,False]:
                try:run_one(s,cfg,air,tag)
                except Exception as ex:
                    directory=ROOT/'cases'/tag/('air' if air else 'mark');directory.mkdir(parents=True,exist_ok=True)
                    jswrite(directory/'FAILED.json',dict(status='failed',settings=cfg,error=str(ex),traceback=traceback.format_exc(),utc=datetime.now(timezone.utc).isoformat()))
                    raise
            sanitize()
        strict_time(s)
        audit_saved_projects(s,['N050_dx10_dy5_t600','N020_dx5_dy2.5_t600','N200_dx20_dy10_t1200_strictmark'])
        sanitize()

if __name__=='__main__':main()
