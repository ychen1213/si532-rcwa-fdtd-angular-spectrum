"""Force the longer mark integration, reusing the genuine full-duration air run."""
import json
from run_fdtd import ROOT,params,hidden_fdtd,run_one,jswrite

def strict_time(s):
    reference='N200_dx20_dy10_t1200/air'
    original=json.loads((ROOT/'cases'/reference/'DONE.json').read_text(encoding='utf-8'))
    cfg=params(200,t_fs=1200,shutoff=1e-12)
    for key in cfg:
        if key!='auto_shutoff' and cfg[key]!=original['settings'][key]:raise ValueError('Air-reference physics differs: '+key)
    actual=original['solver_statistics']['simulated_time_s']
    if actual+1e-17<cfg['t_max_s'] or original['solver_statistics']['terminated_early']:
        raise ValueError('The reference must really finish the full requested integration time.')
    directory=ROOT/'cases/N200_dx20_dy10_t1200_strictmark';directory.mkdir(exist_ok=True)
    jswrite(directory/'AIR_REFERENCE.json',dict(air_reference=reference,reference_reused=True,path_base='fdtd/cases',
      reference_auto_shutoff_min=original['settings']['auto_shutoff'],reference_actual_time_s=actual,
      reference_requested_time_s=original['settings']['t_max_s'],reference_completed_all_requested_time=True,
      new_mark_auto_shutoff_min=cfg['auto_shutoff'],
      reason='The existing air job completed all 1200 fs without satisfying its 1e-7 stopping threshold. Lowering that unused threshold changes no Maxwell update or DFT time step; this full-duration reference is reused explicitly. No new air solve or DONE is claimed.'))
    return run_one(s,cfg,False,directory.name)

if __name__=='__main__':
    from publish_diagnostics import audit_saved_projects,sanitize
    with hidden_fdtd(threads=8) as s:
        strict_time(s)
        audit_saved_projects(s,['N200_dx20_dy10_t1200_strictmark'])
        sanitize()
