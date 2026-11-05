"""Independent angular-norm and zero-accounting audit of saved fine FDTD.

Read-only inputs; no electromagnetic solve, target fitting, or parent-driver edit.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from zero_order_models import (T0,LEAKAGE,KNOWN_ORDERS,angular_weights,b_features,
    candidate_scan,instrument_weights,nearest_order,scan_from_features,stop_cell_field)

ROOT=Path(__file__).resolve().parents[1]


def audit(case='ah53_cap000_fine'):
    p=ROOT/'outputs/spectra'/f'{case}.npz'
    with np.load(p) as z:d={name:z[name] for name in z.files}
    k=d['kx_rad_m'];a=d['A_ref_V'];a0=d['A_zero_degree2_V'];pin=float(d['Pin_W_m'])
    half=.17991979895042
    f=b_features(k,a,pin,selected_orders=KNOWN_ORDERS,stop_halfwidth_rad_um=half,A_zero_ref_V=a0)
    s=np.array([-4.,-3.1,-1.7,0.,.335,2.17,4.])
    compressed=scan_from_features(f,s,convention='K_over_gamma')
    # Independent full-grid norm: one coherent field at every s, then integration.
    w=angular_weights(k,pin);m=nearest_order(k);G=2*np.pi/16e-6
    mask=stop_cell_field(k,half);field=(a-a0)*mask
    iw=np.zeros_like(k)
    for mm,weight in zip(KNOWN_ORDERS,instrument_weights(KNOWN_ORDERS,convention='K_over_gamma')):
        iw[abs(m)==mm]=weight
    direct=[]
    for ss in s:
        shifted=field*np.exp(-1j*m*G*ss*1e-6)
        direct.append(float(.5*np.sum(w*iw*abs(shifted-shifted[::-1])**2)))
    error=float(max(abs(np.asarray(direct)-compressed)))
    zero_rows=[]
    for channel in ('SUM','DIFF'):
        scalar=candidate_scan(f,s,model='scalar_zero_once',convention='K_over_gamma',channel=channel)
        physical=candidate_scan(f,s,model='physical_zero_once',convention='K_over_gamma',channel=channel)
        duplicate=candidate_scan(f,s,model='double_counted_zero',convention='K_over_gamma',channel=channel)
        difference=(1-LEAKAGE[channel])*(f.zero_pass_power-T0*f.zero_raw_power)
        zero_rows.append(dict(channel=channel,
            physical_minus_scalar_expected_constant=difference,
            physical_minus_scalar_identity_error=float(max(abs((physical-scalar)-difference))),
            duplicate_minus_scalar_expected_constant=f.zero_pass_power,
            duplicate_minus_scalar_identity_error=float(max(abs((duplicate-scalar)-f.zero_pass_power))),
            physical_vs_scalar_scan_shape_range=float(np.ptp(physical-scalar))))
    closure=f.zero_pass_power+f.nonzero_pass_power+f.zero_nonzero_cross_pass_power-f.collected_power
    assert error<1e-12 and abs(closure)<1e-12
    assert all(r['physical_minus_scalar_identity_error']<1e-12 and r['duplicate_minus_scalar_identity_error']<1e-12 for r in zero_rows)
    sources={name:hashlib.sha256((ROOT/'analysis'/name).read_bytes()).hexdigest()
        for name in ('prepare_spectra.py','run_analysis.py','zero_order_models.py')}
    report=dict(case=case,saved_spectrum_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source_sha256=sources,
        direct_full_angular_norm_vs_compressed_features_max_error=error,
        modal_power_closure_error=closure,zero_accounts=zero_rows,
        current_symmetric_mark_blocking_math_bugs=[],
        confirmed=[
            'Reflection propagation exp(+i beta*(yref-ymonitor)) and downward incident propagation use consistent opposite signs',
            'Actual Pin enters angular_weights exactly once; gamma(k) and separate K_m/gamma_m are not conflated',
            'physical_zero_once driver records nonzero+Ppass and blocked zero, with fixed tau=0 giving exactly L*Pblocked',
            'A0 extraction is fit to the saved complex field and measured source only, not to a D4C scan target',
            'Power decomposition retains the A0/nonzero cross term as an audit observable',
            'Every principal zero-account candidate changes only a constant if its nonzero fields/stop are otherwise the same; a better one-case RMSE does not establish correct accounting'],
        limitations=[
            'Source projection is a finite polynomial model; total-field L2 can conceal greater weak-order error',
            'Removing A0 and then adding its scalar power deliberately omits A0/nonzero coherent cross terms according to the D4C signal contract',
            'B dark-pair phase gives cos(2mGs-arg(Cm)); target note uses +phi. All eight marks are symmetric so phi is essentially zero; asymmetric marks need the phase convention reconciled',
            'Two leakage-labelled curves do not prove actual SUM/DIFF hardware complementarity or a channel interference sign'],
        fitting_constraint='One parameter vector per declared model across training marks, frozen on holdout; no per-mark offset permitted')
    driver=(ROOT/'analysis/run_analysis.py').read_text(encoding='utf8')
    if 'def periodic_records(' in driver:
        report['confirmed'].append('Current root driver excludes the duplicated periodic endpoint from shared fitting as well as final metrics')
    else:
        report['limitations'].append('Shared fitting includes both periodic endpoints while metrics exclude the duplicate; this is a small weighting difference')
    out=ROOT/'analysis'/'root_implementation_audit.json'
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    return report


if __name__=='__main__':audit()
