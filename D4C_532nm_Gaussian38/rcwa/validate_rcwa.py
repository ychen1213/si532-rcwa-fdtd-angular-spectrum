"""Meaningful checks for the extra cap layer, phase and flux conventions."""
from dataclasses import replace
from pathlib import Path
import importlib.util
import numpy as np
from te_planarized_rcwa import Config,solve_angle,outgoing_sqrt,dump

ROOT=Path(__file__).resolve().parent


def film_r(q,c,thickness):
    ya=outgoing_sqrt(c.air_n**2-(q/c.k0)**2).item()
    yf=outgoing_sqrt(c.fill_n**2-(q/c.k0)**2).item()
    ys=outgoing_sqrt(c.substrate_n**2-(q/c.k0)**2).item()
    raf=(ya-yf)/(ya+yf);rfs=(yf-ys)/(yf+ys)
    e=np.exp(2j*c.k0*yf*thickness)
    # r at the top boundary, then both incident and reflected propagation to y=0.
    return (raf+rfs*e)/(1+raf*rfs*e)*np.exp(-2j*c.k0*ya*thickness)


def main():
    checks=[]
    def ck(name,value,limit,detail=''):
        checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<limit),detail=detail))
    base=Config(ridge_n=1.46+0j,fill_n=1.46+0j,cap_n=1.46+0j,substrate_n=4.15+.05j)
    for cap in (0.,100e-9,200e-9):
        c=replace(base,cap_m=cap)
        for frac in (0.,.04,.3):
            q=frac*c.k0;r=solve_angle(q,16,c);j=16
            ck(f'Uniform_film_Fresnel_cap{cap*1e9:g}_q{frac}',abs(r['r_total_y0'][j]-film_r(q,c,c.top_m)),2e-12)
            ck(f'Uniform_lossless_film_interface_flux_cap{cap*1e9:g}_q{frac}',abs(r['finite_layer_absorption']),2e-12)
            other=np.delete(r['r_total_y0'],j)
            ck(f'Uniform_film_no_spurious_orders_cap{cap*1e9:g}_q{frac}',np.max(abs(other)),2e-12)
    # Check the new code against the previously native-validated single-layer implementation.
    old_path=ROOT/'provenance'/'original_si_rcwa.py'
    spec=importlib.util.spec_from_file_location('archived_si_rcwa',old_path)
    old=importlib.util.module_from_spec(spec)
    import sys;sys.modules[spec.name]=old;spec.loader.exec_module(old)
    for q in (0.,125000.):
        c=Config(cap_m=0.)
        oc=old.Config(wavelength_m=c.wavelength_m,period_m=c.period_m,height_m=c.height_m,
            ridge_width_m=c.ridge_width_m,ridge_centers_m=c.ridge_centers_m,
            air_n=c.air_n,ridge_n=c.ridge_n,substrate_n=c.substrate_n)
        a=solve_angle(q,80,c);b=old.solve_angle(q,80,config=oc)
        ck(f'Archived_single_layer_regression_q{q}',np.max(abs(a['r_total_y0']-b['r_total_y0'])),1e-10)
    c=Config(cap_m=150e-9,fill_n=1.46+0j,cap_n=1.46+0j)
    p=solve_angle(125000.,96,c);m=solve_angle(-125000.,96,c)
    ck('Capped_structure_mirror_symmetry',np.max(abs(p['r_total_y0']-m['r_total_y0'][::-1])),2e-10)
    ck('Capped_structure_passivity',max(-p['finite_layer_absorption'],0.),1e-11)
    ck('Capped_Si_absorption_positive',max(1e-6-p['finite_layer_absorption'],0.),1e-11)
    result=dict(passed=all(v['passed'] for v in checks),checks=checks,
        scope='TE equations; uniform capped-film Fresnel, flux, original solver regression, symmetry, passivity. No FDTD data fitted.')
    dump(ROOT/'solver_validation.json',result)
    print(f"{sum(v['passed'] for v in checks)}/{len(checks)} RCWA checks passed")
    if not result['passed']:raise RuntimeError('RCWA validation failed')


if __name__=='__main__':main()
