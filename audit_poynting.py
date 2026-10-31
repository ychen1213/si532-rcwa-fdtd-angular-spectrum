"""Independent finite-monitor flux audit; never runs a Maxwell solver.

Compare native +y physical Poynting flux of saved scattered Ez/Hx fields
with the outgoing-wave E-only angular-spectrum estimate used by analyze.py.
Records the actual incident Poynting normalization as a second denominator.
"""
from pathlib import Path
import argparse
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
Z0=376.730313668
K0=2*np.pi/.532  # rad / micrometre

def load_line(path,L,dx):
    with np.load(path) as z:
        xx=z['x_m'].ravel()*1e6
        x=np.linspace(-L/2,L/2,int(round(L/dx))+1)
        assert xx[0]<=x[0] and xx[-1]>=x[-1]
        fields=[]
        for key in ['Ez','Hx']:
            a=z[key].ravel()
            fields.append(np.interp(x,xx,a.real)+1j*np.interp(x,xx,a.imag))
    return x,*fields

def audit_case(case):
    cfg=json.loads((case/'mark/DONE.json').read_text(encoding='utf-8'))['settings']
    L=cfg['N']*cfg['pitch_m']*1e6
    dx=cfg['dx_m']*1e6
    x,e,h=load_line(case/'mark/fields_R.npz',L,dx)
    _,ea,ha=load_line(case/'air/fields_R.npz',L,dx)
    _,g,gh=load_line(case/'air/fields_I.npz',L,dx)
    e=e-ea;h=h-ha
    g0=np.trapezoid(g,x)/L
    p_normalization=.5/Z0*L*1e-6*abs(g0)**2
    p_inc_actual=-.5e-6*np.trapezoid(np.real(g*gh.conj()),x)
    p_eh=.5e-6*np.trapezoid(np.real(e*h.conj()),x)
    weights=np.ones(len(x));weights[[0,-1]]=.5
    nfft=2**int(np.ceil(np.log2(len(x)*16)))
    q=2*np.pi*np.fft.fftfreq(nfft,dx)
    # A common x-origin factor cancels in AE*conj(AH), and in |AE|^2.
    AE=dx*np.fft.fft(e*weights,nfft)
    AH=dx*np.fft.fft(h*weights,nfft)
    dq=2*np.pi/(nfft*dx)
    p_eh_spectrum=.5e-6*np.real(np.sum(AE*AH.conj()))*dq/(2*np.pi)
    p_discrete_squaredweights=.5e-6*dx*np.sum(weights**2*np.real(e*h.conj()))
    keep=abs(q)<K0
    beta_over_k=np.sqrt(1-(q[keep]/K0)**2)
    p_eonly=1e-6/(4*np.pi*Z0)*np.sum(beta_over_k*abs(AE[keep])**2)*dq
    # Same trapezoidal q quadrature as the primary analyzer, for exact comparison.
    idx=np.argsort(q[keep]);qs=q[keep][idx]
    p_eonly_trap=1e-6/(4*np.pi*Z0)*np.trapezoid(
        beta_over_k[idx]*abs(AE[keep][idx])**2,qs)
    return dict(case=case.name,N=cfg['N'],L_um=L,
        cropped_monitor_samples=len(x),
        reflected_poynting_W_per_m_slot=p_eh,
        reflected_Eonly_spectrum_W_per_m_slot=p_eonly_trap,
        reflected_EH_spectrum_W_per_m_slot=p_eh_spectrum,
        incident_uniform_mean_normalization_W_per_m_slot=p_normalization,
        actual_incident_poynting_W_per_m_slot=p_inc_actual,
        R_Poynting_uniform_mean_normalization=p_eh/p_normalization,
        R_Eonly_spectrum_uniform_mean_normalization=p_eonly_trap/p_normalization,
        R_Poynting_actual_incident_normalization=p_eh/p_inc_actual,
        R_Eonly_spectrum_actual_incident_normalization=p_eonly_trap/p_inc_actual,
        Eonly_spectrum_minus_poynting_relative=(p_eonly_trap-p_eh)/p_eh,
        actual_incident_over_uniform_mean=p_inc_actual/p_normalization,
        EH_Parseval_relative_error=(p_eh_spectrum-p_discrete_squaredweights)/p_discrete_squaredweights,
        FFT_endpoint_quadrature_relative_difference=(p_eh_spectrum-p_eh)/p_eh,
        Eonly_q_sum_minus_trapezoid_relative=(p_eonly-p_eonly_trap)/p_eh)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--N',type=int,nargs='+',default=[100,400])
    args=parser.parse_args()
    rows=[]
    for N in args.N:
        case=ROOT/'fdtd/cases'/f'N{N:03d}_dx20_dy10_t600'
        if not all((case/s/'DONE.json').exists() for s in ['air','mark']):
            print(f'PENDING: {case.name}',flush=True)
            continue
        result=audit_case(case);rows.append(result)
        print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    out=ROOT/'analysis/poynting_audit.json'
    previous=json.loads(out.read_text(encoding='utf-8')) if out.exists() else []
    table={r['case']:r for r in previous}
    table.update({r['case']:r for r in rows})
    out.write_text(json.dumps(list(table.values()),ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__':main()
