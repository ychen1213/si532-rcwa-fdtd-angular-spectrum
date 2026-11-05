"""Fit measured empty-domain Gaussian spectra and record incident power."""
import json
import numpy as np
from audit_fields import ROOT,Z0,angular,dump

def main():
    out={};analytic=float(np.sqrt(np.pi/2)*19e-6/(2*Z0))
    for p in sorted((ROOT/'cases').glob('air*/DONE.json')):
        case=p.parent.name;I=angular(case,'INC');k=I['k'];a=I['down'];sel=abs(a)>max(abs(a))*.02
        coeff=np.polyfit((k[sel]*1e-6)**2,np.log(abs(a[sel])),1)
        out[case]=dict(Gaussian_waist_fit_um=float(np.sqrt(-4*coeff[0])),Pin_down_W_per_m=I['down_power'],analytic_unit_amplitude_Gaussian_Pin_W_per_m=analytic,measured_vs_nominal_power_relative=I['down_power']/analytic-1,fit_method='ln|incident angular spectrum| vs kx^2 for amplitude>2% peak; Gaussian field exp(-w^2 kx^2/4)',interpretation='2D line Gaussian; 38um intensity1/e^2 diameter corresponds to waist19um')
    dump(ROOT/'incident_beam_audit.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__':main()
