"""Independent physical-x field monitor audit of Fourier phase origin."""
import numpy as np
from run_rcwa import ROOT,session,archive_dataset,write_json,PERIOD,WAVELENGTH,REFERENCE

def simulate():
    with session() as s:
        s.load(str(ROOT/'harmonics_160_dx1.25nm'/'model.fsp'))
        s.switchtolayout();s.addrcwafieldmonitor();s.set('name','origin_audit')
        s.set('monitor type','Linear X')
        s.set('x min',-1e-6);s.set('x max',1e-6);s.set('x number points',1601)
        s.set('y',0);s.set('z',-100e-9)
        s.cd(str(ROOT));s.save('origin_audit.fsp');s.run();s.save('origin_audit.fsp')
        archive_dataset(s.getresult('origin_audit','Es'),ROOT,'spatial_origin_Es')

def project():
    with np.load(ROOT/'spatial_origin_Es.npz') as a:
        x=a['x'].ravel();E=a['E'][:,0,0,0,0,1];z=float(a['z'].ravel()[0])
    k0=2*np.pi/WAVELENGTH;m=np.arange(-3,4);kx=2*np.pi/PERIOD*m
    kz=np.sqrt(k0*k0-kx*kx)
    C=np.array([np.trapezoid(E*np.exp(-1j*k*x),x)/PERIOD for k in kx])
    # Native source s has Ey=-1; subtract incident before r0 projection.
    C[3]+=np.exp(1j*k0*(z-REFERENCE))
    r=-C*np.exp(1j*(kz*z+k0*REFERENCE))
    with np.load(ROOT/'reference_for_analysis.npz') as a:
        native_deembedded=a['r_before_lateral_correction']
    expected=native_deembedded*np.exp(+1j*np.pi*m)
    error=float(max(abs(r-expected)))
    np.savez_compressed(ROOT/'spatial_origin_projection.npz',m=m,r_projected=r,
        r_expected_xmin_origin=expected,C=C,z=np.array(z))
    write_json(ROOT/'spatial_origin_projection.json',{
        'method':'trapezoidal integral Ey(x)*exp(-i*kxm*x)/period at physical monitor coordinates',
        'x_range_nm':[-1000,1000],'points':len(x),'z_nm':z*1e9,
        'deembedding':'r=-C*exp(i*(kzm*zmonitor+k0*zreference)); subtract incident for m=0',
        'prediction':'r_physical_x0=r_raw_z_deembedded*exp(+i*m*pi)',
        'max_complex_error':error,
        'conclusion':'Native RCWA field_amplitudes phase starts at x_min, not at physical x=0.',
        'half_mesh_note':'This test deliberately uses uncorrected coefficients; the finite-mesh lateral phase bias is audited separately.'})
    assert error<1e-10,error
    print('Physical-x origin projection max complex error:',error,flush=True)

if __name__=='__main__':
    if not (ROOT/'spatial_origin_Es.npz').exists():simulate()
    project()
