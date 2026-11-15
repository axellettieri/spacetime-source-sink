"""Linear growth of structures in the source-sink model and comparison with f*sigma8 (RSD) and S8 (lensing).
Two hypotheses on where the energy returned to matter ends up:
  'uniforme' : distributed uniformly (dilutes the contrasts)
  'locale'   : in proportion to the matter already present (does not change the contrasts)"""
import numpy as np, math, json
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from fit3 import *
R=json.load(open('results3.json'))
# f*sigma8: 6dFGS (Beutler 2012), SDSS MGS, BOSS DR12, eBOSS DR16 (Alam et al. 2021)
FS=[(0.067,0.423,0.055),(0.15,0.53,0.16),(0.38,0.497,0.045),(0.51,0.459,0.038),(0.70,0.473,0.041),(0.85,0.315,0.095),(1.48,0.462,0.045)]
N=-U  # ln a increasing
def grow(E,m,g):
    """E,m,g on the U grid. Returns splines of delta and of the velocity (theta) as functions of ln a."""
    sE=CubicSpline(N,E); sm=CubicSpline(N,m); sg=CubicSpline(N,g)
    def f(n,y):
        e=sE(n); de=sE(n,1); gg=sg(n); dg=sg(n,1)
        return [y[1], -(de/e+2+gg/e)*y[1]+(1.5*sm(n)-2*e*gg-e*dg)/e**2*y[0]]
    a0=math.exp(N[0]); s=solve_ivp(f,[N[0],0],[a0,a0],dense_output=True,rtol=1e-8,atol=1e-12)
    return s.sol,sE,sg
out={}
def report(name,E,m,g,hh,Om0,ref=None):
    sol,sE,sg=grow(E,m,g); d0=sol(0)[0]
    out[name]=dict(D0=float(d0),h=hh,Om0=Om0)
    if ref is None: ref=d0
    s8=0.811*d0/ref; out[name]['sigma8']=float(s8); out[name]['S8']=float(s8*math.sqrt(Om0/0.3))
    c2=0; pts=[]
    for z,v,e in FS:
        n=-math.log1p(z); d,dp=sol(n); fs=(dp+sg(n)/sE(n)*d)/d0*s8; c2+=((fs-v)/e)**2; pts.append(round(float(fs),3))
    out[name]['fs8']=pts; out[name]['chi2_fs8']=float(c2)
    print(name, 'sigma8=%.3f S8=%.3f Om0=%.3f chi2_fs8=%.1f'%(s8,out[name]['S8'],Om0,c2), pts)
    return d0
z=np.expm1(U)
p=R['LCDM']['p']; E=e_fluid(p[0],p[2]); h=R['LCDM']['h']
ref=report('LCDM',E,p[0]*(1+z)**3,0*z,h,p[0]/h**2)
p=R['w0wa']['p']; E=e_fluid(p[0],p[2],p[3],p[4]); h=R['w0wa']['h']
report('w0wa',E,p[0]*(1+z)**3,0*z,h,p[0]/h**2,ref)
for name in ('prop_n3','cum_n3','const_n3'):
    p=R[name]['p']; kind=name.split('_')[0]; al=math.exp(p[2]); ga=math.exp(p[3])
    E,o=past_interacting(p[0],al,ga,3,kind); m,x,M=o.T; h=R[name]['h']
    sink=ga*x if kind=='prop' else ga*M if kind=='cum' else ga*np.ones_like(x)
    Q=al*psi(z)*(1+z)**3-sink; g=-Q/m     # g = rate at which matter gains energy, per unit matter
    report(name+'_uniforme',E,m,g,h,float(m[-1]/h**2),ref)
    report(name+'_locale',E,m,0*g,h,float(m[-1]/h**2),ref)
    print('   g/H today = %.3f ; at z=0.5: %.3f ; at z=2: %.3f'%(g[-1]/E[-1], np.interp(math.log(1.5),U[::-1],(g/E)[::-1]), np.interp(math.log(3),U[::-1],(g/E)[::-1])))
print('fs8 data:',[v for _,v,_ in FS],'errors',[e for *_,e in FS])
json.dump(out,open('growth.json','w'),indent=1)
