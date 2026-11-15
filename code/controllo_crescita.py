"""Check C: how wrong is the 'effective fluid' approximation used in the CMB test?
Compares the source of the gravitational potential  Delta = (comoving matter) x (contrast)  in three treatments:
  E = effective fluid (CAMB): the initial matter clusters fully, the rest is smooth
  U = uniform exchange (matter removed and returned diffusely)  -> 'diffuse' version of the model
  L = local exchange (matter removed and returned where matter already is)"""
import numpy as np, math, json
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from fit3 import *
HERE_=os.path.dirname(os.path.abspath(__file__))
N=-U; z=np.expm1(U)
def solveD(E,m,g):
    sE=CubicSpline(N,E); sm=CubicSpline(N,m); sg=CubicSpline(N,g)
    def f(n,y):
        e=sE(n); de=sE(n,1); gg=sg(n); dg=sg(n,1)
        return [y[1], -(de/e+2+gg/e)*y[1]+(1.5*sm(n)-2*e*gg-e*dg)/e**2*y[0]]
    a0=math.exp(N[0]); s=solve_ivp(f,[N[0],0],[a0,a0],dense_output=True,rtol=1e-9,atol=1e-13)
    return s.sol
def Delta_curves(p,kind,n=3):
    al=math.exp(p[2]); ga=math.exp(p[3]); E,o=past_interacting(p[0],al,ga,n,kind); m,x,M=o.T
    sink=ga*x if kind=='prop' else ga*M if kind=='cum' else ga*np.ones_like(x)
    Q=al*psi(z)*(1+z)**n-sink; g=-Q/m
    a3=(1+z)**-3; Mc=m*a3            # comoving matter (in omega units)
    me=p[0]*(1+z)**3
    dU=solveD(E,m,g); dL=solveD(E,m,0*g); dE=solveD(E,me,0*g)
    sM=CubicSpline(N,Mc)
    return dict(E=E,Mc=sM,dU=dU,dL=dL,dE=dE,we=p[0],g=CubicSpline(N,g/E))
if __name__=='__main__':
    R=json.load(open(os.path.join(HERE_,'results3.json'))); f2=os.path.join(HERE_,'cmb2_cum_n3_U.json'); C=json.load(open(f2)) if os.path.exists(f2) else {'p':R['cum_n3']['p']+[0,0,0]}
    pl=R['LCDM']['p']; El=e_fluid(pl[0],pl[2]); dl=solveD(El,pl[0]*(1+z)**3,0*z)
    for lab,p in (('cum_n3 (fit with full CMB)',C['p'][:4]),('cum_n3 (fit with compressed priors)',R['cum_n3']['p'])):
        c=Delta_curves(p,'cum'); print('==',lab)
        print('   z    matter/initial     g/H     Delta_U/Delta_E  Delta_L/Delta_E   Delta_U/LCDM  Delta_L/LCDM  Delta_E/LCDM')
        for zz in (0,0.25,0.5,1,1.5,2,3,5,10):
            n=-math.log1p(zz); Mc=float(c['Mc'](n)); DE=c['we']*c['dE'](n)[0]; DU=Mc*c['dU'](n)[0]; DL=Mc*c['dL'](n)[0]; DLC=pl[0]*dl(n)[0]
            print('  %4.2f     %.3f         %+.3f      %.3f           %.3f            %.3f        %.3f        %.3f'%(zz,Mc/c['we'],float(c['g'](n)),DU/DE,DL/DE,DU/DLC,DL/DLC,DE/DLC))
