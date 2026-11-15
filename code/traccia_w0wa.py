"""The forced relation between w0 and wa predicted by the source-sink model.
For each sink strength (gamma) the other parameters are fitted to the data (quick comparison), then the (w0, wa) pair
that best reproduces the model's expansion history is found (same initial matter, same CMB distance).
The model has one parameter fewer than w0wa: its points lie on a line in the (w0, wa) plane, not on the whole plane.
Usage: python traccia_w0wa.py [cum|prop|const]"""
import sys, json, math, numpy as np
from scipy.optimize import minimize
from fit3 import *
kind=sys.argv[1] if len(sys.argv)>1 and sys.argv[1] in ('cum','prop','const') else 'cum'
z=np.expm1(U); msk=z<=2.5; du=np.gradient(U)
def dist(E): return float(np.sum(np.exp(U)/E*np.abs(du)))
def mappa(E,wm):
    D0=dist(E)
    def f(v):
        try: Ew=e_fluid(wm,v[0],v[1],v[2])
        except Exception: return 1e9
        if not np.all(np.isfinite(Ew)): return 1e9
        return float(np.mean((Ew[msk]/E[msk]-1)**2))*1e6+((dist(Ew)/D0-1)*1e4)**2
    best=None
    for s in ([0.31,-0.9,-0.3],[0.31,-0.7,-1.0]):
        r=minimize(f,s,method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-6,maxiter=4000))
        if best is None or r.fun<best.fun: best=r
    Ew=e_fluid(wm,*best.x); return best.x[1],best.x[2],float(np.max(np.abs(Ew[msk]/E[msk]-1)))
if __name__=='__main__':
    R=json.load(open('results3.json')); L0=R['LCDM']['chi2']; p0=R[kind+'_n3']['p']
    print('sink:',kind,'| best: chi2 vs LCDM %+.2f'%(R[kind+'_n3']['chi2']-L0))
    print(' ln(gamma) | chi2 vs LCDM | H0  | equivalent w0 | equivalent wa | maximum deviation of H(z) up to z=2.5')
    out=[]
    for dg in (-3,-2,-1.5,-1,-0.6,-0.3,0,0.3,0.6,1.0,1.5):
        lg=p0[3]+dg
        def f(v):
            try:
                if not(0.10<v[0]<0.20 and 0.018<v[1]<0.027): return 1e9
                E,o=past_interacting(v[0],math.exp(v[2]),math.exp(lg),3,kind); return sum(observables(E,v[0],v[1])[:3])
            except (ValueError,FloatingPointError,OverflowError): return 1e9
        r=minimize(f,[p0[0],p0[1],p0[2]],method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-3,maxiter=1500))
        r=minimize(f,r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-4,maxiter=1500))
        E,o=past_interacting(r.x[0],math.exp(r.x[2]),math.exp(lg),3,kind); hh=observables(E,r.x[0],r.x[1])[3]
        w0,wa,err=mappa(E,r.x[0]); out.append((lg,r.fun-L0,w0,wa))
        print('  %+6.2f   |   %+6.2f     | %.1f |    %+.3f      |    %+.3f      |   %.2f%%'%(lg,r.fun-L0,100*hh,w0,wa,100*err),flush=True)
    json.dump(out,open('traccia_w0wa_%s.json'%kind,'w'))
    w=R['w0wa']['p']; print('for comparison, best free w0wa on the same data (quick comparison): w0 %+.3f wa %+.3f, chi2 vs LCDM %+.2f'%(w[3],w[4],R['w0wa']['chi2']-L0))
