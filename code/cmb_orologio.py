"""Full Planck (plik-lite TTTEEE + low multipoles) + DESI + supernovae for the source-sink model
with the EMERGENT CLOCK: events are collapsing matter (from the 3D simulation or from theory), not the observed star formation.
  created space   = alpha * matter * (collapse rate)
  absorbed space  = gamma * omega_m * (collapsed fraction)          ('originale' formula of universo3d.py)
Usage: python cmb_orologio.py <clock .json file> <threshold> [iterations]      (threshold = value of 'dc' or of one of the 'soglie' in the file)
Same pipeline as cmb_full.py: effective fluid in CAMB, r_d from CAMB, lensing correction from the growth equations (diffuse exchange)."""
import sys, json, math, time, os
import numpy as np
from scipy.interpolate import PchipInterpolator
_arg=sys.argv[:]; sys.argv=sys.argv[:1]
from cmb_full import *
from cmb_full import _grow
sys.argv=_arg
def carica_orologio(f,soglia):
    J=json.load(open(f)); S=np.array(J['storia']); a=S[:,0]; k=(a<=1.0000001)&((S[:,1]>0)|(S[:,1].max()<=0))
    col=4 if abs(J['par']['dc']-soglia)<1e-9 else 8+[abs(s-soglia)<1e-9 for s in J.get('soglie',[])].index(True)
    u=-np.log(a[k]); F=S[k,col]; i=np.argsort(u); return PchipInterpolator(u[i],F[i])
def passato(wm,al,ga,Fs):
    dF=Fs.derivative()
    m=wm*(1+ZMAX)**3; x=0.0; out=[(m,x)]
    def rhs(u,m,x):
        o=math.exp(u); E2=m+OR*o**4+x
        if E2<=0 or m<=0: raise ValueError
        E=math.sqrt(E2); QE=-al*m*float(dF(u))-ga*wm*float(Fs(u))/E          # Q/E ; dF/du < 0 because u decreases in time
        return 3*m+QE, -QE
    for i in range(NU):
        hh=UL[i+1]-UL[i]; u=UL[i]
        a1,b1=rhs(u,m,x); a2,b2=rhs(u+hh/2,m+hh*a1/2,x+hh*b1/2); a3,b3=rhs(u+hh/2,m+hh*a2/2,x+hh*b2/2); a4,b4=rhs(u+hh,m+hh*a3,x+hh*b3)
        m+=hh*(a1+2*a2+2*a3+a4)/6; x+=hh*(b1+2*b2+2*b3+b4)/6; out.append((m,x))
    o=np.array(out); E=np.sqrt(o[:,0]+OR*(1+ZU)**4+o[:,1])
    Q=E*(-al*o[:,0]*dF(U))-ga*wm*Fs(U)
    return E,o,Q
def alens_orol(wm,e_u,o,Q):
    m=o[:,0]; g=-Q/m
    dM=_grow(e_u,m,g); dE=_grow(e_u,wm*(1+ZU)**3,0*g)
    lna=np.linspace(NLN[0],0,400); Mc=CubicSpline(NLN,m*(1+ZU)**-3)(lna)
    R=CubicSpline(lna,Mc*dM(lna)[0]/(wm*dE(lna)[0]))(-np.log1p(KZ)); return float(np.sum(KW*R**2))
def chi2_orol(p,Fs,parts=False):
    wm,wb,la,lg,logA,ns,tau=p
    try:
        if not(0.12<wm<0.17) or not(0.020<wb<0.025) or not(0.9<ns<1.05) or not(0.01<tau<0.12): return 1e9
        e,o,Q=passato(wm,math.exp(la),math.exp(lg),Fs)
        al=alens_orol(wm,e,o,Q); cl,r=spectra(e,wm,wb,logA,ns,tau,al); rd=float(r.get_derived_params()['rdrag'])
        cb,csn,_,hh=observables(e,wm,wb,rd=rd); cpl=-2*PL.loglike(cl[2:,0],cl[2:,3],cl[2:,1],2)
    except Exception as ex: return 1e9
    ctau=((tau-0.0544)/0.0073)**2
    if parts: return dict(planck=float(cpl),tau=float(ctau),bao=float(cb),sn=float(csn),tot=float(cpl+ctau+cb+csn),h=hh,alens=al,rd=rd,Om0=float(o[-1,0]/hh**2),x0=float(o[-1,1]))
    return cpl+ctau+cb+csn
if __name__=='__main__':
    f=sys.argv[1]; soglia=float(sys.argv[2]); nit=int(sys.argv[3]) if len(sys.argv)>3 else 500
    Fs=carica_orologio(f,soglia); suf='' if SNSET=='pantheon' else '_'+SNSET
    # 1) start: quick fit (compressed Planck) with the same formula
    L3=json.load(open(os.path.join(HERE,'results3%s.json'%suf)))['LCDM']['chi2']
    def rapido(v):
        try:
            if not(0.10<v[0]<0.20 and 0.018<v[1]<0.027): return 1e9
            e,o,Q=passato(v[0],math.exp(v[2]),math.exp(v[3]),Fs); return sum(observables(e,v[0],v[1])[:3])
        except (ValueError,FloatingPointError,OverflowError): return 1e9
    b0=math.log(0.33/(0.143*abs(np.trapezoid(Fs.derivative()(U)*np.exp(3*U),U))))
    best=None
    for g in (-0.5,0.3,1.0):
        r=minimize(rapido,[0.1445,0.02236,b0+0.3,g],method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-3,maxiter=3000,adaptive=True))
        r=minimize(rapido,r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-4,maxiter=3000,adaptive=True))
        if best is None or r.fun<best.fun: best=r
    print('quick fit (compressed Planck): chi2 vs LCDM %+.2f, p=%s'%(best.fun-L3,np.round(best.x,5).tolist()),flush=True)
    # 2) full Planck
    LC=json.load(open(os.path.join(HERE,'cmb2_LCDM_U%s.json'%suf)))
    x0=np.array([0.1425,0.0224,best.x[2],best.x[3],3.0509,0.9665,0.0589]); sc=np.array([3e-4,6e-5,0.02,0.05,0.004,0.002,0.003]); N_=7; t0=time.time(); bb=[1e9,None]; n=[0]
    def F(u):
        v=x0+u*sc; c=chi2_orol(list(v),Fs); n[0]+=1
        if c<bb[0]: bb[0]=c; bb[1]=list(v); print('n=%d t=%.0fs chi2=%.2f (vs LCDM %+.2f) p=%s'%(n[0],time.time()-t0,c,c-LC['tot'],np.round(v,5).tolist()),flush=True)
        return c
    minimize(F,np.zeros(N_),method='Nelder-Mead',options=dict(initial_simplex=np.vstack([np.zeros(N_),np.eye(N_)]),xatol=0.03,fatol=0.01,maxiter=nit,adaptive=True))
    out=chi2_orol(bb[1],Fs,True); out['p']=[float(a) for a in bb[1]]; out['sec']=time.time()-t0; out['delta']=out['tot']-LC['tot']; out['orologio']=os.path.basename(f); out['soglia']=soglia
    print('FINAL',json.dumps(out),flush=True); json.dump(out,open(os.path.join(HERE,'cmb_orologio_%s_%g%s.json'%(os.path.basename(f).replace('.json',''),soglia,suf)),'w'))
