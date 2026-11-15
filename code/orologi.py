"""Check D: is stellar activity really the "clock" preferred by the data?
Same source-sink model, but with alternative sources. Data: supernovae + DESI + compressed Planck priors (as in fit3.py).
Check E: sensitivity to the starting epoch and to high-redshift star formation.
Check F: energy budget (how much matter is converted into space per kg of stars)."""
import numpy as np, math, json, os
from scipy.optimize import minimize
from fit3 import *
def psi_o(o,b=2.9,c=5.6): return o**2.7/(1+(o/b)**c)
def past_gen(wm,al,ga,S,sink,x0=0.0,zstart=1e9):
    m=wm*(1+ZMAX)**3; x=x0; M=0.0; out=[(m,x,M)]
    def rhs(u,m,x,M):
        o=math.exp(u); E2=m+OR*o**4+x
        if E2<=0 or m<=0: raise ValueError
        E=math.sqrt(E2); s=S(o) if o-1<=zstart else 0.0
        sk=ga*x if sink=='prop' else ga*M if sink=='cum' else ga if sink=='const' else 0.0
        Q=al*s-sk
        return 3*m+Q/E, -Q/E, -psi_o(o)/E
    for i in range(NU):
        hh=UL[i+1]-UL[i]; u=UL[i]
        a1,b1,c1=rhs(u,m,x,M); a2,b2,c2=rhs(u+hh/2,m+hh*a1/2,x+hh*b1/2,M+hh*c1/2)
        a3,b3,c3=rhs(u+hh/2,m+hh*a2/2,x+hh*b2/2,M+hh*c2/2); a4,b4,c4=rhs(u+hh,m+hh*a3,x+hh*b3,M+hh*c3)
        m+=hh*(a1+2*a2+2*a3+a4)/6; x+=hh*(b1+2*b2+2*b3+b4)/6; M+=hh*(c1+2*c2+2*c3+c4)/6; out.append((m,x,M))
    o=np.array(out); z=np.expm1(U)
    return np.sqrt(o[:,0]+OR*(1+z)**4+o[:,1]), o
def make(S,sink,mode='src',zstart=1e9,tail=None):
    """mode 'src': p=(wm,wb,ln alpha,ln gamma);  mode 'ini': no source, p=(wm,wb,ln x_initial,ln gamma)."""
    def f(p):
        ga=math.exp(p[3]) if len(p)>3 else 0.0
        if mode=='src':
            al=math.exp(p[2]); x0=al*tail(p[0]) if tail else 0.0
            return past_gen(p[0],al,ga,S,sink,x0=x0,zstart=zstart)
        return past_gen(p[0],0.0,ga,S,sink,x0=math.exp(p[2]))
    return f
def tot(f,p,parts=False):
    try:
        if not(0.10<p[0]<0.20 and 0.018<p[1]<0.027): return 1e9
        e,o=f(p)
        if not np.all(np.isfinite(e)): return 1e9
        cb,csn,ccm,hh=observables(e,p[0],p[1])
    except (ValueError,FloatingPointError,OverflowError): return 1e9
    return (cb,csn,ccm,hh,o) if parts else cb+csn+ccm
def fit(f,starts,npar=4):
    best=None
    for s in starts:
        x0=[0.1435,0.02240]+list(s)[:npar-2]
        r=minimize(lambda v: tot(f,v),x0,method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-3,maxiter=1500,adaptive=True))
        r=minimize(lambda v: tot(f,v),r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-4,maxiter=1500,adaptive=True))
        if best is None or r.fun<best.fun: best=r
    return best
GRID=[(a,g) for a in (-6.5,-4.0,-1.5) for g in (-3.0,-1.0)]
if __name__=='__main__':
    R=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'results3.json'))); L=R['LCDM']['chi2']; out={}
    stel=lambda o: psi_o(o)*o**3
    def tail_inf(wm):   # space created before z=30 if the source is extrapolated to infinity (matter domination)
        from scipy.integrate import quad
        return quad(lambda o: psi_o(o)*o**0.5,1+ZMAX,np.inf)[0]/math.sqrt(wm)
    cases=[
     ('stellar, sink prop.',            make(stel,'prop'),4),
     ('stellar, constant sink',         make(stel,'const'),4),
     ('stellar, sink ∝ accumulated mass',make(stel,'cum'),4),
     ('no source: initial vacuum + sink prop.',   make(stel,'prop','ini'),4),
     ('no source: initial vacuum + constant sink',make(stel,'const','ini'),4),
     ('no source: initial vacuum + sink ∝ stellar mass',make(stel,'cum','ini'),4),
     ('source constant in time (p=0), sink prop.',   make(lambda o:1.0,'prop'),4),
     ('source constant in time (p=0), constant sink',make(lambda o:1.0,'const'),4),
     ('source ∝ (1+z)^1, sink prop.',  make(lambda o:o,'prop'),4),
     ('source ∝ (1+z)^2, sink prop.',  make(lambda o:o*o,'prop'),4),
     ('source ∝ matter density (1+z)^3, sink prop.',make(lambda o:o**3,'prop'),4),
     ('stellar with earlier peak (z~2.6), sink ∝ mass', make(lambda o:psi_o(o,3.6)*o**3,'cum'),4),
     ('stellar with later peak (z~1.2), sink ∝ mass',  make(lambda o:psi_o(o,2.2)*o**3,'cum'),4),
     ('E: stellar, source active only at z<10, sink ∝ mass', make(stel,'cum',zstart=10),4),
     ('E: stellar, source active only at z<20, sink ∝ mass', make(stel,'cum',zstart=20),4),
     ('E: stellar, extrapolated to infinite z, sink ∝ mass',make(stel,'cum',tail=tail_inf),4),
     ('E: stellare, calo piu\' ripido ad alto z (5.6->6.6)', make(lambda o:psi_o(o,2.9,6.6)*o**3,'cum'),4),
     ('E: stellare, calo meno ripido ad alto z (5.6->4.6)',  make(lambda o:psi_o(o,2.9,4.6)*o**3,'cum'),4),
    ]
    print('LCDM chi2 = %.2f'%L)
    for name,f,npar in cases:
        b=fit(f,GRID,npar); cb,csn,ccm,hh,o=tot(f,b.x,True)
        xi=lambda zz: float(np.interp(math.log1p(zz),U[::-1],o[::-1,1]))
        rr=[xi(zz)/xi(0) for zz in (0.5,1,2)]
        out[name]=dict(p=[float(a) for a in b.x],chi2=float(b.fun),d=float(b.fun-L),h=hh,Om0=float(o[-1,0]/hh**2),xr=rr)
        print('%-62s chi2 %.2f  (vs LCDM %+.2f)  h %.3f  Om0 %.3f  x(z)/x(0) at z=0.5,1,2: %.2f %.2f %.2f'%(name,b.fun,b.fun-L,hh,o[-1,0]/hh**2,*rr),flush=True)
    json.dump(out,open('orologi.json','w'),indent=1,ensure_ascii=False)
    # ---- Check F: energy budget for the stellar fit with sink ∝ accumulated mass
    p=out['stellar, sink ∝ accumulated mass']['p']; al=math.exp(p[2]); ga=math.exp(p[3]); f=make(stel,'cum'); e,o=f(p); hh=float(e[-1])
    M0=o[-1,2]; conv=al*M0                 # energy removed from matter per comoving volume (omega units)
    stelle=0.015*M0*9.7779e9/2.775e11      # stellar mass formed per comoving volume (omega units), Madau-Dickinson
    z=np.expm1(U); a3=(1+z)**-3
    back=np.trapezoid(ga*o[:,2]*a3/(-e), U) if hasattr(np,'trapezoid') else np.trapz(ga*o[:,2]*a3/(-e),U)   # energy returned by the sink: int sink a^3 dt, dt=-du/E
    print('\nF) budget per comoving volume, as a fraction of the critical density today:')
    print('   initial matter %.4f | stars formed %.4f | matter converted into space by the source %.4f (%.1f%% of matter)'%(p[0]/hh**2,stelle/hh**2,conv/hh**2,100*conv/p[0]))
    print('   energy returned to matter by the sink %.4f | extra space today %.3f | matter today %.3f'%(back/hh**2,o[-1,1]/hh**2,o[-1,0]/hh**2))
    print('   kg of matter converted into space per kg of stars formed: %.2f'%(conv/stelle))
    for diss in (0.0010,0.0022,0.0035): print('   if each kg of stars dissipates %.4f c^2: trigger factor = %.0f'%(diss,conv/stelle/diss))
