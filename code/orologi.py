"""Controllo D: l'attivita' stellare e' davvero l'"orologio" preferito dai dati?
Stesso modello sorgente-pozzo, ma con sorgenti alternative. Dati: supernovae + DESI + prior compressi di Planck (come fit3.py).
Controllo E: sensibilita' all'epoca di inizio e alla formazione stellare ad alto redshift.
Controllo F: bilancio dell'energia (quanta materia viene convertita in spazio per ogni kg di stelle)."""
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
    """mode 'src': p=(wm,wb,ln alfa,ln gamma);  mode 'ini': nessuna sorgente, p=(wm,wb,ln x_iniziale,ln gamma)."""
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
    def tail_inf(wm):   # spazio creato prima di z=30 se la sorgente viene estrapolata all'infinito (dominio della materia)
        from scipy.integrate import quad
        return quad(lambda o: psi_o(o)*o**0.5,1+ZMAX,np.inf)[0]/math.sqrt(wm)
    cases=[
     ('stellare, pozzo prop.',            make(stel,'prop'),4),
     ('stellare, pozzo costante',         make(stel,'const'),4),
     ('stellare, pozzo ∝ massa accumulata',make(stel,'cum'),4),
     ('nessuna sorgente: vuoto iniziale + pozzo prop.',   make(stel,'prop','ini'),4),
     ('nessuna sorgente: vuoto iniziale + pozzo costante',make(stel,'const','ini'),4),
     ('nessuna sorgente: vuoto iniziale + pozzo ∝ massa stellare',make(stel,'cum','ini'),4),
     ('sorgente costante nel tempo (p=0), pozzo prop.',   make(lambda o:1.0,'prop'),4),
     ('sorgente costante nel tempo (p=0), pozzo costante',make(lambda o:1.0,'const'),4),
     ('sorgente ∝ (1+z)^1, pozzo prop.',  make(lambda o:o,'prop'),4),
     ('sorgente ∝ (1+z)^2, pozzo prop.',  make(lambda o:o*o,'prop'),4),
     ('sorgente ∝ densita\' di materia (1+z)^3, pozzo prop.',make(lambda o:o**3,'prop'),4),
     ('stellare con picco anticipato (z~2.6), pozzo ∝ massa', make(lambda o:psi_o(o,3.6)*o**3,'cum'),4),
     ('stellare con picco ritardato (z~1.2), pozzo ∝ massa',  make(lambda o:psi_o(o,2.2)*o**3,'cum'),4),
     ('E: stellare, sorgente attiva solo da z<10, pozzo ∝ massa', make(stel,'cum',zstart=10),4),
     ('E: stellare, sorgente attiva solo da z<20, pozzo ∝ massa', make(stel,'cum',zstart=20),4),
     ('E: stellare, estrapolata fino a z infinito, pozzo ∝ massa',make(stel,'cum',tail=tail_inf),4),
     ('E: stellare, calo piu\' ripido ad alto z (5.6->6.6)', make(lambda o:psi_o(o,2.9,6.6)*o**3,'cum'),4),
     ('E: stellare, calo meno ripido ad alto z (5.6->4.6)',  make(lambda o:psi_o(o,2.9,4.6)*o**3,'cum'),4),
    ]
    print('LCDM chi2 = %.2f'%L)
    for name,f,npar in cases:
        b=fit(f,GRID,npar); cb,csn,ccm,hh,o=tot(f,b.x,True)
        xi=lambda zz: float(np.interp(math.log1p(zz),U[::-1],o[::-1,1]))
        rr=[xi(zz)/xi(0) for zz in (0.5,1,2)]
        out[name]=dict(p=[float(a) for a in b.x],chi2=float(b.fun),d=float(b.fun-L),h=hh,Om0=float(o[-1,0]/hh**2),xr=rr)
        print('%-62s chi2 %.2f  (da LCDM %+.2f)  h %.3f  Om0 %.3f  x(z)/x(0) a z=0.5,1,2: %.2f %.2f %.2f'%(name,b.fun,b.fun-L,hh,o[-1,0]/hh**2,*rr),flush=True)
    json.dump(out,open('orologi.json','w'),indent=1,ensure_ascii=False)
    # ---- Controllo F: bilancio dell'energia per il fit stellare con pozzo ∝ massa accumulata
    p=out['stellare, pozzo ∝ massa accumulata']['p']; al=math.exp(p[2]); ga=math.exp(p[3]); f=make(stel,'cum'); e,o=f(p); hh=float(e[-1])
    M0=o[-1,2]; conv=al*M0                 # energia tolta alla materia per volume comovente (unita' omega)
    stelle=0.015*M0*9.7779e9/2.775e11      # massa di stelle formata per volume comovente (unita' omega), Madau-Dickinson
    z=np.expm1(U); a3=(1+z)**-3
    back=np.trapezoid(ga*o[:,2]*a3/(-e), U) if hasattr(np,'trapezoid') else np.trapz(ga*o[:,2]*a3/(-e),U)   # energia restituita dal pozzo: int sink a^3 dt, dt=-du/E
    print('\nF) bilancio per volume comovente, in frazione della densita\' critica di oggi:')
    print('   materia iniziale %.4f | stelle formate %.4f | materia convertita in spazio dalla sorgente %.4f (%.1f%% della materia)'%(p[0]/hh**2,stelle/hh**2,conv/hh**2,100*conv/p[0]))
    print('   energia restituita alla materia dal pozzo %.4f | spazio in piu\' oggi %.3f | materia oggi %.3f'%(back/hh**2,o[-1,1]/hh**2,o[-1,0]/hh**2))
    print('   kg di materia convertiti in spazio per ogni kg di stelle formate: %.2f'%(conv/stelle))
    for diss in (0.0010,0.0022,0.0035): print('   se ogni kg di stelle dissipa %.4f c^2: fattore di innesco = %.0f'%(diss,conv/stelle/diss))
