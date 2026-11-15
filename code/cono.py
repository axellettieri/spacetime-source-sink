"""Energia oscura come media degli eventi nel cono di luce passato (traduzione dell'ipotesi nel linguaggio dei causal set).
Nei causal set la costante cosmologica e' coniugata al volume quadridimensionale V del cono di luce passato:  Lambda ~ S/V,
con S somma di contributi dei singoli eventi (casuali in Sorkin). Qui S conta gli eventi della materia:
    rho_x(t) = beta * <n_sorgente>_cono(t)  -  gamma * <n_pozzo>_cono(t)
dove <n>_cono e' la densita' di eventi mediata sul volume quadridimensionale del cono di luce passato dell'osservatore al tempo t.
Nessuna funzione libera: il nucleo che rende 'diffuso' l'effetto di eventi locali e' fissato dalla geometria.
L'energia e' scambiata con la materia (vuoto interagente), come in fit3.py. Dati: supernovae + DESI + prior compressi di Planck."""
import numpy as np, math, json, sys, os
from scipy.optimize import minimize
from fit3 import *
NA=320; LNA=np.linspace(math.log(1/1101.0),0.0,NA); A=np.exp(LNA); Z=1/A-1; DL=LNA[1]-LNA[0]
def psi_z(z): return (1+z)**2.7/(1+((1+z)/2.9)**5.6)
ZCUT=30.0
N_STAR=np.where(Z<=ZCUT,psi_z(Z)*(1+Z)**3,0.0)            # densita' fisica di eventi dissipativi (attivita' stellare)
def solve(wm,beta,gamma=0.0,src='stelle',sink=None,iters=60,full=False):
    x=np.full(NA,0.32); E=np.sqrt(wm/A**3+OR/A**4+x); Mc=np.full(NA,wm)
    for it in range(iters):
        dt=DL/E                                             # dt per passo in ln a
        chi=np.concatenate([np.cumsum((dt/A)[::-1])[::-1][1:],[0.0]])   # distanza comovente da oggi fino ad a
        if src=='stelle': n=N_STAR
        elif src=='materia': n=Mc/A**3                       # ogni cambio di stato: proporzionale alla densita' di energia della materia
        mstar=np.cumsum(np.where(Z<=ZCUT,psi_z(Z),0.0)*dt)   # massa compatta accumulata (comovente)
        ns=None if sink is None else (mstar/A**3 if sink=='diluito' else mstar)
        D=chi[None,:]-chi[:,None]                            # D[i,j] = distanza comovente tra emettitore j e osservatore i (j<=i)
        W=np.where(D>0,(A[None,:]*D)**3*dt[None,:],0.0)       # volume quadridimensionale del cono: raggio fisico^3 * dt
        V=W.sum(1); V[V==0]=1e-300
        xn=beta*(W@n)/V
        if ns is not None: xn=xn-gamma*(W@ns)/V
        xn[0]=xn[1] if NA>1 else 0.0
        x=0.5*x+0.5*xn
        dx=np.diff(x); Mc=wm-np.concatenate([[0],np.cumsum(0.5*(A[1:]**3+A[:-1]**3)*dx)])
        E2=Mc/A**3+OR/A**4+x
        if np.any(E2<=0) or np.any(Mc<=0): raise ValueError
        En=np.sqrt(E2); conv=np.max(np.abs(En/E-1)); E=En
        if it>8 and conv<1e-7: break
    if full: return E,x,Mc,mstar
    return E,x,Mc
def model(p,src,sink):
    wm=p[0]; beta=math.exp(p[2]); gamma=math.exp(p[3]) if sink else 0.0
    E,x,Mc=solve(wm,beta,gamma,src,sink)
    e_u=np.interp(-U,LNA,E)                                  # griglia di fit3 (z da 30 a 0)
    return e_u,E,x,Mc
def tot(p,src,sink,parts=False):
    try:
        if not(0.10<p[0]<0.20 and 0.018<p[1]<0.027): return 1e9
        e_u,E,x,Mc=model(p,src,sink)
        cb,csn,ccm,hh=observables(e_u,p[0],p[1])
    except (ValueError,FloatingPointError,OverflowError): return 1e9
    return (cb,csn,ccm,hh,x,Mc) if parts else cb+csn+ccm
def fit(src,sink,starts):
    best=None
    for s in starts:
        r=minimize(lambda v:tot(v,src,sink),s,method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-3,maxiter=500,adaptive=True))
        r=minimize(lambda v:tot(v,src,sink),r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-4,maxiter=500,adaptive=True))
        if best is None or r.fun<best.fun: best=r
    return best
if __name__=='__main__':
    suf='' if SNSET=='pantheon' else '_'+SNSET
    R=json.load(open('results3%s.json'%suf)); L=R['LCDM']['chi2']
    print('campione supernovae:',SNSET,' LCDM chi2 = %.2f ; w0wa %+.2f ; sorgente-pozzo (cum) %+.2f'%(L,R['w0wa']['chi2']-L,R['cum_n3']['chi2']-L))
    # normalizzazione di partenza: beta tale che x(oggi) ~ 0.32
    E,x,Mc=solve(0.1434,1e-7); b0=math.log(0.32/(x[-1]/1e-7)); print('media sul cono oggi: %.2f ; profilo x(z)/x(0) a z=0.5,1,2,3,5: '%(x[-1]/1e-7)+' '.join('%.2f'%(np.interp(-math.log1p(q),LNA,x)/x[-1]) for q in (0.5,1,2,3,5)))
    out={}
    cases=[('cono: sola sorgente stellare (0 parametri in piu\')','stelle',None,[[0.1434,0.0224,b0],[0.146,0.0223,b0+0.2]]),
           ('cono: sorgente stellare + pozzo diluito (1 in piu\')','stelle','diluito',[[0.1434,0.0224,b0+0.3,b0-1.0],[0.1434,0.0224,b0+0.1,b0-3.0],[0.1434,0.0224,b0+0.6,b0+0.5]]),
           ('cono: sorgente = tutta la materia (0 in piu\')','materia',None,[[0.1434,0.0224,-3.0],[0.1434,0.0224,-6.0]]),
           ('cono: sorgente stellare + pozzo non diluito (1 in piu\')','stelle','nondiluito',[[0.1434,0.0224,b0+0.3,b0-1.0],[0.1434,0.0224,b0+0.1,b0-3.0],[0.1434,0.0224,b0+0.6,b0+0.5]])]
    for name,src,sink,starts in cases:
        b=fit(src,sink,starts); cb,csn,ccm,hh,x,Mc=tot(b.x,src,sink,True)
        xr=[float(np.interp(-math.log1p(zq),LNA,x)/x[-1]) for zq in (0.25,0.5,1,2,3,5)]
        out[name]=dict(p=[float(v) for v in b.x],chi2=float(b.fun),d=float(b.fun-L),h=hh,Om0=float(Mc[-1]/hh**2))
        print('%-58s chi2 %.2f (da LCDM %+.2f) | bao %.1f sn %.1f cmb %.1f | H0 %.1f materia %.3f | x(z)/x(0) a z=0.25,0.5,1,2,3,5: %s'%(name,b.fun,b.fun-L,cb,csn,ccm,100*hh,Mc[-1]/hh**2,' '.join('%.2f'%v for v in xr)),flush=True)
    json.dump(out,open('cono%s.json'%suf,'w'),indent=1,ensure_ascii=False)
