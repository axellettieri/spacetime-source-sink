"""Adatta ai dati (supernovae + DESI + Planck compresso) il modello in cui l'orologio degli eventi viene dalla simulazione 3D.
Uso: python adatta3d.py orologio.json [originale]   -> stampa il confronto con LCDM e scrive parametri3d.json (alpha, gamma, wm, wb)
Stesse formule di universo3d.py: spazio creato = alpha * massa collassata nuova ; assorbito = gamma * massa collassata * tempo."""
import numpy as np, math, json, sys
from scipy.optimize import minimize
from fit3 import *
J=json.load(open(sys.argv[1])); S=np.array(J['storia']); A_=S[:,0]; F_=S[:,4]; LA=np.log(A_)
POZZO=sys.argv[2] if len(sys.argv)>2 else 'materia'      # 'originale' = assorbimento che non si autoalimenta
OROLOGI=[(J['par']['dc'],S[:,4])]+[(th,S[:,8+j]) for j,th in enumerate(J.get('soglie',[]))]
def fondo(wm,alpha,gamma):
    global F_
    x=0.0; Mc=wm; E=np.zeros(len(A_)); X=np.zeros(len(A_))
    for i,a in enumerate(A_):
        E2=Mc/a**3+OR/a**4+x
        if E2<=0: raise ValueError
        E[i]=math.sqrt(E2); X[i]=x
        if i+1<len(A_):
            da=A_[i+1]-a; dt=da/(a*E[i]); dx=alpha*Mc*(F_[i+1]-F_[i])/a**3-gamma*(Mc if POZZO=='materia' else wm)*F_[i]*dt
            x+=dx; Mc-=a**3*dx
            if Mc<=0: raise ValueError
    return E,X,Mc
def tot(p,sink=True,parts=False):
    try:
        if not(0.10<p[0]<0.20 and 0.018<p[1]<0.027): return 1e9
        E,X,Mc=fondo(p[0],math.exp(p[2]),math.exp(p[3]) if sink else 0.0)
        cb,csn,ccm,hh=observables(np.interp(-U,LA,E),p[0],p[1])
    except (ValueError,FloatingPointError,OverflowError): return 1e9
    return (cb,csn,ccm,hh,X,Mc) if parts else cb+csn+ccm
def fit(sink,starts):
    best=None
    for s in starts:
        r=minimize(lambda v:tot(v,sink),s,method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-3,maxiter=3000,adaptive=True))
        r=minimize(lambda v:tot(v,sink),r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-4,maxiter=3000,adaptive=True))
        if best is None or r.fun<best.fun: best=r
    return best
if __name__=='__main__':
    suf='' if SNSET=='pantheon' else '_'+SNSET
    R=json.load(open('results3%s.json'%suf)); L0=R['LCDM']['chi2']
    print('LCDM chi2 = %.2f ; w0wa %+.2f ; sorgente-pozzo con orologio stellare %+.2f'%(L0,R['w0wa']['chi2']-L0,R['cum_n3']['chi2']-L0))
    migliore=None
    for dc,Fc in sorted(OROLOGI):
        F_=Fc
        if F_[-1]<=0: continue
        b0=math.log(0.33/(0.143*np.sum(np.diff(F_)/A_[:-1]**3)))
        r1=fit(False,[[0.1434,0.0224,b0,0],[0.146,0.0223,b0+0.2,0]])
        r2=fit(True,[[0.1434,0.0224,b0+0.3,g] for g in (-1.0,0.5,1.5)]+[[0.1434,0.0224,b0+1.0,2.0]])
        cb,csn,ccm,hh,X,Mc=tot(r2.x,True,True)
        mc=J['par'].get('wm',0.143)*2.775e11*(J['par']['L']/J['par']['ng'])**3*dc if 'L' in J['par'] else 0
        print('soglia di collasso %5.1f (massa nella cella %.1e masse solari) | collassata a z=3,2,1,0.5,0: %s | solo sorgente %+7.2f | sorgente+pozzo %+7.2f (bao %.1f sn %.1f cmb %.1f, H0 %.1f)'%(dc,mc,' '.join('%.3f'%np.interp(math.log(1/(1+z)),LA,F_) for z in (3,2,1,0.5,0)),r1.fun-L0,r2.fun-L0,cb,csn,ccm,100*hh),flush=True)
        if migliore is None or r2.fun<migliore['chi2']: migliore=dict(alpha=math.exp(r2.x[2]),gamma=math.exp(r2.x[3]),wm=float(r2.x[0]),wb=float(r2.x[1]),h=hh,chi2=float(r2.fun),delta=float(r2.fun-L0),dc=dc,pozzo=POZZO)
    OUT='parametri3d.json' if POZZO=='materia' else 'parametri3d_%s.json'%POZZO
    json.dump(migliore,open(OUT,'w'))
    print('migliore: soglia %.1f, chi2 da LCDM %+.2f, alpha %.4g gamma %.4g -> scritto %s'%(migliore['dc'],migliore['delta'],migliore['alpha'],migliore['gamma'],OUT))
