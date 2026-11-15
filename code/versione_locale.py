"""Versione locale: (A) assorbimento locale dello spazio, (B) restituzione locale dell'energia alla materia."""
import numpy as np, math, json
from scipy.interpolate import CubicSpline
from check_growth import *
C=json.load(open('cmb2_cum_n3_U.json')); R=json.load(open('results3.json'))
p=C['p'][:4]; al=math.exp(p[2]); ga=math.exp(p[3]); E,o=past_interacting(p[0],al,ga,3,'cum'); m,x,M=o.T; hh=float(E[-1])
c=Delta_curves(p,'cum'); pl=R['LCDM']['p']; El=e_fluid(pl[0],pl[2]); dl=solveD(El,pl[0]*(1+z)**3,0*z)
lna=np.linspace(N[0],0,3000); zz=np.expm1(-lna); a=np.exp(lna)
# ---------- (A) assorbimento locale: contrasto dello "spazio in piu'" e forza che ne deriva
# pozzo locale: dV/dt = -pozzo_medio*(1+delta)  ->  d(deltaV)/dt = -pozzo_medio*delta ; contrasto kappa = -deltaV/(rho_m*delta)
sE=CubicSpline(N,E); sm=CubicSpline(N,m); sM=CubicSpline(N,M); sx=CubicSpline(N,x)
dU=c['dU'](lna)[0]                       # contrasto di densita' della materia (crescita)
sink=ga*sM(lna); integrand=sink*dU/sE(lna)      # dt = dlna/E
dV=-np.concatenate([[0],np.cumsum((integrand[1:]+integrand[:-1])/2*np.diff(lna))])   # deltaV per unita' di delta iniziale
kappa=-dV/(sm(lna)*dU)
print('(A) contrasto dello spazio in piu\' per unita\' di contrasto di materia (kappa), se il pozzo agisce dove sta la materia:')
for zq in (2,1,0.5,0):
    i=np.argmin(abs(zz-zq)); print('    z=%.1f  kappa=%.3f   (spazio in piu\' / materia = %.2f)'%(zq,kappa[i],sx(lna[i])/sm(lna[i])))
k0=kappa[-1]; Om=m[-1]/hh**2; cH=2997.92/hh
lamJ=2*math.pi*cH*math.sqrt(k0/(1.5*Om))
print('    lunghezza sotto cui la forza del gradiente impedisce alle strutture di crescere, oggi: %.0f Mpc (raggio di Hubble %.0f Mpc)'%(lamJ,cH))
for lam in (1,10,50): print('    perche\' strutture di %d Mpc possano crescere serve kappa < %.1e  (il modello locale da\' %.2f: %.0e volte troppo)'%(lam,1.5*Om*(lam/(2*math.pi*cH))**2,k0,k0/(1.5*Om*(lam/(2*math.pi*cH))**2)))
print('    velocita\' equivalente della "pressione": %.2f c'%math.sqrt(k0))
# ---------- (B) restituzione locale
print('(B1) restituzione locale a tutta la materia, ordinaria compresa: crescita della massa oggi')
g0=float(c['g'](0.0))*hh/9.7779e9     # per anno
print('    dM/dt / M = %.1e per anno ; limite dal sistema solare su d(GM_Sole)/dt: circa 1e-13 per anno -> %.0f volte oltre'%(g0,g0/1e-13))
print('(B2) restituzione locale alla sola materia oscura: effetto Sachs-Wolfe integrato (potenziali che crescono o calano)')
def phi(sol,Mc): return Mc*sol/a
PhiL=phi(c['dL'](lna)[0],c['Mc'](lna)); PhiU=phi(c['dU'](lna)[0],c['Mc'](lna)); PhiC=phi(dl(lna)[0],pl[0])
dL_=c['dL'](lna)[0]; dU_=c['dU'](lna)[0]; dC_=dl(lna)[0]
dPhi=lambda P: np.gradient(P,lna)       # d Phi / d ln a  (negativo = potenziale che decade = ISW normale)
sig=lambda P,d: dPhi(P)*d               # segnale di correlazione ISW-galassie per unita' di redshift (a meno del bias)
print('    variazione del potenziale tra z e oggi:   z     LCDM     diffusa    locale')
for zq in (1,0.5,0.25,0.1):
    i=np.argmin(abs(zz-zq)); print('                                           %.2f   %+.1f%%    %+.1f%%     %+.1f%%'%(zq,100*(PhiC[-1]/PhiC[i]-1),100*(PhiU[-1]/PhiU[i]-1),100*(PhiL[-1]/PhiL[i]-1)))
print('    ampiezza ISW-galassie rispetto a LCDM (1 = come LCDM; negativo = segno opposto):')
for lab,z0,s in (('galassie vicine (z~0,1)',0.1,0.05),('galassie a z~0,3',0.3,0.1),('galassie a z~0,5',0.5,0.15),('radiosorgenti tipo NVSS (z~0,9, larga)',0.9,0.5)):
    w=np.exp(-0.5*((zz-z0)/s)**2)*a      # n(z) dz = n(z) (dz/dlna) dlna, dz/dlna = -(1+z) -> peso ~ n(z)/a ... qui peso per dlna
    w=np.exp(-0.5*((zz-z0)/s)**2)/a
    A=lambda P,d: np.sum(w*sig(P,d))
    print('      %-40s diffusa %+.2f   locale %+.2f'%(lab,A(PhiU,dU_)/A(PhiC,dC_),A(PhiL,dL_)/A(PhiC,dC_)))
print('    crescita della massa nelle strutture (versione locale): da z=0,5 a oggi %+.0f%%, da z=0,25 a oggi %+.0f%%'%(100*(float(c['Mc'](0))/float(c['Mc'](-math.log(1.5)))-1),100*(float(c['Mc'](0))/float(c['Mc'](-math.log(1.25)))-1)))
