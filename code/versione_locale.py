"""Local version: (A) local absorption of space, (B) local return of energy to matter."""
import numpy as np, math, json
from scipy.interpolate import CubicSpline
from check_growth import *
C=json.load(open('cmb2_cum_n3_U.json')); R=json.load(open('results3.json'))
p=C['p'][:4]; al=math.exp(p[2]); ga=math.exp(p[3]); E,o=past_interacting(p[0],al,ga,3,'cum'); m,x,M=o.T; hh=float(E[-1])
c=Delta_curves(p,'cum'); pl=R['LCDM']['p']; El=e_fluid(pl[0],pl[2]); dl=solveD(El,pl[0]*(1+z)**3,0*z)
lna=np.linspace(N[0],0,3000); zz=np.expm1(-lna); a=np.exp(lna)
# ---------- (A) local absorption: contrast of the "extra space" and the resulting force
# local sink: dV/dt = -mean_sink*(1+delta)  ->  d(deltaV)/dt = -mean_sink*delta ; contrast kappa = -deltaV/(rho_m*delta)
sE=CubicSpline(N,E); sm=CubicSpline(N,m); sM=CubicSpline(N,M); sx=CubicSpline(N,x)
dU=c['dU'](lna)[0]                       # matter density contrast (growth)
sink=ga*sM(lna); integrand=sink*dU/sE(lna)      # dt = dlna/E
dV=-np.concatenate([[0],np.cumsum((integrand[1:]+integrand[:-1])/2*np.diff(lna))])   # deltaV per unit initial delta
kappa=-dV/(sm(lna)*dU)
print('(A) contrast of the extra space per unit matter contrast (kappa), if the sink acts where matter is:')
for zq in (2,1,0.5,0):
    i=np.argmin(abs(zz-zq)); print('    z=%.1f  kappa=%.3f   (extra space / matter = %.2f)'%(zq,kappa[i],sx(lna[i])/sm(lna[i])))
k0=kappa[-1]; Om=m[-1]/hh**2; cH=2997.92/hh
lamJ=2*math.pi*cH*math.sqrt(k0/(1.5*Om))
print('    length below which the gradient force prevents structures from growing, today: %.0f Mpc (Hubble radius %.0f Mpc)'%(lamJ,cH))
for lam in (1,10,50): print('    for structures of %d Mpc to grow, kappa < %.1e is needed  (the local model gives %.2f: %.0e times too much)'%(lam,1.5*Om*(lam/(2*math.pi*cH))**2,k0,k0/(1.5*Om*(lam/(2*math.pi*cH))**2)))
print('    equivalent speed of the "pressure": %.2f c'%math.sqrt(k0))
# ---------- (B) local return
print('(B1) local return to all matter, ordinary matter included: mass growth today')
g0=float(c['g'](0.0))*hh/9.7779e9     # per year
print('    dM/dt / M = %.1e per year ; Solar System limit on d(GM_Sun)/dt: about 1e-13 per year -> %.0f times beyond'%(g0,g0/1e-13))
print('(B2) local return to dark matter only: integrated Sachs-Wolfe effect (growing or decaying potentials)')
def phi(sol,Mc): return Mc*sol/a
PhiL=phi(c['dL'](lna)[0],c['Mc'](lna)); PhiU=phi(c['dU'](lna)[0],c['Mc'](lna)); PhiC=phi(dl(lna)[0],pl[0])
dL_=c['dL'](lna)[0]; dU_=c['dU'](lna)[0]; dC_=dl(lna)[0]
dPhi=lambda P: np.gradient(P,lna)       # d Phi / d ln a  (negative = decaying potential = normal ISW)
sig=lambda P,d: dPhi(P)*d               # ISW-galaxy correlation signal per unit redshift (up to the bias)
print('    potential change between z and today:   z     LCDM     diffuse    local')
for zq in (1,0.5,0.25,0.1):
    i=np.argmin(abs(zz-zq)); print('                                           %.2f   %+.1f%%    %+.1f%%     %+.1f%%'%(zq,100*(PhiC[-1]/PhiC[i]-1),100*(PhiU[-1]/PhiU[i]-1),100*(PhiL[-1]/PhiL[i]-1)))
print('    ISW-galaxy amplitude relative to LCDM (1 = like LCDM; negative = opposite sign):')
for lab,z0,s in (('galassie vicine (z~0,1)',0.1,0.05),('galassie a z~0,3',0.3,0.1),('galassie a z~0,5',0.5,0.15),('radiosorgenti tipo NVSS (z~0,9, larga)',0.9,0.5)):
    w=np.exp(-0.5*((zz-z0)/s)**2)*a      # n(z) dz = n(z) (dz/dlna) dlna, dz/dlna = -(1+z) -> weight ~ n(z)/a ... here weight per dlna
    w=np.exp(-0.5*((zz-z0)/s)**2)/a
    A=lambda P,d: np.sum(w*sig(P,d))
    print('      %-40s diffuse %+.2f   local %+.2f'%(lab,A(PhiU,dU_)/A(PhiC,dC_),A(PhiL,dL_)/A(PhiC,dC_)))
print('    mass growth in structures (local version): from z=0.5 to today %+.0f%%, from z=0.25 to today %+.0f%%'%(100*(float(c['Mc'](0))/float(c['Mc'](-math.log(1.5)))-1),100*(float(c['Mc'](0))/float(c['Mc'](-math.log(1.25)))-1)))
