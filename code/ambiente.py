"""Test 2 (environment): does the expansion measured by supernovae depend on the galaxy density along the line of sight?
Map: 2M++ velocity field (Carrick et al. 2015) distributed by pvhub; the density is derived from the divergence."""
import numpy as np, math, os
from scipy.interpolate import RegularGridInterpolator as RGI, interp1d
from astropy.coordinates import SkyCoord
import astropy.units as u
HERE=os.path.dirname(os.path.abspath(__file__)); D=os.path.join(HERE,'data'); C=299792.458; BETA=0.43
NPZ=os.path.join(D,'twompp_vel.npz')
if not os.path.exists(NPZ):   # first run: converts the downloaded map (118 MB) into a compact grid
    import pandas as pd
    d=pd.read_csv(os.path.join(D,'2MPP_redshift.txt')); n=round(len(d)**(1/3))
    np.savez_compressed(NPZ,ax=np.unique(d.sgx.values),V=np.stack([d[c].values.reshape(n,n,n) for c in ('vpred_x','vpred_y','vpred_z')]).astype('f4'),P=d.vproj_2MPP.values.reshape(n,n,n).astype('f4'))
g=np.load(NPZ); ax=g['ax']; V=g['V']; P=g['P']
h=ax[1]-ax[0]
div=np.gradient(V[0],h,axis=0)+np.gradient(V[1],h,axis=1)+np.gradient(V[2],h,axis=2)
delta=-div/BETA; ok=np.isfinite(delta)
fd=RGI((ax,ax,ax),np.where(ok,delta,0.0),bounds_error=False,fill_value=0.0)
fm=RGI((ax,ax,ax),ok.astype(float),bounds_error=False,fill_value=0.0)
fp=RGI((ax,ax,ax),np.where(np.isfinite(P),P,0.0),bounds_error=False,fill_value=0.0)
def sgvec(l,b):
    l,b=np.radians(l),np.radians(b); return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)]).T
print('map check: density in known clusters (expected > 0) and in the Local Void (expected < 0)')
for nm,l,b,cz in [('Virgo',102.9,-2.3,1300),('Coma',89.6,8.1,7000),('Perseo',347.7,-14.1,5200),('Norma',188.5,7.3,4900),('Shapley',149.6,-1.2,14500),('Vuoto Locale',45,60,1500)]:
    print('  %-13s delta = %+.2f'%(nm,fd(sgvec(l,b)*cz)[0]))
print('  mean %.3f  std %.2f'%(delta[ok].mean(),delta[ok].std()))
rows=[l.split() for l in open(os.path.join(D,'pp.dat'))]; hd=rows[0]; rows=rows[1:]
col=lambda n: np.array([float(r[hd.index(n)]) for r in rows])
zhd,zcmb,zhel,mb,cal,ra,dec,vpec=[col(n) for n in ('zHD','zCMB','zHEL','m_b_corr','IS_CALIBRATOR','RA','DEC','VPEC')]
cov=np.loadtxt(os.path.join(D,'pp.cov'),skiprows=1).reshape(1701,1701)
sg=SkyCoord(ra=ra*u.deg,dec=dec*u.deg).supergalactic; nh=sgvec(sg.sgl.deg,sg.sgb.deg)
pos=nh*(C*zcmb)[:,None]
sel=(zhd>0.01)&(C*zcmb<19000)&(fm(pos)>0.99)
print('supernovae used:',sel.sum(),' redshift %.3f-%.3f'%(zhd[sel].min(),zhd[sel].max()))
t=np.linspace(0.02,1,60)
dlos=np.array([fd(p[None,:]*t[:,None]).mean() for p in pos]); dloc=fd(pos); vmap=fp(pos)
zg=np.concatenate([np.linspace(0,0.1,500),np.linspace(0.1,2.5,2000)[1:]]); E=np.sqrt(0.334*(1+zg)**3+0.666); I=np.concatenate([[0],np.cumsum(np.diff(zg)*(1/E[1:]+1/E[:-1])/2)]); DM=interp1d(zg,I,kind='cubic')
mu=lambda z: 5*np.log10((1+zhel)*DM(z))
Ci=np.linalg.inv(cov[np.ix_(sel,sel)]); K=math.log(10)/5
def gls(x,d):
    A=np.c_[np.ones(sel.sum()),x[sel]]; Cv=np.linalg.inv(A.T@Ci@A); b=Cv@(A.T@Ci@d[sel]); return b[1],math.sqrt(Cv[1,1])
rA=mb-mu(zcmb); rB=mb-mu(zhd)
print('\nA) check: UNcorrected residuals against the gravitational prediction (expected slope 1)')
for lab,v in (('velocita\' dal file Pantheon+',vpec),('velocita\' dalla mappa',vmap)):
    s,e=gls(-(5/math.log(10))*v/(C*zcmb),rA); print('   %-28s slope = %.2f +- %.2f'%(lab,s,e))
print('B) gravitational effect as a function of the mean density along the line of sight (uncorrected residuals)')
s,e=gls(dlos,rA); bg=-K*s; print('   d lnH / d delta = %+.3f +- %.3f'%(bg,K*e))
print('C) TEST: residual effect after removing gravity (corrected residuals, zHD)')
s,e=gls(dlos,rB); print('   density along the line of sight: d lnH / d delta = %+.3f +- %.3f'%(-K*s,K*e))
s2,e2=gls(dloc,rB); print('   density at the supernova position:  d lnH / d delta = %+.3f +- %.3f'%(-K*s2,K*e2))
print('   std of the density along the line of sight in the sample: %.2f'%dlos[sel].std())
# voids against dense regions
med=np.median(dlos[sel])
for lab,m in (('less dense lines of sight',dlos<med),('denser lines of sight',dlos>=med)):
    mm=m[sel].astype(float); 
A=np.c_[np.ones(sel.sum()),(dlos<med)[sel].astype(float)]; Cv=np.linalg.inv(A.T@Ci@A); b=Cv@(A.T@Ci@rB[sel])
print('   less dense half minus denser half: expansion difference = %+.2f +- %.2f %%'%(-K*b[1]*100,K*math.sqrt(Cv[1,1])*100))
