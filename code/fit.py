import numpy as np, sys, json
from scipy.optimize import minimize, brentq
from scipy.interpolate import interp1d
import os
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.join(HERE,'data')+os.sep
# ---------- data ----------
import csv
rows=[l.split() for l in open(D+'pp.dat')]
h=rows[0]; rows=rows[1:]
col=lambda n: np.array([float(r[h.index(n)]) for r in rows])
zhd=col('zHD'); mb=col('m_b_corr'); cal=col('IS_CALIBRATOR'); zhel=col('zHEL')
cov=np.loadtxt(D+'pp.cov',skiprows=1).reshape(1701,1701)
m=(zhd>0.01)   # selezione ufficiale Pantheon+ senza SH0ES: 1590 supernovae
zs=zhd[m]; zsh=zhel[m]; mbs=mb[m]; Ci=np.linalg.inv(cov[np.ix_(m,m)])
one=np.ones(m.sum()); Ci1=Ci@one; s11=one@Ci1
bz=[];bv=[];bt=[]
for l in open(D+'desi_mean.txt'):
    if l[0]=='#':continue
    a=l.split(); bz.append(float(a[0]));bv.append(float(a[1]));bt.append(a[2])
bz=np.array(bz);bv=np.array(bv);bCi=np.linalg.inv(np.loadtxt(D+'desi_cov.txt'))
Or=9e-5
zg=np.concatenate([[0],np.expm1(np.linspace(np.log(1.0005),np.log(3.6),700))])
def psi(z): return (1+z)**2.7/(1+((1+z)/2.9)**5.6)   # Madau-Dickinson 2014, normalizzata a meno di costante
# ---------- modelli: restituiscono funzione x(z)=rho_DE/rho_crit0 ----------
def x_lcdm(p): return lambda z: (1-p[0]-Or)+0*z
def x_w0wa(p):
    Om,w0,wa=p; return lambda z:(1-Om-Or)*(1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z))
def x_inst(p):  # B istantaneo: rho_DE proporzionale al tasso di dissipazione
    return lambda z:(1-p[0]-Or)*psi(z)/psi(0)
def make_ss(n,mm,zmax=30,N=500):
    u=np.linspace(np.log(1+zmax),0,N+1)  # u=ln(1+z) decrescente
    def solve(Om,beta,alpha):
        x=0.0; xs=[0.0]
        def f(u,x):
            z=np.expm1(u); E=np.sqrt(Om*(1+z)**3+Or*(1+z)**4+max(x,0))
            return -(alpha*psi(z)*(1+z)**n-beta*(1+z)**mm*x)/E   # dx/du = -dx/dlna
        for i in range(N):
            hh=u[i+1]-u[i]
            k1=f(u[i],x);k2=f(u[i]+hh/2,x+hh*k1/2);k3=f(u[i]+hh/2,x+hh*k2/2);k4=f(u[i+1],x+hh*k3)
            x+=hh*(k1+2*k2+2*k3+k4)/6; xs.append(x)
        return np.array(xs)
    def model(p):
        Om,beta=p; tgt=1-Om-Or
        g=lambda la: solve(Om,beta,np.exp(la))[-1]-tgt
        la=brentq(g,-40,15,xtol=1e-10)
        xs=solve(Om,beta,np.exp(la))
        it=interp1d(u[::-1],xs[::-1],kind='cubic')
        fn=lambda z: it(np.log1p(z))
        fn.alpha=np.exp(la); return fn
    return model
def chi2(xf,Om,q,parts=False):
    E=np.sqrt(Om*(1+zg)**3+Or*(1+zg)**4+xf(zg))
    if not np.all(np.isfinite(E)): return 1e9
    I=np.concatenate([[0],np.cumsum(np.diff(zg)*(1/E[1:]+1/E[:-1])/2)])
    DM=interp1d(zg,I,kind='cubic'); DH=interp1d(zg,1/E,kind='cubic')
    th=np.array([q*DM(z) if t.startswith('DM') else q*DH(z) if t.startswith('DH') else q*(z*DM(z)**2*DH(z))**(1/3) for z,t in zip(bz,bt)])
    r=bv-th; cb=r@bCi@r
    mu=5*np.log10((1+zsh)*DM(zs)); d=mbs-mu
    cs=d@Ci@d-(d@Ci1)**2/s11
    return (cb,cs) if parts else cb+cs
models={
 'LCDM':(x_lcdm,[0.3],[(0.05,0.6)]),
 'w0wa':(x_w0wa,[0.3,-0.9,-0.5],[(0.05,0.6),(-2,0),(-4,2)]),
 'B_istantaneo':(x_inst,[0.3],[(0.02,0.9)]),
}
for n in (0,3):
    for mm in (0,3):
        models[f'SS_n{n}_m{mm}']=(make_ss(n,mm),[0.3,1.0],[(0.05,0.6),(1e-4,30)])
out={}
for name,(mk,p0,bnds) in models.items():
    best=None
    starts=[p0] if not name.startswith('SS') else [[0.3,b] for b in (0.05,0.5,2.0,8.0)]
    for s in starts:
        def F(v):
            try: return chi2(mk(v[:-1]),v[0],v[-1])
            except Exception as e: return 1e9
        r=minimize(F,list(s)+[30.0],method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-3,maxiter=4000))
        if best is None or r.fun<best.fun: best=r
    v=best.x; xf=mk(v[:-1]); cb,cs=chi2(xf,v[0],v[-1],True)
    out[name]=dict(par=[float(a) for a in v],chi2=float(best.fun),bao=float(cb),sn=float(cs),k=len(v))
    if hasattr(xf,'alpha'): out[name]['alpha']=float(xf.alpha)
    zz=np.array([0,0.25,0.5,1,1.5,2,2.33]); out[name]['x_ratio']=[float(a) for a in xf(zz)/xf(0.0)]
    print(name,out[name],flush=True)
json.dump(out,open(os.path.join(HERE,'results.json'),'w'),indent=1)
