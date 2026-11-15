"""Modello sorgente-pozzo con bilancio dell'energia coerente (vuoto interagente),
confronto con Pantheon+ (SN), DESI DR2 (BAO) e prior compressi Planck 2018 (CMB).
Unita': densita' fisiche omega = Omega*h^2, H in unita' di 100 km/s/Mpc, tempo in 9.778 Gyr."""
import numpy as np, math, json, sys
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp, quad
import os
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.join(HERE,'data')+os.sep; CK=2997.92458; TU=9.7779
import os as _os
SNSET=_os.environ.get('SNSET','pantheon')     # pantheon | union3 | desy5 | dovekie
def _load_sn(name):
    """Ritorna z (cosmologico), z eliocentrico, magnitudini o moduli di distanza, inversa della covarianza."""
    S=D+'sn/'
    if name=='pantheon':
        rows=[l.split() for l in open(D+'pp.dat')]; h_=rows[0]; rows=rows[1:]
        col=lambda n: np.array([float(r[h_.index(n)]) for r in rows])
        zhd=col('zHD'); mb=col('m_b_corr'); zhel=col('zHEL')
        cov=np.loadtxt(D+'pp.cov',skiprows=1).reshape(1701,1701)
        msk=(zhd>0.01)   # selezione ufficiale Pantheon+ senza SH0ES: 1590 supernovae
        return zhd[msk],zhel[msk],mb[msk],np.linalg.inv(cov[np.ix_(msk,msk)])
    if name=='union3':     # 22 intervalli di redshift (Rubin et al. 2023)
        d=np.loadtxt(S+'Union3/lcparam_full.txt',usecols=(1,2,4)); c=np.loadtxt(S+'Union3/mag_covmat.txt',skiprows=1).reshape(22,22)
        return d[:,0],d[:,1],d[:,2],np.linalg.inv(c)
    if name=='desy5':      # DES-SN5YR, 1829 supernovae
        d=np.genfromtxt(S+'DESY5/DES-SN5YR_HD.csv',delimiter=',',names=True); f=S+'DESY5/covsys.npy'
        if not _os.path.exists(f): np.save(f,np.loadtxt(S+'DESY5/covsys_000.txt',skiprows=1).reshape(len(d),len(d)))
        c=np.load(f)+np.diag(d['MUERR_FINAL']**2)
        return d['zHD'],d['zHEL'],d['MU'],np.linalg.inv(c)
    if name=='dovekie':    # DES ricalibrato (DES-Dovekie), 1820 supernovae; il file contiene gia' l'inversa della covarianza totale
        d=np.genfromtxt(S+'DES-Dovekie/DES-Dovekie_HD.csv',delimiter=',',names=True); n=len(d)
        p=np.load(S+'DES-Dovekie/covtot_inv_000.npz')['cov']; ci=np.zeros((n,n)); ci[np.triu_indices(n)]=p; ci=ci+ci.T-np.diag(np.diag(ci))
        return d['zHD'],d['zHEL'],d['MU'],ci
    raise ValueError(name)
zs,zsh,mbs,Ci=_load_sn(SNSET)
one=np.ones(len(zs)); Ci1=Ci@one; s11=one@Ci1
bz=[];bv=[];bt=[]
for l in open(D+'desi_mean.txt'):
    if l[0]=='#':continue
    a=l.split(); bz.append(float(a[0]));bv.append(float(a[1]));bt.append(a[2])
bz=np.array(bz);bv=np.array(bv);bCi=np.linalg.inv(np.loadtxt(D+'desi_cov.txt'))
# Prior di distanza Planck 2018 (Chen, Huang & Wang 2019): R, l_a, omega_b
cm=np.array([1.7502,301.471,0.02236]); cs_=np.array([0.0046,0.090,0.00015])
cc=np.array([[1,0.46,-0.66],[0.46,1,-0.33],[-0.66,-0.33,1]]); cCi=np.linalg.inv(cc*np.outer(cs_,cs_))
OR=4.184e-5; ZMAX=30.0; NU=300
U=np.linspace(math.log(1+ZMAX),0,NU+1); UL=U.tolist()
def psi(z): return (1+z)**2.7/(1+((1+z)/2.9)**5.6)
PS0=psi(0.0)
def past_interacting(wm,al,ga,n,kind):
    m=wm*(1+ZMAX)**3; x=0.0; M=0.0; out=[(m,x,M)]
    def rhs(u,m,x,M):
        z=math.expm1(u); o=1+z; E2=m+OR*o**4+x
        if E2<=0 or m<=0: raise ValueError
        E=math.sqrt(E2); ps=o**2.7/(1+(o/2.9)**5.6)
        sink=ga*x if kind=='prop' else ga*M if kind=='cum' else ga
        Q=al*ps*o**n-sink
        return 3*m+Q/E, -Q/E, -ps/E
    for i in range(NU):
        hh=UL[i+1]-UL[i]; u=UL[i]
        a1,b1,c1=rhs(u,m,x,M); a2,b2,c2=rhs(u+hh/2,m+hh*a1/2,x+hh*b1/2,M+hh*c1/2)
        a3,b3,c3=rhs(u+hh/2,m+hh*a2/2,x+hh*b2/2,M+hh*c2/2); a4,b4,c4=rhs(u+hh,m+hh*a3,x+hh*b3,M+hh*c3)
        m+=hh*(a1+2*a2+2*a3+a4)/6; x+=hh*(b1+2*b2+2*b3+b4)/6; M+=hh*(c1+2*c2+2*c3+c4)/6; out.append((m,x,M))
    o=np.array(out); z=np.expm1(U)
    return np.sqrt(o[:,0]+OR*(1+z)**4+o[:,1]), o
def observables(e_u,wm,wb,rd=None):
    """e_u = H/100 sulla griglia U (z da 30 a 0)."""
    eu=interp1d(U[::-1],e_u[::-1],kind='cubic'); hh=float(e_u[-1])
    zg=np.concatenate([[0],np.expm1(np.linspace(math.log(1.0005),math.log(1+ZMAX),900))])
    ie=1/eu(np.log1p(zg)); I=np.concatenate([[0],np.cumsum(np.diff(zg)*(ie[1:]+ie[:-1])/2)])
    DM=interp1d(zg,I,kind='cubic')
    g1=0.0783*wb**-0.238/(1+39.5*wb**0.763); g2=0.560/(1+21.1*wb**1.81)
    zst=1048*(1+0.00124*wb**-0.738)*(1+g1*wm**g2)
    eh=lambda z: math.sqrt(wm*(1+z)**3+OR*(1+z)**4)
    Ist=I[-1]+quad(lambda z:1/eh(z),ZMAX,zst)[0]
    Rb=lambda z: 31500*wb*(2.7255/2.7)**-4/(1+z)
    rs=quad(lambda a: 1/(a*a*eh(1/a-1)*math.sqrt(3*(1+Rb(1/a-1)))),0,1/(1+zst),epsabs=1e-12)[0]
    cmb=np.array([math.sqrt(wm)*Ist, math.pi*Ist/rs, wb])
    if rd is None: rd=147.05*(wb/0.02236)**-0.13*(wm/0.1432)**-0.23   # formula DESI; con omega_m totale riproduce CAMB entro 0,05%
    q=CK/rd
    th=np.array([q*DM(z) if t.startswith('DM') else q*ie_at(eu,z) if t.startswith('DH') else q*(z*DM(z)**2*ie_at(eu,z))**(1/3) for z,t in zip(bz,bt)])
    r=bv-th; cb=float(r@bCi@r)
    d=mbs-5*np.log10((1+zsh)*DM(zs)); csn=float(d@Ci@d-(d@Ci1)**2/s11)
    rc=cmb-cm; ccm=float(rc@cCi@rc)
    return cb,csn,ccm,hh
def ie_at(eu,z): return 1/float(eu(math.log1p(z)))
def e_fluid(wm,wx,w0=-1.0,wa=0.0):
    z=np.expm1(U); return np.sqrt(wm*(1+z)**3+OR*(1+z)**4+wx*(1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z)))
# ---- definizione modelli: p -> (e_u, wm_early, wb, extra)
def M_lcdm(p): return e_fluid(p[0],p[2]),p[0],p[1],None
def M_w0wa(p): return e_fluid(p[0],p[2],p[3],p[4]),p[0],p[1],None
def M_int(n,kind):
    def f(p):
        e,o=past_interacting(p[0],math.exp(p[2]),math.exp(p[3]) if len(p)>3 else 0.0,n,kind); return e,p[0],p[1],o
    return f
MODELS={
 'LCDM':(M_lcdm,[[0.143,0.02236,0.32]]),
 'w0wa':(M_w0wa,[[0.143,0.02236,0.32,-0.85,-0.5],[0.143,0.02236,0.32,-0.7,-1.0]]),
 'solo_rubinetto_n3':(M_int(3,'cum'),[[0.143,0.02236,math.log(0.002)]]),
 'prop_n3':(M_int(3,'prop'),[[0.143,0.02236,math.log(0.0023),math.log(0.17)],[0.143,0.02236,math.log(0.004),math.log(0.5)]]),
 'cum_n3':(M_int(3,'cum'),[[0.143,0.02236,math.log(0.0023),math.log(0.0107)],[0.143,0.02236,math.log(0.003),math.log(0.03)]]),
 'const_n3':(M_int(3,'const'),[[0.143,0.02236,math.log(0.0023),math.log(0.055)],[0.143,0.02236,math.log(0.003),math.log(0.15)]]),
 'cum_n0':(M_int(0,'cum'),[[0.143,0.02236,math.log(0.09),math.log(0.03)],[0.143,0.02236,math.log(0.2),math.log(0.1)]]),
}
def total(mk,p,parts=False):
    try:
        e,wm,wb,o=mk(p)
        if not np.all(np.isfinite(e)) or wm<=0.05 or wm>0.3 or wb<0.015 or wb>0.03: return 1e9
        cb,csn,ccm,hh=observables(e,wm,wb)
    except (ValueError,FloatingPointError,OverflowError): return 1e9
    return (cb,csn,ccm,hh,o) if parts else cb+csn+ccm
def fitmodel(mk,starts,fixed=None):
    best=None
    for s in starts:
        s=list(s)
        F=(lambda v: total(mk,v)) if fixed is None else (lambda v: total(mk,list(v)+[fixed]))
        r=minimize(F,s if fixed is None else s[:-1],method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-3,maxiter=5000,adaptive=True))
        r=minimize(F,r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-4,maxiter=5000,adaptive=True))
        if best is None or r.fun<best.fun: best=r
    return best
def future(p,n,kind,o,hh,tmax=60,efold_gyr=None):
    """Integra avanti nel tempo attraverso l'inversione. Attivita' stellare futura: decadimento esponenziale con il ritmo attuale."""
    al=math.exp(p[2]); ga=math.exp(p[3]); m0,x0,M0=o[-1]; k=2.7*hh if efold_gyr is None else TU/efold_gyr
    def f(t,y):
        a,v,m,x,M=y; ps=PS0*math.exp(-k*t)
        sink=ga*x if kind=='prop' else ga*M if kind=='cum' else ga
        Q=al*ps*a**-n-sink
        return [v, a*(-m/2-OR/a**4+x), -3*(v/a)*m-Q, Q, ps]
    turn=lambda t,y:y[1]; turn.direction=-1
    crunch=lambda t,y:y[0]-0.02; crunch.terminal=True; crunch.direction=-1
    xz=lambda t,y:y[3]; xz.direction=-1
    s=solve_ivp(f,[0,tmax],[1.0,hh,m0,x0,M0],events=[turn,crunch,xz],rtol=1e-9,atol=1e-12,max_step=0.02,dense_output=True)
    g=lambda e: float(e[0])*TU if len(e) else None
    tt=g(s.t_events[0]); res=dict(t_inversione_Gyr=tt,t_crunch_Gyr=g(s.t_events[1]),t_energia_oscura_zero_Gyr=g(s.t_events[2]))
    if tt: res['a_max']=float(s.sol(tt/TU)[0])
    ts=np.linspace(0,s.t[-1],25); res['curva']=[(round(float(t*TU),1),round(float(s.sol(t)[0]),3)) for t in ts]
    return res
if __name__=='__main__':
    R={}
    for name,(mk,starts) in MODELS.items():
        b=fitmodel(mk,starts); cb,csn,ccm,hh,o=total(mk,b.x,True)
        R[name]=dict(p=[float(a) for a in b.x],chi2=float(b.fun),bao=cb,sn=csn,cmb=ccm,h=hh,npar=len(b.x))
        if o is not None:
            R[name].update(Om0=float(o[-1,0]/hh**2),Ox0=float(o[-1,1]/hh**2),x_min=float(o[:,1].min()/hh**2),
                           materia_convertita=float(1-o[-1,0]/b.x[0]))
            zz=[0,0.25,0.5,1,1.5,2,2.33,5]; xi=interp1d(U[::-1],o[::-1,1]); R[name]['x_ratio']=[round(float(xi(math.log1p(z))/o[-1,1]),3) for z in zz]
            if name in('cum_n3','const_n3','cum_n0','prop_n3'):
                kind=name.split('_')[0]; n=int(name[-1]); R[name]['futuro']=future(b.x,n,kind,o,hh)
        print(name,json.dumps(R[name]),flush=True)
    json.dump(R,open('results3.json' if SNSET=='pantheon' else 'results3_%s.json'%SNSET,'w'),indent=1)
