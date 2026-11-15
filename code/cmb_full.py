"""Test 3 (versione 2): spettro completo del fondo cosmico (Planck 2018 plik-lite TTTEEE + TT a bassi multipoli) con CAMB.

Il modello sorgente-pozzo entra come fluido efficace: stessa espansione del modello, energia oscura efficace
rho_eff = rho_m + rho_x - omega_m,iniziale*(1+z)^3 con equazione di stato tabulata w(a).

Novita' rispetto alla versione 1:
  - l'orizzonte sonoro r_d per DESI e' quello calcolato da CAMB (prima: formula approssimata, +0,1%);
  - correzione dell'effetto di lente (pert='U'): nel fluido efficace la materia scambiata resta liscia; nel modello
    con scambio diffuso invece si muove con la materia circostante. La differenza nel potenziale gravitazionale
    (1-2%) e' calcolata dalle equazioni di crescita del modello e applicata come ampiezza di lente A_lens;
  - selezione ufficiale delle supernovae Pantheon+ (1590), tramite fit3.py.
pert='E' riproduce la versione 1 (nessuna correzione di lente)."""
import numpy as np, math, json, sys, os, time
os.environ.setdefault('OMP_NUM_THREADS','1')
import camb
from camb.dark_energy import DarkEnergyPPF
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from fit3 import *
HERE=os.path.dirname(os.path.abspath(__file__))
for d in (os.path.join(HERE,'planck-lite-py'),os.path.join(HERE,'..','ext','planck-lite-py')):
    if os.path.isdir(d): sys.path.insert(0,d); PLDIR=os.path.join(d,'data'); break
from planck_lite_py import PlanckLitePy
PL=PlanckLitePy(data_directory=PLDIR,year=2018,spectra='TTTEEE',use_low_ell_bins=True)
WNU=0.000645; ZT=20.0; NLN=-U; ZU=np.expm1(U)
INFO={'cum_n3':(3,'cum'),'const_n3':(3,'const'),'prop_n3':(3,'prop'),'solo_rubinetto_n3':(3,'cum')}   # (diluizione sorgente, tipo di pozzo)
# ---------- fluido efficace per CAMB
def w_table(e_u,wm):
    rho=e_u**2-wm*(1+ZU)**3-OR*(1+ZU)**4; k=(ZU<=ZT)
    if np.any(rho[k]<=0): raise ValueError('rho_eff negativa')
    lna=-U[k]; s=CubicSpline(lna,np.log(rho[k])); return np.exp(lna),-1-s(lna,1)/3
def spectra(e_u,wm,wb,logA,ns,tau,alens=1.0):
    hh=float(e_u[-1]); p=camb.CAMBparams(); a,w=w_table(e_u,wm); de=DarkEnergyPPF(); de.set_w_a_table(a,w); p.DarkEnergy=de
    p.set_cosmology(H0=100*hh,ombh2=wb,omch2=wm-wb-WNU,mnu=0.06,omk=0,tau=tau,Alens=alens)
    p.InitPower.set_params(As=math.exp(logA)*1e-10,ns=ns); p.set_for_lmax(2600,lens_potential_accuracy=1)
    r=camb.get_results(p); cl=r.get_cmb_power_spectra(p,CMB_unit='muK',spectra=['total'])['total']
    return cl,r
# ---------- crescita lineare: sorgente del potenziale Delta = (materia comovente) x (contrasto di densita')
def _grow(E,m,g):
    sE=CubicSpline(NLN,E); sm=CubicSpline(NLN,m); sg=CubicSpline(NLN,g)
    def f(n,y):
        e=sE(n); de=sE(n,1); gg=sg(n); dg=sg(n,1)
        return [y[1], -(de/e+2+gg/e)*y[1]+(1.5*sm(n)-2*e*gg-e*dg)/e**2*y[0]]
    a0=math.exp(NLN[0]); return solve_ivp(f,[NLN[0],0],[a0,a0],dense_output=True,rtol=1e-8,atol=1e-12).sol
def delta_ratio(name,pm,e_u,o,pert):
    """Rapporto Delta_modello/Delta_fluido-efficace in funzione di ln a (spline)."""
    n,kind=INFO[name]; al=math.exp(pm[2]); ga=math.exp(pm[3]) if len(pm)>3 else 0.0; m,x,M=o.T
    sink=ga*x if kind=='prop' else ga*M if kind=='cum' else ga*np.ones_like(x)
    g=-(al*psi(ZU)*(1+ZU)**n-sink)/m                     # ritmo con cui la materia guadagna energia, per unita' di materia
    dM=_grow(e_u,m,g if pert=='U' else 0*g); dE=_grow(e_u,pm[0]*(1+ZU)**3,0*g)
    lna=np.linspace(NLN[0],0,400); Mc=CubicSpline(NLN,m*(1+ZU)**-3)(lna)
    return CubicSpline(lna,Mc*dM(lna)[0]/(pm[0]*dE(lna)[0]))
# ---------- pesi in redshift della lente del fondo cosmico (Limber, modello LCDM di riferimento)
def _kernel():
    f=os.path.join(HERE,'nucleo_lente.npz')
    if os.path.exists(f):
        d=np.load(f); return d['z'],d['w']
    p=camb.set_params(H0=67.7,ombh2=0.0224,omch2=0.1195,mnu=0.06,As=2.1e-9,ns=0.966,tau=0.055)
    PK=camb.get_matter_power_interpolator(p,zmax=ZMAX,kmax=5.0,nonlinear=False,hubble_units=False,k_hunit=False,zmin=0)
    bg=camb.get_background(p); z=np.linspace(0.01,ZMAX-0.01,1500); chi=bg.comoving_radial_distance(z); chis=bg.comoving_radial_distance(1089.0)
    dchi=np.gradient(chi,z); w=np.zeros_like(z)
    for L in (40,100,200,400):
        wl=dchi*((chis-chi)/chis)**2*(1+z)**2*np.array([PK.P(zz,min((L+0.5)/c,4.9)) for zz,c in zip(z,chi)]); w+=wl/wl.sum()
    w/=w.sum()
    try: np.savez(f+'.%d.tmp.npz'%os.getpid(),z=z,w=w); os.replace(f+'.%d.tmp.npz'%os.getpid(),f)
    except OSError: pass
    return z,w
KZ,KW=_kernel()
def alens_of(name,pm,e_u,o,pert):
    if pert=='E' or name not in INFO: return 1.0
    R=delta_ratio(name,pm,e_u,o,pert)(-np.log1p(KZ)); return float(np.sum(KW*R**2))
# ---------- chi quadro totale
def chi2_all(name,p,parts=False,pert='U'):
    """p = parametri del modello (come in fit3) + [ln(1e10 A_s), n_s, tau]."""
    pm=list(p[:-3]); logA,ns,tau=p[-3:]; mk=MODELS[name][0]
    try:
        e,wm,wb,o=mk(pm)
        if not np.all(np.isfinite(e)) or not(0.12<wm<0.17) or not(0.020<wb<0.025) or not(0.9<ns<1.05) or not(0.01<tau<0.12): return 1e9
        al=alens_of(name,pm,e,o,pert)
        cl,r=spectra(e,wm,wb,logA,ns,tau,al)
        rd=float(r.get_derived_params()['rdrag'])
        cb,csn,_,hh=observables(e,wm,wb,rd=rd)
        cpl=-2*PL.loglike(cl[2:,0],cl[2:,3],cl[2:,1],2)
    except Exception as ex:
        return 1e9
    ctau=((tau-0.0544)/0.0073)**2
    if parts: return dict(planck=float(cpl),tau=float(ctau),bao=float(cb),sn=float(csn),tot=float(cpl+ctau+cb+csn),h=hh,alens=al,rd=rd)
    return cpl+ctau+cb+csn
START={'LCDM':[0.140929,0.0225202,0.325992,3.05048,0.970577,0.0606715],
       'w0wa':[0.142007,0.0224475,0.315483,-0.847838,-0.555816,3.04988,0.967962,0.0589784],
       'cum_n3':[0.142507,0.0223998,-6.05012,-3.54678,3.05094,0.966548,0.0589106],
       'const_n3':[0.1426,0.0224,-5.93,-2.05,3.049,0.9663,0.0577],
       'prop_n3':[0.142287,0.0224111,-5.9587,-0.922367,3.04946,0.967103,0.0584639],
       'solo_rubinetto_n3':[0.1440,0.02235,-6.18,3.049,0.9650,0.0570]}
SCALE={'LCDM':[3e-4,6e-5,1.5e-3,0.004,0.002,0.003],'w0wa':[3e-4,6e-5,2e-3,0.03,0.12,0.004,0.002,0.003],
       'cum_n3':[3e-4,6e-5,0.02,0.05,0.004,0.002,0.003],'const_n3':[3e-4,6e-5,0.02,0.05,0.004,0.002,0.003],
       'prop_n3':[3e-4,6e-5,0.02,0.05,0.004,0.002,0.003],'solo_rubinetto_n3':[3e-4,6e-5,0.02,0.004,0.002,0.003]}
if __name__=='__main__':
    name=sys.argv[1]; pert=sys.argv[2] if len(sys.argv)>2 else 'U'; nit=int(sys.argv[3]) if len(sys.argv)>3 else 400
    x0=np.array(START[name]); sc=np.array(SCALE[name]); N_=len(x0); t0=time.time(); best=[1e9,None]; n=[0]
    suf='' if SNSET=='pantheon' else '_'+SNSET
    if suf:   # con un altro campione di supernovae parte dal fit rapido corrispondente
        r3=json.load(open(os.path.join(HERE,'results3%s.json'%suf)))[name]['p']; x0[:len(r3)]=r3
    def F(u):
        v=x0+u*sc; c=chi2_all(name,list(v),pert=pert); n[0]+=1
        if c<best[0]:
            best[0]=c; best[1]=list(v); print('n=%d t=%.0fs chi2=%.2f p=%s'%(n[0],time.time()-t0,c,np.round(v,5).tolist()),flush=True)
        return c
    minimize(F,np.zeros(N_),method='Nelder-Mead',options=dict(initial_simplex=np.vstack([np.zeros(N_),np.eye(N_)]),xatol=0.03,fatol=0.01,maxiter=nit,adaptive=True))
    out=chi2_all(name,best[1],True,pert=pert); out['p']=[float(a) for a in best[1]]; out['sec']=time.time()-t0; out['pert']=pert
    print('FINALE',name,pert,json.dumps(out),flush=True); json.dump(out,open(os.path.join(HERE,f'cmb2_{name}_{pert}{suf}.json'),'w'))
