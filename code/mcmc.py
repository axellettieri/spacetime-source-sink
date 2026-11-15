"""Test 1: campionamento completo (emcee) del modello con pozzo proporzionale alla massa accumulata o costante."""
import numpy as np, math, json, sys, emcee
from fit3 import *
name=sys.argv[1] if len(sys.argv)>1 else 'cum_n3'; kind=name.split('_')[0]
mk=MODELS[name][0]
R=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'results3.json'))); pb=np.array(R[name]['p'])
LO=np.array([0.12,0.020,-9.0,-9.0]); HI=np.array([0.17,0.025,-3.0,-0.5])
def lnp(p):
    if np.any(p<LO) or np.any(p>HI): return -np.inf
    c=total(mk,list(p)); return -0.5*c if c<1e8 else -np.inf
nw=24; nst=int(sys.argv[2]) if len(sys.argv)>2 else 2500
rng=np.random.default_rng(1); p0=pb+rng.normal(size=(nw,4))*np.array([5e-4,1e-4,0.03,0.08])
s=emcee.EnsembleSampler(nw,4,lnp); s.run_mcmc(p0,nst,progress=False)
try: tau=s.get_autocorr_time(tol=0); print('autocorr',tau)
except Exception as e: tau=None; print(e)
ch=s.get_chain(discard=nst//4,flat=True); lp=s.get_log_prob(discard=nst//4,flat=True)
np.save(f'chain_{name}.npy',np.c_[ch,lp])
sel=ch[rng.choice(len(ch),600,replace=False)]; der=[]
for p in sel:
    cb,csn,ccm,hh,o=total(mk,list(p),True); f=future(list(p),3,kind,o,hh,tmax=400)
    der.append([hh*100,o[-1,0]/hh**2,f['t_inversione_Gyr'] or np.nan,f['t_crunch_Gyr'] or np.nan,f.get('a_max',np.nan)])
der=np.array(der); q=lambda x: [round(float(v),2) for v in np.nanpercentile(x,[2.5,16,50,84,97.5])]
out=dict(n=len(ch),accept=float(np.mean(s.acceptance_fraction)),best_chi2=float(-2*lp.max()),
  ln_gamma=q(ch[:,3]),H0=q(der[:,0]),Om0=q(der[:,1]),t_inversione=q(der[:,2]),t_crunch=q(der[:,3]),a_max=q(der[:,4]),
  frazione_con_crunch=float(np.mean(np.isfinite(der[:,3]))))
print(json.dumps(out)); json.dump(out,open(f'mcmc_{name}.json','w'),indent=1)
