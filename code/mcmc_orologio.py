"""Full uncertainty chain (emcee) for the model with the EMERGENT CLOCK:
full Planck spectrum + DESI + supernovae, as in cmb_orologio.py.
Usage:  python mcmc_orologio.py <clock.json> <threshold> <steps> <processes>
      e.g.  python mcmc_orologio.py fine_L250_s1.json 10 3000 14
Parameters: omega_m, omega_b, ln(alpha), ln(gamma), ln(1e10 A_s), n_s, tau.
Saves every 20 steps in catena_orologio_<clock>_<threshold>[_<supernovae>].npy; relaunched with the same command it resumes from there.
Supernova sample: SNSET variable (pantheon default, union3, dovekie)."""
import numpy as np, sys, os, json, time
os.environ['OMP_NUM_THREADS']='1'
from multiprocessing import Pool
import emcee
_arg=sys.argv[:]; sys.argv=sys.argv[:1]
import cmb_orologio as co
sys.argv=_arg
f_or=sys.argv[1]; soglia=float(sys.argv[2]); nst=int(sys.argv[3]) if len(sys.argv)>3 else 3000; npr=int(sys.argv[4]) if len(sys.argv)>4 else 8
nw=int(os.environ.get('CAMMINATORI','32'))
FS=co.carica_orologio(f_or,soglia)
def lnp(p):
    c=co.chi2_orol(list(p),FS); return -0.5*c if c<1e8 else -np.inf
if __name__=='__main__':
    suf='' if co.SNSET=='pantheon' else '_'+co.SNSET
    nome=os.path.basename(f_or).replace('.json','')
    start=os.path.join(co.HERE,'cmb_orologio_%s_%g%s.json'%(nome,soglia,suf))
    x0=np.array(json.load(open(start))['p']) if os.path.exists(start) else np.array([0.1424,0.0224,-2.18,0.38,3.0497,0.9667,0.0584])
    sg=np.array([6e-4,1.2e-4,0.04,0.10,0.008,0.004,0.006]); nd=len(x0)
    f=os.path.join(co.HERE,'catena_orologio_%s_%g%s.npy'%(nome,soglia,suf)); rng=np.random.default_rng(3)
    if os.path.exists(f): old=np.load(f); p0=old[-1,:,:nd]; print('resuming from',f,old.shape,flush=True)
    else: old=None; p0=x0+rng.normal(size=(nw,nd))*sg
    print('start:',np.round(x0,5).tolist(),'| walkers',nw,'| processes',npr,flush=True)
    with Pool(npr) as pool:
        s=emcee.EnsembleSampler(nw,nd,lnp,pool=pool); t0=time.time()
        for i,_ in enumerate(s.sample(p0,iterations=nst)):
            if (i+1)%20==0 or i+1==nst:
                ch=np.concatenate([s.get_chain(),s.get_log_prob()[:,:,None]],axis=2)
                np.save(f+'.tmp.npy',ch if old is None else np.concatenate([old,ch])); os.replace(f+'.tmp.npy',f)
                print('step %d/%d  %.0f s (%.1f s per step)  minimum chi2 %.2f'%(i+1,nst,time.time()-t0,(time.time()-t0)/(i+1),-2*s.get_log_prob().max()),flush=True)
