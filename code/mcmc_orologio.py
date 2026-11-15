"""Catena completa delle incertezze (emcee) per il modello con OROLOGIO EMERGENTE:
spettro completo di Planck + DESI + supernovae, come cmb_orologio.py.
Uso:  python mcmc_orologio.py <orologio.json> <soglia> <passi> <processi>
      es.  python mcmc_orologio.py fine_L250_s1.json 10 3000 14
Parametri: omega_m, omega_b, ln(alfa), ln(gamma), ln(1e10 A_s), n_s, tau.
Salva ogni 20 passi in catena_orologio_<orologio>_<soglia>[_<supernovae>].npy; rilanciato con lo stesso comando riprende da li'.
Campione di supernovae: variabile SNSET (pantheon predefinito, union3, dovekie)."""
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
    if os.path.exists(f): old=np.load(f); p0=old[-1,:,:nd]; print('riprendo da',f,old.shape,flush=True)
    else: old=None; p0=x0+rng.normal(size=(nw,nd))*sg
    print('partenza:',np.round(x0,5).tolist(),'| camminatori',nw,'| processi',npr,flush=True)
    with Pool(npr) as pool:
        s=emcee.EnsembleSampler(nw,nd,lnp,pool=pool); t0=time.time()
        for i,_ in enumerate(s.sample(p0,iterations=nst)):
            if (i+1)%20==0 or i+1==nst:
                ch=np.concatenate([s.get_chain(),s.get_log_prob()[:,:,None]],axis=2)
                np.save(f+'.tmp.npy',ch if old is None else np.concatenate([old,ch])); os.replace(f+'.tmp.npy',f)
                print('passo %d/%d  %.0f s (%.1f s a passo)  chi2 minimo %.2f'%(i+1,nst,time.time()-t0,(time.time()-t0)/(i+1),-2*s.get_log_prob().max()),flush=True)
