"""Campionamento completo con lo spettro del fondo cosmico (versione 2, da eseguire sul PC).
Uso:  python mcmc_cmb.py cum_n3 3000 8 [U|E]     (modello, passi, processi, trattamento della lente: U = corretto, predefinito)
Modelli: cum_n3 (pozzo proporzionale alla massa accumulata), prop_n3 (pozzo proporzionale al livello), LCDM, w0wa.
Salva ogni 20 passi in catena_v2_<modello>.npy e, rilanciato con lo stesso comando, riprende da li'.
Con 8 processi: circa 4 secondi per valutazione, 32 camminatori -> circa 17 s a passo, 3000 passi in circa 14 ore."""
import numpy as np, sys, os, json, time
os.environ['OMP_NUM_THREADS']='1'
from multiprocessing import Pool
import emcee
from cmb_full import *
name=sys.argv[1] if len(sys.argv)>1 else 'cum_n3'; nst=int(sys.argv[2]) if len(sys.argv)>2 else 3000
npr=int(sys.argv[3]) if len(sys.argv)>3 else 8; pert=sys.argv[4] if len(sys.argv)>4 else 'U'
def lnp(p):
    c=chi2_all(name,list(p),pert=pert); return -0.5*c if c<1e8 else -np.inf
if __name__=='__main__':
    x0=np.array(START[name]); sg=1.5*np.array(SCALE[name]); nd=len(x0); nw=32
    f=os.path.join(HERE,f'catena_v2_{name}.npy'); rng=np.random.default_rng(2)
    if os.path.exists(f): old=np.load(f); p0=old[-1,:,:nd]; print('riprendo da',f,old.shape)
    else: old=None; p0=x0+rng.normal(size=(nw,nd))*sg
    with Pool(npr) as pool:
        s=emcee.EnsembleSampler(nw,nd,lnp,pool=pool); t0=time.time()
        for i,_ in enumerate(s.sample(p0,iterations=nst)):
            if (i+1)%20==0 or i+1==nst:
                ch=np.concatenate([s.get_chain(),s.get_log_prob()[:,:,None]],axis=2)
                np.save(f,ch if old is None else np.concatenate([old,ch]))
                print('passo %d/%d  %.0f s  chi2 minimo %.2f'%(i+1,nst,time.time()-t0,-2*s.get_log_prob().max()),flush=True)
