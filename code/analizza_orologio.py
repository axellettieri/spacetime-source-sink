"""Summary of the mcmc_orologio.py chain.
Usage:  python analizza_orologio.py <clock.json> <threshold> [burn-in in steps]
Prints: central values and 68% and 95% intervals, H0, matter fraction today, equivalent (w0, wa) point,
autocorrelation time (the chain is reliable if the kept steps are at least 30-50 times this time)."""
import numpy as np, sys, os, json, math
_arg=sys.argv[:]; sys.argv=sys.argv[:1]
import cmb_orologio as co
from traccia_w0wa import mappa
sys.argv=_arg
f_or=sys.argv[1]; soglia=float(sys.argv[2]); suf='' if co.SNSET=='pantheon' else '_'+co.SNSET
nome=os.path.basename(f_or).replace('.json','')
C=np.load(os.path.join(co.HERE,'catena_orologio_%s_%g%s.npy'%(nome,soglia,suf))); nst=C.shape[0]
burn=int(sys.argv[3]) if len(sys.argv)>3 else nst//3
try:
    import emcee; tau=emcee.autocorr.integrated_time(C[burn:,:,:7],quiet=True); print('total steps %d, discarded %d; autocorrelation time per parameter: %s'%(nst,burn,np.round(tau,0).tolist()))
    print('kept steps / maximum autocorrelation time = %.0f  (at least 30-50 needed)'%((nst-burn)/tau.max()))
except Exception as e: print('autocorrelation not estimated:',e)
S=C[burn:].reshape(-1,C.shape[2]); P,L=S[:,:7],S[:,7]
LC=json.load(open(os.path.join(co.HERE,'cmb2_LCDM_U%s.json'%suf)))['tot']
print('minimum chi2 in the chain %.2f, difference from LCDM %+.2f'%(-2*L.max(),-2*L.max()-LC))
def riga(nome,v): q=np.percentile(v,[2.5,16,50,84,97.5]); print('%-28s %.5g   68%%: %.5g - %.5g   95%%: %.5g - %.5g'%(nome,q[2],q[1],q[3],q[0],q[4]))
for i,n in enumerate(['omega_m','omega_b','ln alfa','ln gamma','ln(1e10 As)','n_s','tau']): riga(n,P[:,i])
FS=co.carica_orologio(f_or,soglia); rng=np.random.default_rng(1); idx=rng.choice(len(P),min(300,len(P)),replace=False)
H0=[];Om=[];W0=[];WA=[]
for k in idx:
    p=P[k]
    try:
        e,o,Q=co.passato(p[0],math.exp(p[2]),math.exp(p[3]),FS); hh=co.observables(e,p[0],p[1])[3]
        w0,wa,_=mappa(e,p[0]); H0.append(100*hh); Om.append(o[-1,0]/hh**2); W0.append(w0); WA.append(wa)
    except Exception: pass
riga('H0 (km/s/Mpc)',np.array(H0)); riga('matter today (fraction)',np.array(Om)); riga('equivalent w0',np.array(W0)); riga('equivalent wa',np.array(WA))
json.dump(dict(H0=H0,Om=Om,w0=W0,wa=WA),open(os.path.join(co.HERE,'derivati_orologio_%s_%g%s.json'%(nome,soglia,suf)),'w'))
