import numpy as np, math, json
from fit3 import *
R=json.load(open('results3.json')); out={}
for name in ('cum_n3','const_n3'):
    kind=name.split('_')[0]; mk=MODELS[name][0]; pb=R[name]['p']; rows=[]
    for dl in np.arange(-1.8,1.21,0.15):
        lg=pb[3]+dl; b=fitmodel(mk,[pb[:3]+[lg]],fixed=lg); p=list(b.x)+[lg]
        cb,csn,ccm,hh,o=total(mk,p,True); f=future(p,3,kind,o,hh,tmax=200)
        rows.append(dict(lg=lg,chi2=float(b.fun),dchi2=float(b.fun-R[name]['chi2']),h=hh,Om0=float(o[-1,0]/hh**2),t_inv=f['t_inversione_Gyr'],t_crunch=f['t_crunch_Gyr'],a_max=f.get('a_max')))
        print(name,'dlg=%.2f dchi2=%.2f h=%.4f Om0=%.3f inv=%s crunch=%s amax=%s'%(dl,rows[-1]['dchi2'],hh,rows[-1]['Om0'],f['t_inversione_Gyr'],f['t_crunch_Gyr'],f.get('a_max')),flush=True)
    out[name]=rows
json.dump(out,open('profile.json','w'),indent=1)
