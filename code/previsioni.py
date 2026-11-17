"""Predictions of the source-sink model (emergent clock) for future DESI and Euclid measurements.
For a random subset of the uncertainty chain it computes, at the redshifts of the DESI tracers:
  D_M/r_d, D_H/r_d, D_V/r_d (baryon acoustic oscillations),
  f*sigma8 (growth rate, redshift-space distortions), plus sigma8, S8 and H0.
The same quantities are computed for the best-fit LCDM model (full Planck + DESI DR2 + Pantheon+), to show
where the two models differ.
Usage: python previsioni.py <chain .npy> <clock .json> <threshold> [samples] [processes]
Writes previsioni_<clock>_<threshold>.json"""
import sys, os, json, math, time
import numpy as np
from multiprocessing import Pool
_arg=sys.argv[:]; sys.argv=sys.argv[:1]
import cmb_orologio as co
sys.argv=_arg
import camb
from scipy.interpolate import CubicSpline, interp1d

ZB=[(0.295,'BGS'),(0.510,'LRG1'),(0.706,'LRG2'),(0.934,'LRG3+ELG1'),(1.321,'ELG2'),(1.484,'QSO'),(2.330,'Lya')]
ZF=[0.295,0.51,0.706,0.934,1.0,1.2,1.321,1.4,1.484,1.65]     # f*sigma8: DESI tracers and Euclid spectroscopic bins
CK=co.CK

def distances(e_u,rd):
    U=co.U; eu=interp1d(U[::-1],e_u[::-1],kind='cubic')
    zg=np.concatenate([[0],np.expm1(np.linspace(math.log(1.0005),math.log(4.0),900))])
    ie=1/eu(np.log1p(zg)); I=np.concatenate([[0],np.cumsum(np.diff(zg)*(ie[1:]+ie[:-1])/2)]); DM=interp1d(zg,I,kind='cubic')
    out={}
    for z,lab in ZB:
        dm=CK*float(DM(z))/rd; dh=CK/float(eu(math.log1p(z)))/rd; dv=(z*dm*dm*dh)**(1/3)
        out[lab]=dict(z=z,DM=dm,DH=dh,DV=dv)
    return out

def sigma8_camb(e_u,wm,wb,logA,ns,tau,zs):
    """sigma8(z) of the clustering matter in CAMB, with the effective fluid of the model."""
    hh=float(e_u[-1]); p=camb.CAMBparams(); a,w=co.w_table(e_u,wm); de=co.DarkEnergyPPF(); de.set_w_a_table(a,w); p.DarkEnergy=de
    p.set_cosmology(H0=100*hh,ombh2=wb,omch2=wm-wb-co.WNU,mnu=0.06,omk=0,tau=tau)
    p.InitPower.set_params(As=math.exp(logA)*1e-10,ns=ns); p.set_matter_power(redshifts=sorted(set([0.0]+list(zs)),reverse=True),kmax=2.0)
    p.WantCls=False; r=camb.get_results(p)
    zz=np.array(sorted(set([0.0]+list(zs)),reverse=True)); s8=r.get_sigma8()
    return dict(zip([float(z) for z in zz],[float(v) for v in s8])), float(r.get_derived_params()['rdrag'])

def growth_total(e_u,m,g,wm):
    """ratio (total-matter contrast in the model)/(contrast of the clustering matter in CAMB) and growth rate, vs ln a"""
    N=co.NLN; dM=co._grow(e_u,m,g); dE=co._grow(e_u,wm*(1+co.ZU)**3,0*g)
    sm=CubicSpline(N,m); se=CubicSpline(N,e_u); sg=CubicSpline(N,g)
    return dM,dE,sm,se,sg

def one(p,Fs,lcdm=False):
    wm,wb=p[0],p[1]; logA,ns,tau=p[-3],p[-2],p[-1]
    if lcdm:
        e=co.e_fluid(wm,p[2]); m=wm*(1+co.ZU)**3; g=0*m
    else:
        e,o,Q=co.passato(wm,math.exp(p[2]),math.exp(p[3]),Fs); m=o[:,0]; g=-Q/m
    hh=float(e[-1]); s8c,rd=sigma8_camb(e,wm,wb,logA,ns,tau,ZF)
    dM,dE,sm,se,sg=growth_total(e,m,g,wm)
    def s8m(z):   # sigma8 of the total matter in the model
        n=-math.log1p(z)
        return s8c[z if z in s8c else min(s8c,key=lambda k:abs(k-z))]*dM(n)[0]/dE(n)[0] if not lcdm else s8c[z if z in s8c else min(s8c,key=lambda k:abs(k-z))]
    fs8={}
    for z in ZF:
        n=-math.log1p(z); d,dp=dM(n); dd,ddp=dE(n)
        # f*sigma8 = d sigma8 / d ln a of the total matter, with the velocity term of the exchange (as in growth.py)
        f=(dp+float(sg(n))/float(se(n))*d)/d if not lcdm else ddp/dd
        fs8[str(z)]=f*s8m(z)
    Om0=float(m[-1])/hh**2; s80=s8m(0.0)
    return dict(H0=100*hh,Om0=Om0,rd=rd,sigma8=s80,S8=s80*math.sqrt(Om0/0.3),s8Om=s80*Om0,bao=distances(e,rd),fs8=fs8)

def worker(args):
    p,f,s,lcdm=args
    try: return one(p,co.carica_orologio(f,s),lcdm)
    except Exception as ex: return None

if __name__=='__main__':
    fc,f,s=sys.argv[1],sys.argv[2],float(sys.argv[3]); ns_=int(sys.argv[4]) if len(sys.argv)>4 else 300; npr=int(sys.argv[5]) if len(sys.argv)>5 else 2
    C=np.load(fc); burn=C.shape[0]//3; X=C[burn:,:,:7].reshape(-1,7); rng=np.random.default_rng(7); X=X[rng.choice(len(X),ns_,replace=False)]
    L=json.load(open('cmb2_LCDM_U.json'))['p']
    t0=time.time()
    ref=one(L,None,lcdm=True); print('LCDM reference done in %.0f s'%(time.time()-t0),flush=True)
    with Pool(npr) as pool: R=[r for r in pool.map(worker,[(x,f,s,False) for x in X]) if r is not None]
    print('samples ok: %d of %d in %.0f s'%(len(R),ns_,time.time()-t0),flush=True)
    def summ(vals):
        v=np.array(vals); q=np.percentile(v,[2.5,16,50,84,97.5]); return [float(x) for x in q]
    out=dict(n=len(R),lcdm=ref,model={})
    for k in ('H0','Om0','rd','sigma8','S8','s8Om'): out['model'][k]=summ([r[k] for r in R])
    out['model']['bao']={lab:{q:summ([r['bao'][lab][q] for r in R]) for q in ('DM','DH','DV')} for _,lab in ZB}
    out['model']['fs8']={str(z):summ([r['fs8'][str(z)] for r in R]) for z in ZF}
    nome=os.path.basename(f).replace('.json','')
    json.dump(out,open('previsioni_%s_%g.json'%(nome,s),'w'),indent=1)
    D={}
    Cv=np.loadtxt('data/desi_cov.txt'); i=0
    for l in open('data/desi_mean.txt'):
        if l.startswith('#'): continue
        z,v,q=l.split(); D[(round(float(z),3),q[:2])]=(float(v),math.sqrt(Cv[i,i])); i+=1
    out['desi_dr2']={'%s_%s'%k:v for k,v in D.items()}
    json.dump(out,open('previsioni_%s_%g.json'%(nome,s),'w'),indent=1)
    print('quantity | model: median [68%] | LCDM best fit | model-LCDM % | DESI DR2 measured +- error')
    for _,lab in ZB:
        for q in ('DM','DH','DV'):
            a=out['model']['bao'][lab][q]; b=ref['bao'][lab][q]
            m=D.get((round(out['model']['bao'][lab]['DM'][2] and [zz for zz,l2 in ZB if l2==lab][0],3),q),None)
            print('%-10s %s/rd  %.3f [%.3f-%.3f] | %.3f | %+.2f%% | %s'%(lab,q,a[2],a[1],a[3],b,100*(a[2]/b-1),'%.3f +- %.3f'%m if m else '-'))
    for z in ZF:
        a=out['model']['fs8'][str(z)]; b=ref['fs8'][str(z)]
        print('fs8 z=%.3f  %.4f [%.4f-%.4f] | %.4f | %+.2f%%'%(z,a[2],a[1],a[3],b,100*(a[2]/b-1)))
    for k in ('H0','Om0','sigma8','S8','s8Om'):
        a=out['model'][k]; print('%-7s %.4f [%.4f-%.4f] | %.4f'%(k,a[2],a[1],a[3],ref[k]))
