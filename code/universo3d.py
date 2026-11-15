"""Universo in tre dimensioni con gravita': lo spazio nasce dagli eventi e viene assorbito dalle masse.
Simulazione a N corpi (particelle su griglia, metodo PM) su GPU con CuPy (senza GPU usa NumPy: solo per prove piccole).
- La materia si addensa per gravita'. EVENTO = una particella entra per la prima volta in una regione collassata
  (densita' > DC volte la media): e' l'orologio degli eventi, non piu' imposto ma prodotto dalla simulazione.
- Ogni evento crea energia del vuoto (spazio nuovo):  alpha * massa collassata nuova.
- Le masse collassate la riassorbono:                 gamma * massa collassata totale, per unita' di tempo.
- LAM = raggio (Mpc) su cui lo spazio creato si sparge attorno all'evento. 'inf' = sparso ovunque (versione diffusa).
  Con LAM finito lo spazio in piu' resta vicino agli eventi e spinge via la materia: espansione locale diversa.
- L'espansione globale NON e' imposta: viene dall'equazione di Friedmann con la materia e lo spazio creato medi.
Unita': distanze in Mpc, H in 100 km/s/Mpc, tempo in 9.78 miliardi di anni, densita' in omega = Omega*h^2.
Uscite: <nome>.json (storia dell'espansione, orologio degli eventi, espansione locale attorno a osservatori)."""
import numpy as np, math, json, sys, os, time, argparse
try:
    if os.environ.get('CPU'): raise ImportError
    import cupy as xp; GPU=True
except ImportError:
    xp=np; GPU=False
def tonp(a): return xp.asnumpy(a) if GPU else a
OR=4.184e-5
trapz=getattr(np,"trapezoid",None) or np.trapz
ap=argparse.ArgumentParser()
ap.add_argument('nome'); ap.add_argument('--np',type=int,default=512); ap.add_argument('--ng',type=int,default=512)
ap.add_argument('--L',type=float,default=1000.); ap.add_argument('--passi',type=int,default=600)
ap.add_argument('--alpha',type=float,default=0.); ap.add_argument('--gamma',type=float,default=0.)
ap.add_argument('--lam',default='inf'); ap.add_argument('--xinit',type=float,default=0.)
ap.add_argument('--wm',type=float,default=0.1430); ap.add_argument('--wb',type=float,default=0.02236)
ap.add_argument('--dc',type=float,default=20.); ap.add_argument('--ainiz',type=float,default=0.02); ap.add_argument('--afin',type=float,default=1.0)
ap.add_argument('--soglie',default='',help='altre soglie di collasso da registrare come orologi, es. 2,3,5,10,50')
ap.add_argument('--pozzo',default='materia',help="materia: assorbimento proporzionale alla materia collassata attuale (si autoalimenta); originale: alla materia collassata iniziale, come nel modello sorgente-pozzo adattato ai dati")
ap.add_argument('--collasso',action='store_true',help="prosegue oltre l'inversione, fino al collasso (a < acoll)")
ap.add_argument('--acoll',type=float,default=0.1); ap.add_argument('--dtmax',type=float,default=0.004)
ap.add_argument('--seme',type=int,default=1); ap.add_argument('--par',default=None); ap.add_argument('--salva',type=int,default=100)
A=ap.parse_args()
if A.par:
    P=json.load(open(A.par)); A.alpha,A.gamma,A.wm,A.wb,A.dc=P['alpha'],P['gamma'],P['wm'],P['wb'],P.get('dc',A.dc); A.pozzo=P.get('pozzo',A.pozzo)
NP,NG,L=A.np,A.ng,A.L; N=NP**3; M=NG**3; LAM=None if A.lam=='inf' else float(A.lam); CH=1<<24
f32=xp.float32
# ---------- spettro iniziale (Eisenstein-Hu senza oscillazioni), normalizzato a sigma8=0.811 oggi in LCDM ----------
H_FID=0.674; S8=0.811; NS=0.965
def transfer(k,wm,wb,h=H_FID):
    fb=wb/wm; s=44.5*math.log(9.83/wm)/math.sqrt(1+10*wb**0.75); ag=1-0.328*math.log(431*wm)*fb+0.38*math.log(22.3*wm)*fb**2
    G=wm/h*(ag+(1-ag)/(1+(0.43*k*s)**4)); q=k/h*(2.7255/2.7)**2/G
    L0=np.log(2*math.e+1.8*q); C0=14.2+731/(1+62.5*q); return L0/(L0+C0*q*q)
def growth_lcdm(a,wm,h=H_FID):
    om=wm/h**2; aa=np.linspace(1e-4,a,4000); E=np.sqrt(om/aa**3+1-om); return E[-1]*trapz(1/(aa*E)**3,aa)
def pk_norm(wm,wb):
    k=np.logspace(-4,2,4000); R=8/H_FID; x=k*R; W=3*(np.sin(x)-x*np.cos(x))/x**3
    return S8**2/trapz(k**NS*transfer(k,wm,wb)**2*W**2*k**3/(2*math.pi**2),np.log(k))
def kgrid(n):
    k1=(2*math.pi*xp.fft.fftfreq(n,d=L/n)).astype(f32); kz=(2*math.pi*xp.fft.rfftfreq(n,d=L/n)).astype(f32)
    return k1[:,None,None],k1[None,:,None],kz[None,None,:]
# ---------- particelle <-> griglia (CIC, a blocchi per non riempire la memoria) ----------
def corners(p):
    g=p*(NG/L); i0=xp.floor(g).astype(xp.int32); d=(g-i0).astype(f32); i0%=NG; i1=(i0+1)%NG
    for cx in (0,1):
        ix=(i1 if cx else i0)[:,0].astype(xp.int64); wx=d[:,0] if cx else 1-d[:,0]
        for cy in (0,1):
            iy=(i1 if cy else i0)[:,1]; wy=d[:,1] if cy else 1-d[:,1]
            for cz in (0,1):
                iz=(i1 if cz else i0)[:,2]; wz=d[:,2] if cz else 1-d[:,2]
                yield (ix*NG+iy)*NG+iz, wx*wy*wz
def deposit(pos,sel=None):
    rho=xp.zeros(M,xp.float64)
    for s in range(0,N,CH):
        p=pos[s:s+CH]
        if sel is not None:
            p=p[sel[s:s+CH]]
            if len(p)==0: continue
        for idx,w in corners(p): rho+=xp.bincount(idx,weights=w,minlength=M)
    return (rho*(M/N)).astype(f32)                     # media 1 se si depositano tutte le particelle
def gather(field,p):
    out=xp.zeros(len(p),f32)
    for idx,w in corners(p): out+=field[idx]*w
    return out
def ngp(p):
    i=xp.floor(p*(NG/L)).astype(xp.int64)%NG; return (i[:,0]*NG+i[:,1])*NG+i[:,2]
# ---------- espansione locale attorno a osservatori (il test della tensione di Hubble) ----------
def hubble_locale(pos,mom,flag,a,E,raggi,nobs=300,ntr=300000,seme=7):
    rng=np.random.default_rng(seme); ii=tonp(xp.nonzero(flag)[0])
    if len(ii)<1000: return None
    ii=np.sort(rng.choice(ii,min(ntr,len(ii)),replace=False)); gi=xp.asarray(ii)
    X=tonp(pos[gi]).astype(np.float64); V=tonp(mom[gi]).astype(np.float64)/a      # velocita' peculiari (x100 km/s)
    out={}; res={R:[] for R in raggi}
    for o in rng.choice(len(X),min(nobs,len(X)),replace=False):
        d=X-X[o]; d-=L*np.round(d/L); r2=(d*d).sum(1); vr=((V-V[o])*d).sum(1)
        for R in raggi:
            m=(r2>20.**2)&(r2<R*R)
            if m.sum()>30: res[R].append(vr[m].sum()/r2[m].sum()/E)                 # dH/H
    for R in raggi:
        v=np.array(res[R]); out[str(R)]=dict(media_pct=float(100*v.mean()),scarto_pct=float(100*v.std()),n=len(v))
    return out
# ---------- condizioni iniziali (Zel'dovich) oppure ripresa ----------
ck=A.nome+'_stato.npz'; hist=[]; hub={}
dl=math.log(A.afin/A.ainiz)/A.passi
if os.path.exists(ck):
    Z=np.load(ck); pos=xp.asarray(Z['pos']); mom=xp.asarray(Z['mom']); flag=xp.asarray(Z['flag']); dxf=xp.asarray(Z['dxf']) if LAM else None
    n0,xbar,Mc,F=int(Z['n']),float(Z['xbar']),float(Z['Mc']),float(Z['F']); hist=Z['hist'].tolist(); hub=json.loads(str(Z['hub']))
    RIP=dict(modo=str(Z['modo']),adot=float(Z['adot']),aT=float(Z['aT']),tcos=float(Z['tcos']),tempi=Z['tempi'].tolist(),bits=Z['bits'],Fk=Z['Fk'].tolist()) if 'modo' in Z else None
    print('ripresa dal passo',n0,flush=True)
else:
    t0=time.time(); a=A.ainiz; amp=pk_norm(A.wm,A.wb)*(growth_lcdm(a,A.wm)/growth_lcdm(1.0,A.wm))**2
    kx,ky,kz=kgrid(NP); k2=kx*kx+ky*ky+kz*kz; k2[0,0,0]=1
    kk=tonp(xp.sqrt(k2)).astype(np.float64); Pk=amp*kk**NS*transfer(kk,A.wm,A.wb)**2; del kk
    w=xp.asarray(np.random.default_rng(A.seme).standard_normal((NP,NP,NP),dtype=np.float32))
    dk=(xp.fft.rfftn(w)*xp.asarray(np.sqrt(Pk/(L/NP)**3).astype(np.float32))).astype(xp.complex64); dk[0,0,0]=0; del w,Pk
    E=math.sqrt(A.wm/a**3+OR/a**4+A.xinit)
    pos=xp.empty((N,3),f32); mom=xp.empty((N,3),f32); q=((xp.arange(NP,dtype=f32)+0.5)*(L/NP))
    for ax,kax in enumerate((kx,ky,kz)):
        psi=xp.fft.irfftn(1j*kax/k2*dk,s=(NP,)*3).astype(f32).ravel()
        sh=[1,1,1]; sh[ax]=NP; pos[:,ax]=xp.broadcast_to(q.reshape(sh),(NP,)*3).ravel()+psi; mom[:,ax]=(a*a*E)*psi
    print('condizioni iniziali: spostamento rms %.3f Mpc (cella %.2f Mpc) in %.0f s'%(float(xp.sqrt((mom**2).mean()*3))/(a*a*E),L/NG,time.time()-t0),flush=True)
    del dk,psi,k2; pos%=L; flag=xp.zeros(N,bool); dxf=xp.zeros(M,f32) if LAM else None
    n0=0; xbar=A.xinit; Mc=A.wm; F=0.0; RIP=None
kx,ky,kz=kgrid(NG); k2=kx*kx+ky*ky+kz*kz; k2[0,0,0]=1; ik2=(1/k2).astype(f32); ik2[0,0,0]=0
smooth=xp.exp(-k2*(LAM**2/2)).astype(f32) if LAM else None; del k2
TH=[float(v) for v in A.soglie.split(',') if v]; bits=xp.zeros(N,xp.uint8) if TH else None; Fk=[0.0]*len(TH)
raggi=[R for R in (75,150,300,450) if R<L/2]
print('GPU' if GPU else 'CPU (NumPy)','| particelle %d^3, griglia %d^3, scatola %.0f Mpc, alpha %.4g gamma %.4g lam %s'%(NP,NG,L,A.alpha,A.gamma,A.lam),flush=True)
# ---------- ciclo principale ----------
# Fase F: il passo e' fissato in ln(a) e H viene dall'equazione di Friedmann (come nelle corse precedenti).
# Fase T (solo con --collasso): vicino all'inversione e durante la contrazione si integra nel tempo
#   l'equazione dell'accelerazione  a'' = -a/2 (materia + 2 radiazione - 2 spazio creato), cosi' H puo' cambiare segno.
t0=time.time(); fine='afin'; modo='F'; adot=0.0; aT=0.0; tcos=0.0; tempi=[]; HSW=0.12
if RIP:
    modo,adot,aT,tcos,tempi,Fk=RIP['modo'],RIP['adot'],RIP['aT'],RIP['tcos'],RIP['tempi'],RIP['Fk']
    if TH: bits=xp.asarray(RIP['bits'])
NMAX=A.passi*10 if A.collasso else A.passi; n=n0; amax=0.0
while True:
    if modo=='F':
        a=A.ainiz*math.exp(n*dl); E2=Mc/a**3+OR/a**4+xbar
        if E2<=0 or (A.collasso and E2<HSW**2):
            if not A.collasso: fine='espansione ferma (inversione)'; break
            modo='T'; aT=a; adot=a*math.sqrt(max(E2,1e-12)); print('  -> fase T (integrazione nel tempo) da a=%.3f'%a,flush=True)
    if modo=='F': H=math.sqrt(E2); da=a*(math.exp(dl)-1); dt=da/(a*H)
    else: a=aT; H=adot/a; dt=min(dl/max(abs(H),1e-9),A.dtmax)
    amax=max(amax,a)
    rho=deposit(pos); rms=float(rho.std())
    if os.environ.get('DEBUG'): rms=float(xp.sqrt((xp.abs(xp.fft.rfftn(rho.reshape(NG,NG,NG))[:4,:4,1:4])**2).sum()))   # solo modi a grande scala (controllo crescita lineare)
    # eventi: particelle che entrano per la prima volta in una regione collassata
    nnew=0; new=xp.zeros(N,bool) if LAM else None
    for s in range(0,N,CH):
        r=rho[ngp(pos[s:s+CH])]; c=(~flag[s:s+CH])&(r>A.dc); nnew+=int(c.sum())
        for j,th in enumerate(TH):
            nb=(r>th)&(((bits[s:s+CH]>>j)&1)==0); Fk[j]+=int(nb.sum())/N; bits[s:s+CH]|=(nb.astype(xp.uint8)<<j)
        del r
        if LAM: new[s:s+CH]=c
        flag[s:s+CH]|=c
    dF=nnew/N; Fold=F; F+=dF
    Ms=Mc if A.pozzo=='materia' else A.wm
    dx=A.alpha*Mc*dF/a**3-A.gamma*Ms*Fold*dt                           # spazio creato dagli eventi - spazio assorbito dalle masse
    if LAM and (A.alpha or A.gamma):
        inc=(A.alpha*Mc/a**3)*deposit(pos,new)-(A.gamma*Ms*dt)*(deposit(pos,flag) if Fold>0 else 0)
        ik=xp.fft.rfftn(inc.reshape(NG,NG,NG))*smooth; ik[0,0,0]=0; dxf+=xp.fft.irfftn(ik,s=(NG,)*3).astype(f32).ravel(); del inc,ik
    del new
    xbar+=dx; Mc-=a**3*dx
    hist.append([a,H,xbar,Mc,F,dF/dt,rms,float(dxf.std()) if LAM else 0.0]+list(Fk)); tempi.append(tcos)
    if modo=='F' and (abs(math.log(a))<dl/2 or n==A.passi):
        h=hubble_locale(pos,mom,flag,a,H,raggi)
        if h: hub['a=%.3f'%a]=h; print('  espansione locale a a=%.3f: '%a+' | '.join('R<%s Mpc: %+.2f%% +- %.2f%%'%(R,v['media_pct'],v['scarto_pct']) for R,v in h.items()),flush=True)
    if modo=='F' and n==A.passi: break
    if modo=='T' and a<=A.acoll: fine='collasso'; break
    if n>=NMAX: fine='passi esauriti'; break
    # gravita': materia + spinta dello spazio in piu' (il vuoto respinge: -2 a^2 dx)
    S=(1.5*Mc/a)*(rho-1)
    if LAM: S-=(3*a*a)*dxf
    phik=-xp.fft.rfftn(S.reshape(NG,NG,NG))*ik2; del S,rho
    for ax,kax in enumerate((kx,ky,kz)):
        g=xp.fft.irfftn(-1j*kax*phik,s=(NG,)*3).astype(f32).ravel()
        for s in range(0,N,CH): mom[s:s+CH,ax]+=gather(g,pos[s:s+CH])*dt
        del g
    del phik
    if modo=='F':
        am=a*math.exp(dl/2); Em2=Mc/am**3+OR/am**4+xbar
        if Em2<=0 and not A.collasso: fine='espansione ferma (inversione)'; break
        drift=da/(am**3*math.sqrt(Em2)) if Em2>0 else dt/am**2
    else:
        adot+=-0.5*a*(Mc/a**3+2*OR/a**4-2*xbar)*dt; am=a+adot*dt/2; drift=dt/am**2; aT=a+adot*dt
    for s in range(0,N,CH): pos[s:s+CH]=(pos[s:s+CH]+mom[s:s+CH]*drift)%L
    tcos+=dt
    if n%max(1,A.passi//40)==0:
        print('passo %5d a=%.4f t=%6.2f Gyr H=%+.4f spazio creato=%+.4f materia=%.4f collassata=%.4f rms densita\'=%.2f  [%.0f s]'%(n,a,tcos*9.7779,H,xbar,Mc,F,rms,time.time()-t0),flush=True)
    n+=1
    if A.salva and n%A.salva==0:
        np.savez(ck+'.tmp.npz',pos=tonp(pos),mom=tonp(mom),flag=tonp(flag),dxf=tonp(dxf) if LAM else 0,n=n,xbar=xbar,Mc=Mc,F=F,hist=np.array(hist),hub=json.dumps(hub),
                 modo=modo,adot=adot,aT=aT,tcos=tcos,tempi=np.array(tempi),bits=tonp(bits) if TH else 0,Fk=np.array(Fk))
        os.replace(ck+'.tmp.npz',ck)
json.dump(dict(par=vars(A),soglie=TH,fine=fine,colonne=['a','H','spazio_creato','materia','frazione_collassata','eventi_per_tempo','rms_densita','rms_spazio_locale'],storia=hist,tempo_Gyr=[t*9.7779 for t in tempi],hubble_locale=hub),open(A.nome+'.json','w'))
if os.path.exists(ck): os.remove(ck)
print('fine:',fine,'| a=%.3f H=%+.4f | dimensione massima %.2f | durata cosmica %.1f Gyr | tempo totale %.0f s'%(hist[-1][0],hist[-1][1],amax,tcos*9.7779,time.time()-t0))
