"""Experiment 5: like experiment 4, plus DECAY: a composite object (mass >= 2) that becomes hot enough again
(energy per unit mass >= TC) is unstable and splits in two, with probability PF at each step.
At the start single particles cannot split; during expansion merged objects cool down and are stable;
during collapse they reabsorb space, heat up and, past the threshold, decay: new objects, new events, new space.
The run does not stop: we look at whether the universe oscillates (cycles) and whether the cycles resemble each other.
--- rules inherited from experiment 4 ---
Space = closed network of triangles (surface of a sphere), starting from the smallest possible graph: 4 vertices ("singularity").
- NP hot particles, all in the 4 initial vertices. Each particle jumps to a neighbouring vertex with probability ~ energy/mass.
- EVENT: two objects in the same vertex. They lose a fraction (1-C) of their thermal energy; that energy becomes SPACE:
  one new vertex per EPS of energy. With probability PM the two merge into a single mass.
- ABSORPTION: each object removes neighbouring vertices at rate G per unit mass and takes back their energy (EPS per vertex).
  The matter that was in the removed vertex falls onto the object.
- RELAXATION: edge flips around the touched vertices, towards the regular network.
Nothing is created and nothing is destroyed: total mass and total energy (thermal + space) are conserved exactly.
Usage: python esperimento5.py NP EPS C PM G steps TC PF [seed]     [seed MCRIT TCRIT]   (TC=0: decay only at the singularity; TC=99: never;
     MCRIT: fraction of the total mass above which an object decays completely; TCRIT: minimum energy per unit mass)"""
import random, math, sys, time
from esperimento3 import Rete
from esperimento4 import Sfera, dimensione
def esegui(NP=300,EPS=0.02,C=0.9,PM=0.05,G=0.005,passi=400000,TC=0.9,PF=0.01,seme=1,MCRIT=0.0,TCRIT=0.0,RF=20,TG=0.3,verbose=True,osserva=None):
    random.seed(seme); R=Sfera()
    pos=[random.randrange(4) for _ in range(NP)]; m=[1]*NP; e=[1.0]*NP; vivo=[True]*NP
    occ={v:[] for v in range(4)}
    for i,v in enumerate(pos): occ[v].append(i)
    liberi=[]; nfr=0; critici=[]; vmin_dopo=None
    resto=0.0; Etot0=sum(e); storia=[]; ev=0; t0=time.time(); vmax=0; dim_max=float('nan')
    def rilassa(vv):
        for v in vv:
            if v not in R.nbr: continue
            loc=[v]+list(R.nbr[v])
            for _ in range(RF):
                x=random.choice(loc); y=random.choice(tuple(R.nbr[x])); R.flip(x,y,TG)
    for t in range(passi+1):
        V=len(R.nbr); ogg=[i for i in range(len(vivo)) if vivo[i]]
        if t%max(1,passi//2000)==0:
            espazio=(V-4)*EPS+resto; storia.append((t,V,len(ogg),max(m[i] for i in ogg),sum(e[i] for i in ogg),espazio,ev))
            if osserva: osserva(t,R,[(pos[i],m[i],e[i]) for i in ogg])   # optional measurements (esperimento9.py); does not change the dynamics
        if V>vmax: vmax=V
        if verbose and t%max(1,passi//20)==0:
            print('step %6d: vertices %6d, objects %4d, maximum mass %4d, thermal energy %.1f + space %.1f = %.1f, events %d, dimension %.2f  [%.0f s]'%(t,V,len(ogg),max(m[i] for i in ogg),sum(e[i] for i in ogg),(V-4)*EPS+resto,sum(e[i] for i in ogg)+(V-4)*EPS+resto,ev,dimensione(R),time.time()-t0),flush=True)
        # 1) motion
        for i in ogg:
            if random.random()<min(0.5,0.5*e[i]/m[i]):
                v=pos[i]; w=random.choice(tuple(R.nbr[v])); occ[v].remove(i); occ.setdefault(w,[]).append(i); pos[i]=w
        # 2) events
        toccati=[]
        for v in [v for v,l in occ.items() if len(l)>=2]:
            l=occ[v]; i,j=l[0],l[1]; ev+=1
            dE=(1-C)*(e[i]+e[j]); e[i]*=C; e[j]*=C; resto+=dE
            n=int(resto/EPS); resto-=n*EPS; nuovi=[]
            for _ in range(n):
                nuovi.append(R.split(v,random.choice(tuple(R.nbr[v]))))
            for w in nuovi: occ[w]=[]
            if random.random()<PM:
                m[i]+=m[j]; e[i]+=e[j]; vivo[j]=False; l.remove(j); liberi.append(j)
            elif nuovi:
                w=nuovi[-1]; l.remove(j); occ[w].append(j); pos[j]=w        # created space separates the two
            toccati+=nuovi+[v]
        # 3) absorption
        for i in ogg:
            if not vivo[i]: continue
            k=G*m[i]; n=int(k)+(random.random()<k-int(k))
            for _ in range(n):
                v=pos[i]
                if len(R.nbr)<=4: break
                u=random.choice(tuple(R.nbr[v])); vicini=list(R.nbr[u])
                if R.togli(u):
                    for j in occ.pop(u,[]): occ[v].append(j); pos[j]=v
                    e[i]+=EPS; toccati+=vicini
        # 4) decay. TC<=0: singularity rule -> only when space has returned to the minimal graph, and then every
        #    composite object decays completely into single particles (same energy per unit mass)
        if TC<=0 and len(R.nbr)<=4:
            for i in ogg:
                if vivo[i] and m[i]>=2:
                    eu=e[i]/m[i]; n=m[i]-1; m[i]=1; e[i]=eu; nfr+=1
                    for _ in range(n):
                        w=random.choice(tuple(R.nbr))
                        if liberi: j=liberi.pop(); m[j]=1; e[j]=eu; vivo[j]=True; pos[j]=w
                        else: j=len(vivo); m.append(1); e.append(eu); vivo.append(True); pos.append(w)
                        occ.setdefault(w,[]).append(j)
        # CRITICAL MASS: an object exceeding the fraction MCRIT of the total mass (and at least TCRIT hot) decays completely
        if MCRIT>0:
            for i in ogg:
                if vivo[i] and m[i]>=MCRIT*NP and e[i]>=TCRIT*m[i]:
                    eu=e[i]/m[i]; n=m[i]-1; m[i]=1; e[i]=eu; nfr+=1; v=pos[i]; qui=[v]+list(R.nbr[v]); critici.append((t,len(R.nbr),n+1,eu))
                    for _ in range(n):
                        w=random.choice(qui)
                        if liberi: j=liberi.pop(); m[j]=1; e[j]=eu; vivo[j]=True; pos[j]=w
                        else: j=len(vivo); m.append(1); e.append(eu); vivo.append(True); pos.append(w)
                        occ.setdefault(w,[]).append(j)
        # otherwise: decay of composite objects that became hot again
        for i in ogg:
            if TC>0 and vivo[i] and m[i]>=2 and e[i]>=TC*m[i] and random.random()<PF:
                m2=m[i]//2; e2=e[i]*m2/m[i]; m[i]-=m2; e[i]-=e2; v=pos[i]; w=random.choice(tuple(R.nbr[v])); nfr+=1
                if liberi: j=liberi.pop(); m[j]=m2; e[j]=e2; vivo[j]=True; pos[j]=w
                else: j=len(vivo); m.append(m2); e.append(e2); vivo.append(True); pos.append(w)
                occ.setdefault(w,[]).append(j)
        rilassa(toccati)
    ogg=[i for i in range(len(vivo)) if vivo[i]]; V=len(R.nbr)
    res=dict(decadimenti=nfr,storia=storia,vmax=vmax,V=V,oggetti=len(ogg),massa=sum(m[i] for i in ogg),energia=sum(e[i] for i in ogg)+(V-4)*EPS+resto,E0=Etot0,passi=t,eventi=ev)
    # cycles: maxima and minima of space (smoothed series)
    Vs=[q[1] for q in storia]; w=5; sm=[sum(Vs[max(0,a-w):a+w+1])/len(Vs[max(0,a-w):a+w+1]) for a in range(len(Vs))]
    amp=max(sm)-min(sm); estremi=[]; verso=1; ref=sm[0]; iref=0
    for a,x in enumerate(sm):
        if verso==1:
            if x>ref: ref=x; iref=a
            elif x<ref-0.25*amp: estremi.append(('max',storia[iref][0],int(ref),storia[iref][2])); verso=-1; ref=x; iref=a
        else:
            if x<ref: ref=x; iref=a
            elif x>ref+0.25*amp: estremi.append(('min',storia[iref][0],int(ref),storia[iref][2])); verso=1; ref=x; iref=a
    res['estremi']=estremi
    if verbose:
        print('final mass %d, final energy %.1f (initial %.1f), decays %d, events %d'%(res['massa'],res['energia'],Etot0,nfr,ev))
        if critici:
            print('DECAYS BY CRITICAL MASS (step, vertices at that moment, mass, energy per unit mass):')
            for c in critici: print('   step %7d: %6d vertices, mass %3d, energy/mass %.2f'%c)
            print('   minimum space after the first maximum: %d vertices (4 = singularity)'%min(q[1] for q in storia[5:]))
        print('CYCLES (smoothed space): '+('no clear oscillation' if len(estremi)<3 else ''))
        for tipo,tt,vv,oo in estremi: print('   %s at step %7d: %6d vertices, %4d objects'%('MAXIMUM' if tipo=='max' else 'minimum ',tt,vv,oo))
    return res
if __name__=='__main__':
    a=sys.argv[1:]; f=lambda k,d,ty: ty(a[k]) if len(a)>k else d
    esegui(f(0,300,int),f(1,0.02,float),f(2,0.9,float),f(3,0.05,float),f(4,0.005,float),f(5,400000,int),f(6,0.9,float),f(7,0.01,float),f(8,1,int),f(9,0.0,float),f(10,0.0,float))
