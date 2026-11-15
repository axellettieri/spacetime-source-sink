"""Experiment 4: matter and events on the graph without coordinates, with exact energy bookkeeping.
Space = closed network of triangles (surface of a sphere), starting from the smallest possible graph: 4 vertices ("singularity").
- NP hot particles, all in the 4 initial vertices. Each particle jumps to a neighbouring vertex with probability ~ energy/mass.
- EVENT: two objects in the same vertex. They lose a fraction (1-C) of their thermal energy; that energy becomes SPACE:
  one new vertex per EPS of energy. With probability PM the two merge into a single mass.
- ABSORPTION: each object removes neighbouring vertices at rate G per unit mass and takes back their energy (EPS per vertex).
  The matter that was in the removed vertex falls onto the object.
- RELAXATION: edge flips around the touched vertices, towards the regular network.
Nothing is created and nothing is destroyed: total mass and total energy (thermal + space) are conserved exactly.
Usage: python esperimento4.py NP EPS C PM G steps [seed]"""
import random, math, sys, time
from esperimento3 import Rete
class Sfera(Rete):
    def __init__(s):
        s.tri={}; s.edge={}; s.nbr={v:set() for v in range(4)}; s.nt=0; s.nv=4
        for t in ((0,1,2),(0,3,1),(0,2,3),(1,3,2)): s.add(*t)
    def togli(s,v):
        """removes vertex v (first brings it to 3 neighbours with edge flips). False if not possible."""
        if len(s.nbr)<=4: return False
        guard=0
        while len(s.nbr[v])>3:
            ok=False
            for u in random.sample(tuple(s.nbr[v]),len(s.nbr[v])):
                if len(s.nbr[u])<=3: continue
                c=s.opp(v,u); d=s.opp(u,v)
                if c==d or d in s.nbr[c]: continue
                s.rem(s.edge[(v,u)]); s.rem(s.edge[(u,v)]); s.nbr[v].discard(u); s.nbr[u].discard(v)
                s.add(v,d,c); s.add(d,u,c); ok=True; break
            guard+=1
            if not ok or guard>200: return False
        p=next(iter(s.nbr[v])); q=s.opp(v,p)            # triangle (v,p,q)
        r=[x for x in s.nbr[v] if x!=p and x!=q][0]
        if len(s.nbr[p])<=3 or len(s.nbr[q])<=3 or len(s.nbr[r])<=3: return False
        s.rem(s.edge[(v,p)]); s.rem(s.edge[(v,q)]); s.rem(s.edge[(v,r)])
        for x in (p,q,r): s.nbr[x].discard(v)
        del s.nbr[v]; s.add(p,q,r); return True
def dimensione(R,nrad=12):
    vs=list(R.nbr); V=len(vs)
    if V<200: return float('nan')
    r2=max(4,int(0.35*math.sqrt(V))); r1=max(2,r2//2); n1=n2=0
    for r0 in random.sample(vs,nrad):
        d=R.bfs(r0,r2); n1+=sum(1 for x in d.values() if x<=r1); n2+=len(d)
    return math.log(n2/n1)/math.log(r2/r1)
def esegui(NP=300,EPS=0.02,C=0.9,PM=0.05,G=0.0005,passi=20000,seme=1,RF=20,TG=0.3,verbose=True):
    random.seed(seme); R=Sfera()
    pos=[random.randrange(4) for _ in range(NP)]; m=[1]*NP; e=[1.0]*NP; vivo=[True]*NP
    occ={v:[] for v in range(4)}
    for i,v in enumerate(pos): occ[v].append(i)
    resto=0.0; Etot0=sum(e); storia=[]; ev=0; t0=time.time(); vmax=0; dim_max=float('nan')
    def rilassa(vv):
        for v in vv:
            if v not in R.nbr: continue
            loc=[v]+list(R.nbr[v])
            for _ in range(RF):
                x=random.choice(loc); y=random.choice(tuple(R.nbr[x])); R.flip(x,y,TG)
    for t in range(passi+1):
        V=len(R.nbr); ogg=[i for i in range(NP) if vivo[i]]
        if t%max(1,passi//200)==0 or (V<=4 and len(ogg)==1):
            espazio=(V-4)*EPS+resto; storia.append((t,V,len(ogg),max(m[i] for i in ogg),sum(e[i] for i in ogg),espazio,ev))
        if V>vmax: vmax=V
        if V<=4 and len(ogg)==1 and t>10: break
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
                m[i]+=m[j]; e[i]+=e[j]; vivo[j]=False; l.remove(j)
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
        rilassa(toccati)
    ogg=[i for i in range(NP) if vivo[i]]; V=len(R.nbr)
    res=dict(storia=storia,vmax=vmax,V=V,oggetti=len(ogg),massa=sum(m[i] for i in ogg),energia=sum(e[i] for i in ogg)+(V-4)*EPS+resto,E0=Etot0,passi=t,eventi=ev)
    if verbose:
        j=max(range(len(storia)),key=lambda q:storia[q][1])
        print('START: 4 vertices, %d objects, mass %d, energy %.1f'%(NP,NP,Etot0))
        print('MAXIMUM: %d vertices at step %d (%d objects)'%(storia[j][1],storia[j][0],storia[j][2]))
        print('END (step %d): %d vertices, %d objects, mass %d, energy %.1f (thermal %.1f, space %.1f)'%(t,V,len(ogg),res['massa'],res['energia'],sum(e[i] for i in ogg),(V-4)*EPS+resto))
    return res
if __name__=='__main__':
    a=sys.argv[1:]; f=lambda k,d,ty: ty(a[k]) if len(a)>k else d
    esegui(f(0,300,int),f(1,0.02,float),f(2,0.9,float),f(3,0.05,float),f(4,0.0005,float),f(5,20000,int),f(6,1,int))
