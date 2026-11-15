"""Experiment 3: space as a graph without coordinates (two space dimensions).
Space is a closed network of triangles: each vertex is an atom of space, distances are counted in steps on the graph.
- Local CREATION: in an active region (where events happen) vertices are added, one at a time, splitting a triangle into three.
- RELAXATION: edge flips accepted if they bring the network closer to the regular one (6 neighbours per vertex = flat space),
  with some thermal agitation T. R = flip attempts per created vertex.
Questions: does space created at a point stay there or spread? Do distances really grow? Does space stay two-dimensional?
Usage: python esperimento3.py n growth R T [seed] [ovunque]     (n = side of the initial network, growth = fraction of vertices added)"""
import random, math, sys, time
from collections import deque
class Rete:
    def __init__(s,n):
        s.tri={}; s.edge={}; s.nbr={}; s.nt=0; s.n=n
        vid=lambda i,j:(i%n)*n+(j%n)
        for v in range(n*n): s.nbr[v]=set()
        for i in range(n):
            for j in range(n):
                a,b,c,d=vid(i,j),vid(i+1,j),vid(i,j+1),vid(i+1,j+1)
                s.add(a,b,c); s.add(b,d,c)
        s.nv=n*n
    def add(s,a,b,c):
        t=s.nt; s.nt+=1; s.tri[t]=(a,b,c)
        s.edge[(a,b)]=t; s.edge[(b,c)]=t; s.edge[(c,a)]=t
        s.nbr[a].update((b,c)); s.nbr[b].update((a,c)); s.nbr[c].update((a,b))
    def rem(s,t):
        a,b,c=s.tri.pop(t); del s.edge[(a,b)],s.edge[(b,c)],s.edge[(c,a)]
    def opp(s,a,b):
        for v in s.tri[s.edge[(a,b)]]:
            if v!=a and v!=b: return v
    def split(s,a,b):
        """adds a vertex inside the triangle containing the oriented edge a->b"""
        t=s.edge[(a,b)]; c=s.opp(a,b); v=s.nv; s.nv+=1; s.nbr[v]=set(); s.rem(t)
        s.add(a,b,v); s.add(b,c,v); s.add(c,a,v); return v
    def flip(s,a,b,T):
        c=s.opp(a,b); d=s.opp(b,a)
        if c==d or d in s.nbr[c] or len(s.nbr[a])<=3 or len(s.nbr[b])<=3: return False
        dE=2*(len(s.nbr[c])+len(s.nbr[d])-len(s.nbr[a])-len(s.nbr[b]))+4
        if dE>0 and (T<=0 or random.random()>math.exp(-dE/T)): return False
        s.rem(s.edge[(a,b)]); s.rem(s.edge[(b,a)]); s.nbr[a].discard(b); s.nbr[b].discard(a)
        s.add(a,d,c); s.add(d,b,c); return True
    def bfs(s,r0,rmax=10**9):
        d={r0:0}; q=deque([r0])
        while q:
            u=q.popleft(); du=d[u]
            if du>=rmax: continue
            for w in s.nbr[u]:
                if w not in d: d[w]=du+1; q.append(w)
        return d
def dimensione(R,radici,r1,r2):
    """dimension from the growth of the number of vertices within distance r:  N(r) ~ r^d"""
    n1=n2=0
    for r0 in radici:
        d=R.bfs(r0,r2); n1+=sum(1 for x in d.values() if x<=r1); n2+=len(d)
    return math.log(n2/n1)/math.log(r2/r1)
def geometria(R,nrad,r):
    """measures that do not depend on which vertices are looked at: mean distance between random vertices, and volume of balls of radius r
    divided by that of flat space (1 = flat; the spread between centres tells whether there are swollen or pinched regions)"""
    vs=list(R.nbr); rad=random.sample(vs,nrad); dm=0; pal=[]; piatto=3*r*r+3*r+1
    for r0 in rad:
        d=R.bfs(r0); dm+=sum(d.values())/len(d); pal.append(sum(1 for x in d.values() if x<=r)/piatto)
    m=sum(pal)/len(pal); return dm/nrad,m,math.sqrt(sum((x-m)**2 for x in pal)/len(pal))
def esegui(n=40,crescita=0.5,Rf=100,T=0.3,seme=1,verbose=True,ovunque=False):
    random.seed(seme); R=Rete(n); V0=R.nv; centro=(n//2)*n+n//2
    d0c=R.bfs(centro); rz=n//5
    zona=[v for v,x in d0c.items() if x<=rz or ovunque]                       # active region: events happen here
    trac=[(i*n+j) for i in range(0,n,4) for j in range(0,n,4)]      # reference vertices (like distant galaxies)
    quieti=[v for v in trac if d0c[v]>=rz+4]; dentro=[v for v in trac if d0c[v]<=rz-1]
    lontani0=[v for v,x in d0c.items() if x>=rz+6]
    D0={v:R.bfs(v) for v in trac}; g0=geometria(R,60,n//5)
    dim0=dimensione(R,trac[::5],3,n//4)
    nadd=int(crescita*V0); nf=0; t0=time.time(); vs=list(R.nbr)
    for k in range(nadd):
        a=random.choice(zona); b=random.choice(tuple(R.nbr[a])); v=R.split(a,b); zona.append(v); vs.append(v)
        for _ in range(Rf):
            x=random.choice(vs); y=random.choice(tuple(R.nbr[x])); nf+=R.flip(x,y,T)
    D1={v:R.bfs(v) for v in trac}; g1=geometria(R,120,n//5)
    def rapporto(A,B,lo,hi):
        r=[D1[u][w]/D0[u][w] for u in A for w in B if u<w and lo<=D0[u][w]<=hi]
        if not r: return float('nan'),float('nan'),0
        m=sum(r)/len(r); return m,math.sqrt(sum((x-m)**2 for x in r)/len(r)),len(r)
    # spreading: how many NEW vertices there are around old vertices far from the active region (1 = spread everywhere)
    fr=[]
    for v in random.sample(lontani0,min(150,len(lontani0))):
        d=R.bfs(v,3); fr.append(sum(1 for w in d if w>=V0)/len(d))
    sparso=(sum(fr)/len(fr))/((R.nv-V0)/R.nv)
    gr=[len(x) for x in R.nbr.values()]; difetti=sum(1 for g in gr if g!=6)/len(gr)
    atteso=math.sqrt(R.nv/V0)
    res=dict(sparso=sparso,dist_media=g1[0]/g0[0],palle=g1[1],palle_scarto=g1[2],n=n,R=Rf,T=T,V0=V0,V1=R.nv,atteso=atteso,flip=nf,difetti=difetti,grado_max=max(gr),
             quieti=rapporto(quieti,quieti,4,n//3),zona=rapporto(dentro,dentro,2,n),tutti=rapporto(trac,trac,n//4,n//2),
             vicini=rapporto(quieti,quieti,4,7),lontani=rapporto(quieti,quieti,n//3,n),dim0=dim0,dim1=dimensione(R,quieti[::3],3,n//4),sec=time.time()-t0)
    if verbose:
        print('CREATION EVERYWHERE (comparison)' if ovunque else 'CREATION IN ONE REGION (13% of space)')
        print('network %dx%d, vertices %d -> %d (if created space spread everywhere, distances would grow by %.3f times)'%(n,n,V0,R.nv,atteso))
        print('  flips attempted per created vertex: %d (succeeded %d), temperature %.2f, time %.0f s'%(Rf,nf,T,res['sec']))
        print('  SPREADING of created space far from the active region: %.0f%% (0 = stays where it was born, 100 = spread everywhere)'%(100*sparso))
        print('  mean distance between any two points: x%.3f | volume of balls relative to flat space: %.2f, spread between regions %.0f%%'%(res['dist_media'],g1[1],100*g1[2]/g1[1]))
        print('  distances between references FAR from the active region: x%.3f (spread %.3f, %d pairs)'%res['quieti'])
        print('  distances between references INSIDE the active region:   x%.3f (spread %.3f, %d pairs)'%res['zona'])
        print('  Hubble law outside the region: near pairs x%.3f, far pairs x%.3f'%(res['vicini'][0],res['lontani'][0]))
        print('  dimension of space: before %.2f, after %.2f | irregular vertices %.0f%%, maximum number of neighbours %d'%(dim0,res['dim1'],100*difetti,max(gr)))
    return res
if __name__=='__main__':
    a=sys.argv[1:]
    esegui(int(a[0]) if a else 40,float(a[1]) if len(a)>1 else 0.5,int(a[2]) if len(a)>2 else 100,float(a[3]) if len(a)>3 else 0.3,int(a[4]) if len(a)>4 else 1,True,len(a)>5 and a[5]=='ovunque')
