"""Experiment 6: the maximum speed on the graph.
A signal advances one link per step (it cannot do more: every local rule has this limit).
PART A - is the limit the same in all directions and at all points?
   A regular network (all triangles equal) is compared with a random one (random points joined by triangles).
   Checking requires the point coordinates, used ONLY to measure, never by the signal.
PART B - what happens to the signal if space is born along the way?
   Closed network where vertices are created everywhere at a constant rate (with relaxation). We measure when the signal reaches
   targets at various distances and compare with the law of light in an expanding universe.
Usage: python esperimento6.py A      or      python esperimento6.py B [rate] [seed] [side] [steps]"""
import sys, math, random, time
import numpy as np
from collections import deque
def bfs(nbr,r0):
    d={r0:0}; q=deque([r0])
    while q:
        u=q.popleft()
        for w in nbr[u]:
            if w not in d: d[w]=d[u]+1; q.append(w)
    return d
def analizza(nome,nbr,spost,passo,nsorg=40,lo=15,hi=25,nb=36):
    """spost(a,b) = vector from a to b (minimum image); passo = mean distance between neighbours (unit of length)"""
    N=len(nbr); bins=[[] for _ in range(nb)]; per_s=[]
    for s in random.sample(range(N),nsorg):
        d=bfs(nbr,s); vs=[]
        for t,h in d.items():
            if h==0: continue
            dx,dy=spost(s,t); de=math.hypot(dx,dy)/passo
            if lo<de<hi:
                v=de/h; vs.append(v); bins[int(((math.atan2(dy,dx)+math.pi)/(2*math.pi))*nb)%nb].append(v)
        per_s.append(sum(vs)/len(vs))
    mb=[sum(b)/len(b) for b in bins]; m=sum(mb)/nb
    print('%-16s mean speed %.3f | between directions: from %.3f to %.3f (difference %.1f%%) | between starting points: spread %.1f%%'%(nome,m,min(mb),max(mb),100*(max(mb)-min(mb))/m,100*np.std(per_s)/np.mean(per_s)))
    return mb
def parte_A(seme=1):
    random.seed(seme); rng=np.random.default_rng(seme)
    # regular network
    n=70; idx=lambda i,j:(i%n)*n+(j%n); nbr={idx(i,j):[idx(i+1,j),idx(i-1,j),idx(i,j+1),idx(i,j-1),idx(i+1,j-1),idx(i-1,j+1)] for i in range(n) for j in range(n)}
    def sp(a,b):
        di=(b//n-a//n+n//2)%n-n//2; dj=(b%n-a%n+n//2)%n-n//2; best=None
        for oi in (-n,0,n):
            for oj in (-n,0,n):
                x=(di+oi)+(dj+oj)/2; y=(dj+oj)*math.sqrt(3)/2
                if best is None or x*x+y*y<best[0]: best=(x*x+y*y,x,y)
        return best[1],best[2]
    analizza('rete regolare',nbr,sp,1.0)
    # random network: random points, Delaunay triangulation on a torus
    from scipy.spatial import Delaunay
    N=5000; P=rng.random((N,2)); tiles=np.concatenate([P+np.array([ox,oy]) for ox in (-1,0,1) for oy in (-1,0,1)]); T=Delaunay(tiles)
    nb={i:set() for i in range(N)}
    for tri in T.simplices:
        for a in range(3):
            u,w=tri[a],tri[(a+1)%3]
            if 4*N<=u<5*N or 4*N<=w<5*N:
                nb[u%N].add(w%N); nb[w%N].add(u%N)
    for i in nb: nb[i].discard(i)
    lung=np.mean([math.hypot(*(((P[w]-P[u]+0.5)%1)-0.5)) for u in range(0,N,10) for w in nb[u]])
    sp2=lambda a,b: tuple(((P[b]-P[a]+0.5)%1)-0.5)
    analizza('rete casuale',nb,sp2,lung)
def parte_B(H=0.01,seme=1,n=50,passi=140,RF=30,TG=0.3):
    sys.path.insert(0,'.'); from esperimento3 import Rete
    random.seed(seme); R=Rete(n); V0=R.nv; sorg=(n//2)*n+n//2; d0=R.bfs(sorg)
    bers={}
    for v,x in d0.items():
        if 3<=x<=n//2-1 and x%4==0 and v%5==0: bers.setdefault(x,[]).append(v)
    trac=random.sample(range(V0),40); D0=[R.bfs(v) for v in trac[:12]]
    coppie=[(i,w) for i in range(12) for w in trac if 6<=D0[i][w]<=n//3]
    raggiunti={sorg}; fronte={sorg}; arrivo={}; scala=[1.0]; tutti=list(R.nbr); t0=time.time(); acc=0.0
    for t in range(1,passi+1):
        nuovo=set()
        for u in fronte:
            for w in R.nbr[u]:
                if w not in raggiunti: nuovo.add(w)
        raggiunti|=nuovo
        fronte={u for u in (fronte|nuovo) if any(w not in raggiunti for w in R.nbr[u])}
        for x,l in bers.items():
            for v in l:
                if v in raggiunti and v not in arrivo: arrivo[v]=t
        acc+=2*H*R.nv; k=int(acc); acc-=k                          # space (number of vertices) grows by 2H per step: lengths by H
        for _ in range(k):
            a=random.choice(tutti); b=random.choice(tuple(R.nbr[a])); c=R.opp(a,b); v=R.split(a,b); tutti.append(v)
            if a in raggiunti and b in raggiunti and c in raggiunti: raggiunti.add(v)
            loc=[v,a,b,c]
            for _ in range(RF):
                x=random.choice(loc); y=random.choice(tuple(R.nbr[x])); R.flip(x,y,TG)
        if t%10==0 or t==passi:
            D1=[R.bfs(trac[i]) for i in range(12)]; scala.append(sum(D1[i][w]/D0[i][w] for i,w in coppie)/len(coppie))
        else: scala.append(None)
    # effective scale factor, interpolated
    ts=[i for i,s in enumerate(scala) if s is not None]; a_eff=np.interp(range(passi+1),ts,[scala[i] for i in ts])
    eta=np.concatenate([[0],np.cumsum(1/a_eff[1:])])                 # initial distance that can be covered by step t
    print('creation rate H=%.3f per step, vertices %d -> %d, effective lengths x%.2f (expected if uniform: x%.2f)  [%.0f s]'%(H,V0,R.nv,a_eff[-1],math.sqrt(R.nv/V0),time.time()-t0))
    print('horizon: in %d steps light covers at most %.1f steps of initial distance; without creation it would cover %d'%(passi,eta[-1],passi))
    print(' initial distance | measured arrival (mean) | arrival predicted by the expansion law | without space creation')
    for x in sorted(bers):
        ar=[arrivo[v] for v in bers[x] if v in arrivo]; frac=len(ar)/len(bers[x])
        prev=int(np.searchsorted(eta,x)) if x<=eta[-1] else None
        print('      %3d         |  %s  |  %s  |  %d'%(x,('%6.1f (reached %3.0f%%)'%(sum(ar)/len(ar),100*frac)) if ar else '  never reached      ',('%4d'%prev) if prev else 'never (beyond the horizon)',x))
if __name__=='__main__':
    a=sys.argv[1:]
    if not a or a[0]=='A': parte_A()
    if not a or a[0]=="B": parte_B(float(a[1]) if len(a)>1 else 0.02,int(a[2]) if len(a)>2 else 1,int(a[3]) if len(a)>3 else 120,int(a[4]) if len(a)>4 else 70)
