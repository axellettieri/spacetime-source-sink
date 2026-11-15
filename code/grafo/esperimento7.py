"""Esperimento 7: il tempo di un oggetto in moto, contato in eventi.
Idea: un oggetto e' un brulicare di eventi interni. Qui l'oggetto piu' semplice possibile: due "specchi" che si scambiano lampi.
Ogni lampo si allarga sul grafo di un legame per passo, in tutte le direzioni (nessuna coordinata usata dal segnale).
Un TIC dell'orologio = un lampo va dallo specchio A allo specchio B e un altro torna indietro: sono gli eventi interni dell'oggetto.
Lo sperimentatore sposta l'oggetto a velocita' costante (qui si usano le coordinate, solo per spostarlo).
Domanda: un oggetto in moto compie meno eventi interni nello stesso numero di passi del grafo? E quanto meno?
La relativita' prevede: orologio messo di traverso al moto -> rallenta di 1/sqrt(1-v^2/c^2).
Uso: python esperimento7.py [N punti] [tic per misura] [seme] [distanza fra gli specchi]"""
import sys, math, time
import numpy as np
from collections import deque
from scipy.spatial import Delaunay, cKDTree
def rete_casuale(N,rng):
    P=rng.random((N,2)); tiles=np.concatenate([P+np.array([ox,oy]) for ox in (-1,0,1) for oy in (-1,0,1)]); T=Delaunay(tiles)
    nb=[set() for _ in range(N)]
    for tri in T.simplices:
        for a in range(3):
            u,w=int(tri[a]),int(tri[(a+1)%3])
            if 4*N<=u<5*N or 4*N<=w<5*N: nb[u%N].add(w%N); nb[w%N].add(u%N)
    return P,[sorted(s-{i}) for i,s in enumerate(nb)]
def bfs(nbr,r0,N,cap=10**9):
    d=np.full(N,10**8,np.int32); d[r0]=0; q=deque([r0])
    while q:
        u=q.popleft(); du=d[u]+1
        if du>cap: break
        for w in nbr[u]:
            if d[w]>du: d[w]=du; q.append(w)
    return d
def orologio(P,nbr,tree,rng,LE,vE,traverso,ntic,cap=10**9):
    """ritorna il numero medio di passi per tic. vE = velocita' dell'oggetto (lunghezza per passo), LE = distanza fra gli specchi."""
    N=len(P); tot=0
    for k in range(ntic):
        th=rng.random()*2*math.pi; dirv=np.array([math.cos(th),math.sin(th)]); perp=np.array([-dirv[1],dirv[0]])
        off=(perp if traverso else dirv)*LE; xA=rng.random(2); t=0
        nodo=lambda x: int(tree.query(x%1.0)[1])
        for verso in (0,1):                                   # andata A->B, ritorno B->A
            src=xA+dirv*vE*t+(off if verso else 0); d=bfs(nbr,nodo(src),N,cap); s=0
            while True:
                s+=1; dst=xA+dirv*vE*(t+s)+(0 if verso else off)
                if d[nodo(dst)]<=s: break
            t+=s
        tot+=t
    return tot/ntic
if __name__=='__main__':
    a=sys.argv[1:]; N=int(a[0]) if a else 120000; ntic=int(a[1]) if len(a)>1 else 100; rng=np.random.default_rng(int(a[2]) if len(a)>2 else 1)
    t0=time.time(); P,nbr=rete_casuale(N,rng); tree=cKDTree(P,boxsize=1.0); LE=float(a[3]) if len(a)>3 else 0.06
    T0=orologio(P,nbr,tree,rng,LE,0.0,True,2*ntic); c=2*LE/T0
    print('rete casuale di %d punti; specchi a distanza %.1f legami; velocita\' massima c = %.5f per passo; orologio fermo: %.2f passi per tic'%(N,T0/2,c,T0))
    print(' velocita\' v/c | orologio di TRAVERSO: rallentamento misurato | previsto dalla relativita\' | orologio LUNGO il moto: misurato | previsto senza contrazione | previsto con contrazione')
    for b in (0.2,0.4,0.6,0.8,0.9):
        g=1/math.sqrt(1-b*b); Tt=orologio(P,nbr,tree,rng,LE,b*c,True,ntic,int(T0*g)+10)
        Tl=orologio(P,nbr,tree,rng,LE,b*c,False,ntic,int(T0/(1-b))+10) if LE/(1-b)<0.45 else float('nan')     # oltre: il lampo farebbe il giro dell'universo chiuso
        print('     %.1f      |              %.3f                         |          %.3f              |            %.3f              |          %.3f            |        %.3f'%(b,Tt/T0,g,Tl/T0,g*g,g),flush=True)
    print('tempo %.0f s'%(time.time()-t0))
