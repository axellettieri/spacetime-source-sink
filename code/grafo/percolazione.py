"""Esperimento 0 sul grafo causale: crescita sequenziale con la regola piu' semplice (percolazione transitiva).
Ogni nuovo evento sceglie come cause dirette gli eventi precedenti, ciascuno con probabilita' p = 1/2^k; le cause delle cause sono cause.
'tempo' di un evento = lunghezza della catena di cause piu' lunga che lo raggiunge;
'volume spaziale' a un dato tempo = numero di eventi con quel tempo;
'post' = evento che ha per cause TUTTI i precedenti e da cui discendono TUTTI i successivi: l'universo passa per un solo evento."""
import random, sys, statistics
from collections import Counter
def run(N,k,seed=1):
    rng=random.Random(seed); anc=[0]*N; lev=[0]*N
    for n in range(1,N):
        ch=rng.getrandbits(n)
        for _ in range(k-1): ch&=rng.getrandbits(n)
        a=0; l=0
        while ch:
            j=ch.bit_length()-1; a|=anc[j]|(1<<j)
            if lev[j]+1>l: l=lev[j]+1
            ch&=~a
        anc[n]=a; lev[n]=l
    common=(1<<N)-1; posts=[]
    for n in range(N-1,-1,-1):
        if anc[n]==(1<<n)-1 and (common>>n)&1 and n<N-50: posts.append(n)
        common&=anc[n]|(1<<n)
    posts.reverse(); return lev,posts
N=int(float(sys.argv[1])) if len(sys.argv)>1 else 40000
for k in (1,2,3):
    lev,posts=run(N,k); p=0.5**k
    gaps=[b-a for a,b in zip(posts[:-1],posts[1:])]
    print('p=%.3f  eventi=%d  post=%d  eventi per ciclo: mediana %s, medio %s, massimo %s'%(p,N,len(posts),statistics.median(gaps) if gaps else '-',round(statistics.mean(gaps),1) if gaps else '-',max(gaps) if gaps else '-'))
    c=Counter(lev); w=[c[t] for t in range(max(c)+1)]; print('   volume spaziale medio per istante: %.2f  (massimo %d)'%(statistics.mean(w),max(w)))
    if gaps:
        i=gaps.index(max(gaps)); a,b=posts[i],posts[i+1]
        cc=Counter(lev[n]-lev[a] for n in range(a,b+1)); T=max(cc); ww=[cc[t] for t in range(T+1)]
        s=max(1,len(ww)//10); print('   ciclo piu\' lungo: %d istanti; volume lungo il ciclo:'%T,[round(sum(ww[i:i+s])/len(ww[i:i+s]),1) for i in range(0,len(ww),s)][:11])
