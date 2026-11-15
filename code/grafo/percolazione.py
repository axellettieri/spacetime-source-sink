"""Experiment 0 on the causal graph: sequential growth with the simplest rule (transitive percolation).
Each new event picks earlier events as direct causes, each with probability p = 1/2^k; causes of causes are causes.
'time' of an event = length of the longest chain of causes reaching it;
'spatial volume' at a given time = number of events with that time;
'post' = event whose causes are ALL the earlier ones and from which ALL later ones descend: the universe passes through a single event."""
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
    print('p=%.3f  events=%d  post=%d  events per cycle: median %s, mean %s, maximum %s'%(p,N,len(posts),statistics.median(gaps) if gaps else '-',round(statistics.mean(gaps),1) if gaps else '-',max(gaps) if gaps else '-'))
    c=Counter(lev); w=[c[t] for t in range(max(c)+1)]; print('   mean spatial volume per instant: %.2f  (maximum %d)'%(statistics.mean(w),max(w)))
    if gaps:
        i=gaps.index(max(gaps)); a,b=posts[i],posts[i+1]
        cc=Counter(lev[n]-lev[a] for n in range(a,b+1)); T=max(cc); ww=[cc[t] for t in range(T+1)]
        s=max(1,len(ww)//10); print('   longest cycle: %d instants; volume along the cycle:'%T,[round(sum(ww[i:i+s])/len(ww[i:i+s]),1) for i in range(0,len(ww),s)][:11])
