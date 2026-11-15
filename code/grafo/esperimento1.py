"""Experiment 1 on the causal graph: link probability varying according to a source-sink budget.
Each new event links to each earlier event with probability p = p0/S, where S is the current 'space':
   S grows with activity (fuel that runs out) and falls with the accumulated mass (remnants of the activity).
Only the 'times' of events are simulated (longest chain of causes), which is enough to measure the spatial volume:
   the highest parent of a new event is the m-th event from the top, with m geometric with parameter p.
Measures: spatial volume per instant, comparison with the 0.6/p law, fluctuations, shape of the collapse."""
import math, random, sys, json
def run(N=5_000_000,W0=300.0,tau_frac=0.1,a_tot=3.0,seed=1,sink=True):
    rng=random.Random(seed); p0=0.6/W0; tau=N*tau_frac; a=a_tot/tau; g=(1+a_tot)/(tau*N) if sink else 0.0
    cnt=[1]; top=0; S=1.0; F=1.0; M=0.0; decay=1-1/tau
    rec=[]   # per level: (events, sum S, n)
    lev_events=[1]; lev_S=[1.0]; lev_n=[0]
    n=1; rnd=rng.random; log=math.log
    while n<2*N and S>0 and (sink or n<N):
        p=p0/S
        if p>=1: break
        m=int(log(1.0-rnd())/log(1.0-p))       # how many events, from the top, are skipped before the first parent
        L=top; c=cnt[L]
        while m>=c and L>0:
            m-=c; L-=1; c=cnt[L]
        l=L+1 if m<c else 0
        if l>top:
            cnt.append(0); lev_events.append(0); lev_S.append(0.0); lev_n.append(n); top=l
        cnt[l]+=1; lev_events[l]+=1; lev_S[l]+=S
        S+=a*F-g*M; M+=F; F*=decay; n+=1
    return dict(n=n,top=top,cnt=cnt,lev_S=lev_S,lev_n=lev_n,p0=p0,S_end=S)
def report(r,label,W0=300.0):
    cnt=r['cnt']; T=r['top']; print('=== %s: events %d, instants (emergent time) %d, maximum volume %d'%(label,r['n'],T,max(cnt)))
    k=20; step=T//k
    print('   fraction of time | measured volume | expected volume 0.6/p | ratio | relative fluctuation x sqrt(volume)')
    for i in range(k):
        a,b=i*step,(i+1)*step; w=cnt[a:b]; mean=sum(w)/len(w)
        Sm=sum(r['lev_S'][a:b])/max(1,sum(w)); exp=W0*Sm
        var=sum((x-mean)**2 for x in w)/len(w)
        # remove the rough linear trend: use first differences
        d=[w[j+1]-w[j] for j in range(len(w)-1)]; sd=math.sqrt(sum(x*x for x in d)/len(d)/2)
        print('        %.2f-%.2f      %8.1f          %8.1f         %.2f        %.2f'%(i/k,(i+1)/k,mean,exp,mean/exp,sd/math.sqrt(mean)))
    imax=max(range(T),key=lambda i:cnt[i]); half=max(cnt)/2
    t_up=next(i for i in range(T) if cnt[i]>=half and i>5); t_dn=max(i for i in range(T) if cnt[i]>=half)
    print('   maximum at %.0f%% of the emergent time; half height: rise at %.0f%%, fall at %.0f%%'%(100*imax/T,100*t_up/T,100*t_dn/T))
    if r['S_end']<=0 or r['p0']/max(r['S_end'],1e-12)>=1:
        m=max(cnt); ih=max(i for i in range(T) if cnt[i]>=m/10); print('   full collapse: space runs out after %d events; from 10%% of the maximum to the end there are %d instants out of %d (%.1f%%)'%(r['n'],T-ih,T,100*(T-ih)/T))
        print('   volume in the last 300 instants, in blocks of 30:',[round(sum(cnt[i:i+30])/30,1) for i in range(T-300,T,30)])
    else: print('   no collapse: the run stops at %d events with volume still full (the last instants are incomplete and should be ignored)'%r['n'])
r=run(); report(r,'source and sink')
r2=run(sink=False); report(r2,'source only (no sink)')
json.dump(dict(cnt=r['cnt'][::max(1,r['top']//400)]),open('esperimento1.json','w'))
