"""Esperimento 1 sul grafo causale: probabilita' di legame variabile secondo un bilancio sorgente-pozzo.
Ogni nuovo evento si lega a ciascun evento precedente con probabilita' p = p0/S, dove S e' lo 'spazio' corrente:
   S aumenta con l'attivita' (carburante che si esaurisce) e cala con la massa accumulata (resti dell'attivita').
Si simulano solo i 'tempi' degli eventi (catena di cause piu' lunga), che bastano per misurare il volume spaziale:
   il genitore piu' alto di un nuovo evento e' l'm-esimo evento dall'alto, con m geometrico di parametro p.
Misure: volume spaziale per istante, confronto con la legge 0,6/p, fluttuazioni, forma del collasso."""
import math, random, sys, json
def run(N=5_000_000,W0=300.0,tau_frac=0.1,a_tot=3.0,seed=1,sink=True):
    rng=random.Random(seed); p0=0.6/W0; tau=N*tau_frac; a=a_tot/tau; g=(1+a_tot)/(tau*N) if sink else 0.0
    cnt=[1]; top=0; S=1.0; F=1.0; M=0.0; decay=1-1/tau
    rec=[]   # per livello: (eventi, somma S, n)
    lev_events=[1]; lev_S=[1.0]; lev_n=[0]
    n=1; rnd=rng.random; log=math.log
    while n<2*N and S>0 and (sink or n<N):
        p=p0/S
        if p>=1: break
        m=int(log(1.0-rnd())/log(1.0-p))       # quanti eventi, dall'alto, vengono saltati prima del primo genitore
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
    cnt=r['cnt']; T=r['top']; print('=== %s: eventi %d, istanti (tempo emergente) %d, volume massimo %d'%(label,r['n'],T,max(cnt)))
    k=20; step=T//k
    print('   frazione del tempo | volume misurato | volume atteso 0,6/p | rapporto | fluttuazione relativa x radice(volume)')
    for i in range(k):
        a,b=i*step,(i+1)*step; w=cnt[a:b]; mean=sum(w)/len(w)
        Sm=sum(r['lev_S'][a:b])/max(1,sum(w)); exp=W0*Sm
        var=sum((x-mean)**2 for x in w)/len(w)
        # togli la tendenza lineare grossolana: usa differenze prime
        d=[w[j+1]-w[j] for j in range(len(w)-1)]; sd=math.sqrt(sum(x*x for x in d)/len(d)/2)
        print('        %.2f-%.2f      %8.1f          %8.1f         %.2f        %.2f'%(i/k,(i+1)/k,mean,exp,mean/exp,sd/math.sqrt(mean)))
    imax=max(range(T),key=lambda i:cnt[i]); half=max(cnt)/2
    t_up=next(i for i in range(T) if cnt[i]>=half and i>5); t_dn=max(i for i in range(T) if cnt[i]>=half)
    print('   massimo al %.0f%% del tempo emergente; meta\' altezza: salita al %.0f%%, discesa al %.0f%%'%(100*imax/T,100*t_up/T,100*t_dn/T))
    if r['S_end']<=0 or r['p0']/max(r['S_end'],1e-12)>=1:
        m=max(cnt); ih=max(i for i in range(T) if cnt[i]>=m/10); print('   collasso completo: lo spazio si esaurisce dopo %d eventi; dal 10%% del massimo alla fine passano %d istanti su %d (%.1f%%)'%(r['n'],T-ih,T,100*(T-ih)/T))
        print('   volume negli ultimi 300 istanti, a blocchi di 30:',[round(sum(cnt[i:i+30])/30,1) for i in range(T-300,T,30)])
    else: print('   nessun collasso: la corsa si ferma a %d eventi con volume ancora pieno (gli ultimi istanti sono incompleti e vanno ignorati)'%r['n'])
r=run(); report(r,'sorgente e pozzo')
r2=run(sink=False); report(r2,'solo sorgente (nessun pozzo)')
json.dump(dict(cnt=r['cnt'][::max(1,r['top']//400)]),open('esperimento1.json','w'))
