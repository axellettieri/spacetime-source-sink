"""Esperimento 2: lo spazio nasce dagli eventi, direttamente sul grafo (nessuna equazione di espansione imposta).
Universo = anello di celle (una dimensione di spazio) che evolve a passi; il grafo causale e' l'anello nel tempo.
- N oggetti (particelle) si muovono a caso: la probabilita' di salto e' la 'temperatura' divisa per la massa.
- Quando due oggetti si toccano c'e' un EVENTO: crea K celle nuove fra i due e li raffredda (fattore C).
  Con probabilita' PM i due si fondono (nasce una massa piu' grande).
- Ogni oggetto ASSORBE le celle accanto a se', al ritmo G per unita' di massa.
Si misura: dimensione dell'universo nel tempo, ritmo degli eventi, e quanto l'espansione varia da zona a zona.
Uso: python esperimento2.py N0 PASSI K C PM G [seme]"""
import numpy as np, sys, math, json
def run(N0=200000,steps=20000,K=4,C=0.98,PM=0.02,G=0.0,seed=1,verbose=True,snap=None):
    rng=np.random.default_rng(seed)
    gap=np.ones(N0,np.int64); m=np.ones(N0,np.int64); T=np.ones(N0); ident=np.arange(N0)
    hist=[]; ev_tot=0; snaps={}
    for t in range(steps+1):
        N=len(gap)
        if snap and t in snap: snaps[t]=(ident.copy(),gap.copy(),m.copy())
        if t%max(1,steps//400)==0 or N<=2:
            W=int(gap.sum())+N; hist.append((t,W,N,ev_tot,float((T*m).sum()/m.sum()),int(m.max())))
            if N<=2 or gap.sum()==0: break
        # 1) moto: due mezzi passi (oggetti pari, poi dispari) per evitare conflitti
        h=np.minimum(0.5,0.5*T/m)
        for par in (0,1):
            idx=np.arange(par,N-(N%2 if par==0 else 0),2)
            if par==1 and N%2==1: idx=idx[:-1] if len(idx) and idx[-1]==N-1 else idx
            r=rng.random(len(idx)); hh=h[idx]
            right=(r<hh)&(gap[idx]>0); left=(r>1-hh)&(gap[idx-1]>0)
            gap[idx[right]]-=1; gap[idx[right]-1]+=1
            gap[idx[left]]+=1; gap[idx[left]-1]-=1
        # 2) eventi: coppie a contatto (alternando la parita' per non avere catene)
        c=np.nonzero(gap[:-1]==0)[0]; c=c[c%2==t%2]
        if len(c):
            ev_tot+=len(c); gap[c]+=K; T[c]*=C; T[c+1]*=C
            mg=c[rng.random(len(c))<PM]
            if len(mg):
                T[mg]=(T[mg]*m[mg]+T[mg+1]*m[mg+1])/(m[mg]+m[mg+1]); m[mg]+=m[mg+1]; gap[mg]+=gap[mg+1]
                keep=np.ones(N,bool); keep[mg+1]=False
                gap=gap[keep]; m=m[keep]; T=T[keep]; ident=ident[keep]
        # 3) assorbimento: ogni oggetto toglie celle ai due lati
        if G>0:
            d=rng.poisson(G*m/2.0); dem=d+np.roll(rng.poisson(G*m/2.0),-1)
            gap=np.maximum(gap-dem,0)
    return np.array(hist,float),snaps
def local_scatter(s1,s2,dt,sizes=(10,100,1000,10000)):
    """variazione dell'espansione da zona a zona: finestre di n oggetti consecutivi, stessi oggetti di confine ai due istanti"""
    id1,g1,m1=s1; id2,g2,m2=s2
    x1=np.concatenate([[0],np.cumsum(g1+1)[:-1]]); x2=np.concatenate([[0],np.cumsum(g2+1)[:-1]])
    pos2=dict(zip(id2.tolist(),x2.tolist())); out={}
    for n in sizes:
        b=np.arange(0,len(id1)-n,n); H=[];dens=[]
        for i in b:
            a,bb=int(id1[i]),int(id1[i+n])
            if a in pos2 and bb in pos2 and pos2[bb]>pos2[a]:
                L1=x1[i+n]-x1[i]; H.append(math.log((pos2[bb]-pos2[a])/L1)/dt); dens.append(n/L1)
        H=np.array(H); dens=np.array(dens)
        if len(H)>5: out[n]=(float(H.mean()),float(H.std()/abs(H.mean())),float(np.corrcoef(H,dens)[0,1]),len(H))
    return out
if __name__=='__main__':
    a=sys.argv[1:]
    N0=int(a[0]) if len(a)>0 else 200000; steps=int(a[1]) if len(a)>1 else 20000
    K=int(a[2]) if len(a)>2 else 4; C=float(a[3]) if len(a)>3 else 0.98; PM=float(a[4]) if len(a)>4 else 0.02
    G=float(a[5]) if len(a)>5 else 0.0; seed=int(a[6]) if len(a)>6 else 1
    t1,t2=steps//20,steps//10
    h,sn=run(N0,steps,K,C,PM,G,seed,snap=(t1,t2))
    i=int(np.argmax(h[:,1])); print('N0=%d K=%d C=%.3f PM=%.3f G=%.4f'%(N0,K,C,PM,G))
    print('dimensione: inizio %d -> massimo %d (x%.1f) al passo %d -> fine %d al passo %d | oggetti finali %d, massa massima %d, eventi %d'%(h[0,1],h[i,1],h[i,1]/h[0,1],h[i,0],h[-1,1],h[-1,0],h[-1,2],h[-1,5],h[-1,3]))
    q=np.linspace(0,len(h)-1,9).astype(int)
    print('passo:dimensione/iniziale  '+'  '.join('%d:%.2f'%(h[j,0],h[j,1]/h[0,1]) for j in q))
    if t1 in sn and t2 in sn:
        for n,(Hm,sc,cc,nw) in local_scatter(sn[t1],sn[t2],t2-t1).items():
            print('zone di %5d oggetti: espansione media %.2e per passo, variazione fra zone %.1f%%, correlazione con la densita\' %+.2f (%d zone)'%(n,Hm,100*sc,cc,nw))
    np.save('esp2_N%d_C%.3f_G%.4f.npy'%(N0,C,G),h)
