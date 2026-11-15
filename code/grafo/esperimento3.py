"""Esperimento 3: lo spazio come grafo senza coordinate (due dimensioni di spazio).
Lo spazio e' una rete di triangoli chiusa su se stessa: ogni vertice e' un atomo di spazio, le distanze si contano in passi sul grafo.
- CREAZIONE locale: in una zona attiva (dove avvengono gli eventi) si aggiungono vertici, uno alla volta, spezzando un triangolo in tre.
- RILASSAMENTO: scambi di lati (flip) accettati se avvicinano la rete a quella regolare (6 vicini per vertice = spazio piatto),
  con un po' di agitazione termica T. R = tentativi di scambio per ogni vertice creato.
Domande: lo spazio creato in un punto resta li' o si sparge? Le distanze crescono davvero? Lo spazio resta a due dimensioni?
Uso: python esperimento3.py n crescita R T [seme] [ovunque]     (n = lato della rete iniziale, crescita = frazione di vertici aggiunti)"""
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
        """aggiunge un vertice dentro il triangolo che contiene il lato orientato a->b"""
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
    """dimensione dalla crescita del numero di vertici entro distanza r:  N(r) ~ r^d"""
    n1=n2=0
    for r0 in radici:
        d=R.bfs(r0,r2); n1+=sum(1 for x in d.values() if x<=r1); n2+=len(d)
    return math.log(n2/n1)/math.log(r2/r1)
def geometria(R,nrad,r):
    """misure che non dipendono da quali vertici si guardano: distanza media fra vertici a caso, e volume delle palle di raggio r
    diviso quello dello spazio piatto (1 = piatto; lo scarto fra radici dice se ci sono zone gonfie o strozzate)"""
    vs=list(R.nbr); rad=random.sample(vs,nrad); dm=0; pal=[]; piatto=3*r*r+3*r+1
    for r0 in rad:
        d=R.bfs(r0); dm+=sum(d.values())/len(d); pal.append(sum(1 for x in d.values() if x<=r)/piatto)
    m=sum(pal)/len(pal); return dm/nrad,m,math.sqrt(sum((x-m)**2 for x in pal)/len(pal))
def esegui(n=40,crescita=0.5,Rf=100,T=0.3,seme=1,verbose=True,ovunque=False):
    random.seed(seme); R=Rete(n); V0=R.nv; centro=(n//2)*n+n//2
    d0c=R.bfs(centro); rz=n//5
    zona=[v for v,x in d0c.items() if x<=rz or ovunque]                       # zona attiva: qui avvengono gli eventi
    trac=[(i*n+j) for i in range(0,n,4) for j in range(0,n,4)]      # vertici di riferimento (come galassie lontane)
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
    # sparpagliamento: quanti vertici NUOVI ci sono attorno ai vertici vecchi che stavano lontano dalla zona attiva (1 = sparsi ovunque)
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
        print('CREAZIONE OVUNQUE (confronto)' if ovunque else 'CREAZIONE IN UNA ZONA (13% dello spazio)')
        print('rete %dx%d, vertici %d -> %d (se lo spazio creato si spargesse ovunque le distanze crescerebbero di %.3f volte)'%(n,n,V0,R.nv,atteso))
        print('  scambi tentati per vertice creato: %d (riusciti %d), temperatura %.2f, tempo %.0f s'%(Rf,nf,T,res['sec']))
        print('  SPARPAGLIAMENTO dello spazio creato lontano dalla zona attiva: %.0f%% (0 = resta dove e\' nato, 100 = sparso ovunque)'%(100*sparso))
        print('  distanza media fra due punti qualsiasi: x%.3f | volume delle palle rispetto allo spazio piatto: %.2f, scarto fra zone %.0f%%'%(res['dist_media'],g1[1],100*g1[2]/g1[1]))
        print('  distanze fra riferimenti LONTANI dalla zona attiva: x%.3f (scarto %.3f, %d coppie)'%res['quieti'])
        print('  distanze fra riferimenti DENTRO la zona attiva:      x%.3f (scarto %.3f, %d coppie)'%res['zona'])
        print('  legge di Hubble fuori zona: coppie vicine x%.3f, coppie lontane x%.3f'%(res['vicini'][0],res['lontani'][0]))
        print('  dimensione dello spazio: prima %.2f, dopo %.2f | vertici non regolari %.0f%%, massimo numero di vicini %d'%(dim0,res['dim1'],100*difetti,max(gr)))
    return res
if __name__=='__main__':
    a=sys.argv[1:]
    esegui(int(a[0]) if a else 40,float(a[1]) if len(a)>1 else 0.5,int(a[2]) if len(a)>2 else 100,float(a[3]) if len(a)>3 else 0.3,int(a[4]) if len(a)>4 else 1,True,len(a)>5 and a[5]=='ovunque')
