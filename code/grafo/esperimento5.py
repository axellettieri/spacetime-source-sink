"""Esperimento 5: come l'esperimento 4, piu' il DECADIMENTO: un oggetto composto (massa >= 2) che torna abbastanza caldo
(energia per unita' di massa >= TC) e' instabile e si spezza in due, con probabilita' PF a ogni passo.
All'inizio le particelle singole non possono spezzarsi; durante l'espansione gli oggetti fusi si raffreddano e sono stabili;
durante il collasso riassorbono lo spazio, si riscaldano e, superata la soglia, decadono: nuovi oggetti, nuovi eventi, nuovo spazio.
La corsa non si ferma: si guarda se l'universo oscilla (cicli) e se i cicli si somigliano.
--- regole ereditate dall'esperimento 4 ---
Spazio = rete di triangoli chiusa (superficie di una sfera), si parte dal grafo piu' piccolo possibile: 4 vertici ("singolarita'").
- NP particelle calde, tutte nei 4 vertici iniziali. Ogni particella salta a un vertice vicino con probabilita' ~ energia/massa.
- EVENTO: due oggetti nello stesso vertice. Perdono una frazione (1-C) della loro energia termica; quell'energia diventa SPAZIO:
  un vertice nuovo ogni EPS di energia. Con probabilita' PM i due si fondono in una massa unica.
- ASSORBIMENTO: ogni oggetto toglie vertici vicini al ritmo G per unita' di massa e ne riprende l'energia (EPS per vertice).
  La materia che stava nel vertice tolto cade sull'oggetto.
- RILASSAMENTO: scambi di lati attorno ai vertici toccati, verso la rete regolare.
Nulla si crea e nulla si distrugge: massa totale ed energia totale (termica + spazio) sono conservate esattamente.
Uso: python esperimento5.py NP EPS C PM G passi TC PF [seme]     [seme MCRIT TCRIT]   (TC=0: decadimento solo alla singolarita'; TC=99: mai;
     MCRIT: frazione della massa totale oltre la quale un oggetto decade del tutto; TCRIT: energia minima per unita' di massa)"""
import random, math, sys, time
from esperimento3 import Rete
from esperimento4 import Sfera, dimensione
def esegui(NP=300,EPS=0.02,C=0.9,PM=0.05,G=0.005,passi=400000,TC=0.9,PF=0.01,seme=1,MCRIT=0.0,TCRIT=0.0,RF=20,TG=0.3,verbose=True):
    random.seed(seme); R=Sfera()
    pos=[random.randrange(4) for _ in range(NP)]; m=[1]*NP; e=[1.0]*NP; vivo=[True]*NP
    occ={v:[] for v in range(4)}
    for i,v in enumerate(pos): occ[v].append(i)
    liberi=[]; nfr=0; critici=[]; vmin_dopo=None
    resto=0.0; Etot0=sum(e); storia=[]; ev=0; t0=time.time(); vmax=0; dim_max=float('nan')
    def rilassa(vv):
        for v in vv:
            if v not in R.nbr: continue
            loc=[v]+list(R.nbr[v])
            for _ in range(RF):
                x=random.choice(loc); y=random.choice(tuple(R.nbr[x])); R.flip(x,y,TG)
    for t in range(passi+1):
        V=len(R.nbr); ogg=[i for i in range(len(vivo)) if vivo[i]]
        if t%max(1,passi//2000)==0:
            espazio=(V-4)*EPS+resto; storia.append((t,V,len(ogg),max(m[i] for i in ogg),sum(e[i] for i in ogg),espazio,ev))
        if V>vmax: vmax=V
        if verbose and t%max(1,passi//20)==0:
            print('passo %6d: vertici %6d, oggetti %4d, massa massima %4d, energia termica %.1f + spazio %.1f = %.1f, eventi %d, dimensione %.2f  [%.0f s]'%(t,V,len(ogg),max(m[i] for i in ogg),sum(e[i] for i in ogg),(V-4)*EPS+resto,sum(e[i] for i in ogg)+(V-4)*EPS+resto,ev,dimensione(R),time.time()-t0),flush=True)
        # 1) moto
        for i in ogg:
            if random.random()<min(0.5,0.5*e[i]/m[i]):
                v=pos[i]; w=random.choice(tuple(R.nbr[v])); occ[v].remove(i); occ.setdefault(w,[]).append(i); pos[i]=w
        # 2) eventi
        toccati=[]
        for v in [v for v,l in occ.items() if len(l)>=2]:
            l=occ[v]; i,j=l[0],l[1]; ev+=1
            dE=(1-C)*(e[i]+e[j]); e[i]*=C; e[j]*=C; resto+=dE
            n=int(resto/EPS); resto-=n*EPS; nuovi=[]
            for _ in range(n):
                nuovi.append(R.split(v,random.choice(tuple(R.nbr[v]))))
            for w in nuovi: occ[w]=[]
            if random.random()<PM:
                m[i]+=m[j]; e[i]+=e[j]; vivo[j]=False; l.remove(j); liberi.append(j)
            elif nuovi:
                w=nuovi[-1]; l.remove(j); occ[w].append(j); pos[j]=w        # lo spazio creato separa i due
            toccati+=nuovi+[v]
        # 3) assorbimento
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
        # 4) decadimento. TC<=0: regola della singolarita' -> solo quando lo spazio e' tornato al grafo minimo, e allora ogni oggetto
        #    composto decade completamente in particelle singole (stessa energia per unita' di massa)
        if TC<=0 and len(R.nbr)<=4:
            for i in ogg:
                if vivo[i] and m[i]>=2:
                    eu=e[i]/m[i]; n=m[i]-1; m[i]=1; e[i]=eu; nfr+=1
                    for _ in range(n):
                        w=random.choice(tuple(R.nbr))
                        if liberi: j=liberi.pop(); m[j]=1; e[j]=eu; vivo[j]=True; pos[j]=w
                        else: j=len(vivo); m.append(1); e.append(eu); vivo.append(True); pos.append(w)
                        occ.setdefault(w,[]).append(j)
        # MASSA CRITICA: un oggetto che supera la frazione MCRIT della massa totale (ed e' caldo almeno TCRIT) decade del tutto
        if MCRIT>0:
            for i in ogg:
                if vivo[i] and m[i]>=MCRIT*NP and e[i]>=TCRIT*m[i]:
                    eu=e[i]/m[i]; n=m[i]-1; m[i]=1; e[i]=eu; nfr+=1; v=pos[i]; qui=[v]+list(R.nbr[v]); critici.append((t,len(R.nbr),n+1,eu))
                    for _ in range(n):
                        w=random.choice(qui)
                        if liberi: j=liberi.pop(); m[j]=1; e[j]=eu; vivo[j]=True; pos[j]=w
                        else: j=len(vivo); m.append(1); e.append(eu); vivo.append(True); pos.append(w)
                        occ.setdefault(w,[]).append(j)
        # altrimenti: decadimento degli oggetti composti tornati caldi
        for i in ogg:
            if TC>0 and vivo[i] and m[i]>=2 and e[i]>=TC*m[i] and random.random()<PF:
                m2=m[i]//2; e2=e[i]*m2/m[i]; m[i]-=m2; e[i]-=e2; v=pos[i]; w=random.choice(tuple(R.nbr[v])); nfr+=1
                if liberi: j=liberi.pop(); m[j]=m2; e[j]=e2; vivo[j]=True; pos[j]=w
                else: j=len(vivo); m.append(m2); e.append(e2); vivo.append(True); pos.append(w)
                occ.setdefault(w,[]).append(j)
        rilassa(toccati)
    ogg=[i for i in range(len(vivo)) if vivo[i]]; V=len(R.nbr)
    res=dict(decadimenti=nfr,storia=storia,vmax=vmax,V=V,oggetti=len(ogg),massa=sum(m[i] for i in ogg),energia=sum(e[i] for i in ogg)+(V-4)*EPS+resto,E0=Etot0,passi=t,eventi=ev)
    # cicli: massimi e minimi dello spazio (serie lisciata)
    Vs=[q[1] for q in storia]; w=5; sm=[sum(Vs[max(0,a-w):a+w+1])/len(Vs[max(0,a-w):a+w+1]) for a in range(len(Vs))]
    amp=max(sm)-min(sm); estremi=[]; verso=1; ref=sm[0]; iref=0
    for a,x in enumerate(sm):
        if verso==1:
            if x>ref: ref=x; iref=a
            elif x<ref-0.25*amp: estremi.append(('max',storia[iref][0],int(ref),storia[iref][2])); verso=-1; ref=x; iref=a
        else:
            if x<ref: ref=x; iref=a
            elif x>ref+0.25*amp: estremi.append(('min',storia[iref][0],int(ref),storia[iref][2])); verso=1; ref=x; iref=a
    res['estremi']=estremi
    if verbose:
        print('massa finale %d, energia finale %.1f (iniziale %.1f), decadimenti %d, eventi %d'%(res['massa'],res['energia'],Etot0,nfr,ev))
        if critici:
            print('DECADIMENTI PER MASSA CRITICA (passo, vertici in quel momento, massa, energia per unita\' di massa):')
            for c in critici: print('   passo %7d: %6d vertici, massa %3d, energia/massa %.2f'%c)
            print('   spazio minimo dopo il primo massimo: %d vertici (4 = singolarita\')'%min(q[1] for q in storia[5:]))
        print('CICLI (spazio lisciato): '+('nessuna oscillazione netta' if len(estremi)<3 else ''))
        for tipo,tt,vv,oo in estremi: print('   %s al passo %7d: %6d vertici, %4d oggetti'%('MASSIMO' if tipo=='max' else 'minimo ',tt,vv,oo))
    return res
if __name__=='__main__':
    a=sys.argv[1:]; f=lambda k,d,ty: ty(a[k]) if len(a)>k else d
    esegui(f(0,300,int),f(1,0.02,float),f(2,0.9,float),f(3,0.05,float),f(4,0.005,float),f(5,400000,int),f(6,0.9,float),f(7,0.01,float),f(8,1,int),f(9,0.0,float),f(10,0.0,float))
