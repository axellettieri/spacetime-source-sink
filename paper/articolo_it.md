# Spazio creato dagli eventi al posto dell'energia oscura: un modello sorgente–pozzo e il confronto con i dati

**Pasquale Axel Lettieri** ([ORCID 0009-0001-0756-6467](https://orcid.org/0009-0001-0756-6467)) · Ricercatore indipendente, Novedrate (CO), Italia · 15 novembre 2026 · bozza

## Riassunto

Un modello che spiega l'accelerazione dell'espansione senza un'energia oscura come sostanza a sé: al suo posto c'è spazio creato dagli eventi della materia e riassorbito dalle masse. Il modello si accorda con i dati cosmologici meglio del modello standard ΛCDM, con una preferenza nominale fra 2,9 e 3,6 sigma secondo il campione di supernovae. Il modello usa un parametro in più di ΛCDM e uno in meno della parametrizzazione w0wa, che eguaglia.

Il ritmo degli eventi è stato preso in due modi indipendenti: dalla storia osservata della formazione stellare e dal collasso gravitazionale della materia, calcolato con simulazioni a N corpi e con la teoria di Press e Schechter. I due orologi danno lo stesso accordo con i dati di Planck 2018, DESI DR2 e tre campioni di supernovae (Pantheon+, Union3, DES-Dovekie).

Se le masse continuano ad assorbire spazio anche dopo aver esaurito quello creato dagli eventi recenti, il modello prevede che l'espansione si fermi fra circa 30 miliardi di anni e che l'universo ricollassi fra circa 45–50 miliardi di anni. Il modello non spiega la tensione di Hubble, e la sua versione locale è esclusa dai dati.

Una serie di esperimenti numerici esplorativi su grafi mostra che una regola locale di creazione può dare un effetto diffuso, che con un bilancio chiuso dell'energia si ottengono cicli di espansione e collasso, e che un oggetto in moto fatto di segnali interni rallenta il proprio tempo e si accorcia, come nella relatività, con scarti che calano al crescere dell'oggetto rispetto alla grana del grafo. Questi esperimenti sono giocattoli a due dimensioni e non costituiscono una derivazione del modello.

L'analisi non è stata verificata in modo indipendente.

## 1. Introduzione

L'ipotesi di partenza è che lo spaziotempo non sia un contenitore dato ma il prodotto degli eventi fisici: ogni evento che dissipa energia crea spazio, e le masse lo riassorbono. In questa lettura l'energia oscura non è una sostanza: è il conto dello spazio creato e non ancora riassorbito.

L'ipotesi nasce in forma qualitativa dall'immagine dell'universo come grafo causale, in cui i nodi sono eventi e gli archi relazioni di causa ed effetto (Bombelli, Lee, Meyer e Sorkin, 1987). Questo lavoro la traduce in un modello cosmologico minimo, confrontabile con i dati, all'interno della relatività generale.

Il contesto osservativo è favorevole a modelli di questo tipo. Nel 2025 la collaborazione DESI ha riportato una preferenza per un'energia oscura che varia nel tempo, fra 2,8 e 4,2 sigma secondo il campione di supernovae usato, descritta con la parametrizzazione w0wa. La forma preferita dai dati è un'energia oscura che prima cresce e poi cala. Un modello con una sorgente che si esaurisce e un assorbimento che continua produce proprio questa forma.

Gli obiettivi iniziali erano tre: spiegare l'accelerazione dell'espansione, spiegare la tensione di Hubble con un'espansione localmente non uniforme, e ottenere un universo che ricollassa. Il primo e il terzo trovano riscontro, con le condizioni discusse nel testo; il secondo no.

## 2. Il modello

Il modello è un vuoto che scambia energia con la materia, in un universo piatto descritto dalla relatività generale. Lo spazio creato ha densità di energia ρx e si comporta come energia del vuoto; la materia ha densità ρm. Lo scambio Q conserva l'energia totale:

```latex
\dot\rho_x = Q, \qquad \dot\rho_m + 3H\rho_m = -Q, \qquad H^2 = \frac{8\pi G}{3}\,(\rho_m+\rho_r+\rho_x)
```

Lo scambio è la differenza fra una sorgente, proporzionale al ritmo degli eventi, e un pozzo, proporzionale alla massa che assorbe:

```latex
Q = \alpha\,\psi(z)\,(1+z)^3 \;-\; \gamma\,M(z), \qquad M(z)=\int \psi\,dt
```

Qui ψ è il ritmo degli eventi per volume comovente, α lo spazio creato per evento e γ lo spazio assorbito per unità di massa e di tempo. All'inizio ρx è zero: tutta l'energia oscura è creata dagli eventi. I parametri liberi dell'energia oscura sono due (α e γ), contro uno di ΛCDM e tre di w0wa.

**Due orologi degli eventi.** Il ritmo ψ è stato preso in due modi indipendenti.

| Orologio | Che cosa conta come evento | Da dove viene |
|---|---|---|
| Stellare | La formazione delle stelle | Storia osservata (Madau e Dickinson 2014) |
| Emergente | La materia che entra per la prima volta in una regione collassata | Simulazioni a N corpi di questo lavoro e teoria di Press e Schechter |

Con l'orologio emergente la sorgente è proporzionale alla materia che collassa nell'unità di tempo e il pozzo alla frazione di materia già collassata.

**Varianti provate.** Una sorgente non diluita dall'espansione è esclusa dai dati. Per il pozzo sono state provate tre forme (proporzionale allo spazio creato, costante, proporzionale alla massa accumulata): danno accordi equivalenti con i dati di oggi ma destini diversi. La forma proporzionale allo spazio creato dà espansione eterna; le altre due portano al collasso.

**Perturbazioni.** Lo scambio è trattato come diffuso: la materia perde o riceve energia in modo uniforme, senza scambio di quantità di moto nel proprio sistema di riferimento. La crescita delle strutture segue l'equazione corrispondente; nello spettro del fondo cosmico il modello entra come fluido efficace, con una correzione all'ampiezza della lente gravitazionale calcolata dalle equazioni di crescita (circa 0,98).

## 3. Dati e metodo

Il modello è confrontato con tre famiglie di dati: supernovae, oscillazioni acustiche dei barioni e fondo cosmico a microonde.

| Dati | Dettaglio |
|---|---|
| Supernovae Pantheon+ | 1590 curve di luce di supernovae con la selezione ufficiale (redshift sopra 0,01), matrice di covarianza completa |
| Supernovae Union3 e DES-Dovekie | Campioni alternativi, usati come controllo |
| DESI DR2 | 13 misure di oscillazioni acustiche dei barioni |
| Planck 2018 | Spettro completo TT, TE, EE (verosimiglianza plik-lite) con vincolo sulla profondità ottica; in alternativa i vincoli compressi sulle distanze per i confronti rapidi |
| Crescita delle strutture | 7 misure di fσ8 (6dFGS, SDSS MGS, BOSS DR12, eBOSS DR16), dalla compilazione di Alam e altri 2021, usate come controllo |

**Due livelli di confronto.** Il confronto rapido usa i vincoli compressi di Planck e richiede secondi; serve a esplorare le varianti. Il confronto completo calcola lo spettro del fondo cosmico con CAMB e usa l'orizzonte sonoro calcolato da CAMB; è quello a cui si riferiscono i risultati principali.

**Statistica.** Per ogni modello si cerca il minimo del χ² e lo si confronta con quello di ΛCDM sugli stessi dati. La differenza è tradotta in sigma nominali per il numero di parametri in più. I modelli non sono annidati in modo esatto, quindi i sigma sono indicativi. Per il modello con orologio stellare è stata fatta anche una catena di Markov completa (32 catene per 3000 passi) per ricavare le incertezze sul futuro.

**Simulazioni.** L'orologio emergente viene da simulazioni a N corpi con metodo a griglia: 512³ particelle, scatole da 125 a 1000 Mpc, su una scheda grafica. Senza eventi, la crescita delle strutture a grande scala riproduce la teoria lineare entro il 2%.

## 4. Risultati

### 4.1 Confronto principale

Con Planck completo, DESI DR2 e Pantheon+ il modello sorgente–pozzo migliora il χ² di ΛCDM di 8,3–8,8 punti, quanto w0wa, con un parametro in meno di quest'ultimo.

| Modello | Parametri in più di ΛCDM | χ² | Differenza da ΛCDM |
|---|---|---|---|
| ΛCDM | 0 | 2008,70 | 0 |
| w0wa | 2 | 2000,01 | −8,7 |
| Sorgente–pozzo, orologio stellare, pozzo sulla massa accumulata | 1 | 2000,36 | −8,3 |
| Sorgente–pozzo, orologio stellare, pozzo sullo spazio creato | 1 | 1999,94 | −8,8 |
| Sorgente–pozzo, orologio emergente (simulazione) | 1, più la soglia dell'orologio | 2000,03 | −8,7 |
| Sorgente–pozzo, orologio emergente (teoria) | 1, più la massa minima | 2000,31 | −8,4 |

La storia dell'espansione del modello coincide con quella del miglior w0wa entro lo 0,2% fino a redshift 5. La costante di Hubble risultante è 67,7 km/s/Mpc.

### 4.2 Servono entrambi gli ingredienti

Nel confronto rapido la sola sorgente peggiora l'accordo (+9,1) e il solo pozzo lo migliora poco (da −1,4 a −2,1); insieme danno da −7,3 a −7,7. La forma esatta dell'orologio non è determinante: una sorgente generica che cala con l'espansione dà −6,7.

### 4.3 L'orologio emergente

L'orologio calcolato dal collasso della materia vale quanto quello osservato, a condizione che gli eventi comincino presto: fra il 7 e l'11% della materia collassata a redshift 3 e fra il 55 e il 61% oggi. Con la teoria di Press e Schechter questo corrisponde ad aloni sopra 10^10–10^11 masse solari, cioè quelli che ospitano le galassie; con masse minime di 10^12 o più l'accordo si perde. Nelle simulazioni il risultato si ritrova a tutte le risoluzioni provate e con semi diversi, ma la soglia di densità che lo produce cambia con la risoluzione: il metodo a griglia non basta a fissare la massa minima.

### 4.4 Altri campioni di supernovae

| Supernovae | Orologio stellare | Orologio emergente | Sigma nominali |
|---|---|---|---|
| Pantheon+ | −8,3 | −8,7 | 2,9 |
| Union3 | −13,0 | −12,7 | 3,6 |
| DES-Dovekie | −10,9 | −11,0 | 3,3 |

Contando come parametro anche la soglia dell'orologio emergente, i sigma scendono a 2,5, 3,1 e 2,9.

### 4.5 Risultati negativi

- **Tensione di Hubble.** Non è spiegata: il modello dà H0 = 67,7, come ΛCDM.
- **Espansione locale più veloce.** Nei dati delle velocità delle galassie vicine non c'è traccia di un'espansione locale maggiore, entro l'1%; l'8% necessario per la tensione è escluso a oltre 15 sigma. Nessuna dipendenza dall'ambiente (−0,6 ± 0,6%), misurato con la mappa di densità e velocità 2M++ (Carrick e altri 2015).
- **Versione locale in simulazione.** Se lo spazio creato resta vicino agli eventi, l'espansione locale aumenta, ma al massimo di 3,6 punti percentuali entro 75 Mpc e 0,2 entro 300 Mpc: almeno dieci volte meno del necessario.
- **Scambio locale.** L'assorbimento locale di spazio da parte delle masse è escluso; la restituzione locale di energia a tutta la materia è esclusa dai limiti sulla variazione di massa del Sole; quella alla sola materia oscura è sfavorita a circa 3 sigma.
- **Media degli eventi invece dell'accumulo.** Un'energia oscura proporzionale alla densità media di eventi nel passato visibile peggiora l'accordo di 104 punti.

### 4.6 Bilancio dell'energia

La sorgente converte circa il 9% della materia in spazio, pari a circa 4,5 kg per ogni kg di stelle formate. È circa 3000 volte l'energia che le stelle dissipano: l'evento agisce da innesco, non da fonte dell'energia, e il modello non spiega questo fattore.

### 4.7 Materia, lenti gravitazionali e crescita delle strutture

Il modello prevede che oggi la materia sia il 39–40% del totale, contro il 30% di ΛCDM, ma questo non entra in conflitto con le misure indipendenti di lente gravitazionale e di crescita. La materia restituita dal pozzo è distribuita in modo uniforme: aumenta la densità media e diluisce i contrasti nella stessa proporzione, così il prodotto dei due, che è ciò che la lente misura, resta quello di ΛCDM.

| Redshift | Materia rispetto a ΛCDM | Contrasto rispetto a ΛCDM | Sorgente della lente rispetto a ΛCDM |
|---|---|---|---|
| 0 | 1,31 | 0,78 | 1,03 |
| 0,2 | 1,10 | 0,92 | 1,01 |
| 0,5 | 0,99 | 1,01 | 1,00 |
| 1 | 0,97 | 1,04 | 1,01 |
| 2 | 0,99 | 1,03 | 1,01 |

Valori per l'orologio stellare; con l'orologio emergente differiscono di meno dello 0,5% nell'ultima colonna. Pesata sulle distanze tipiche della lente debole, l'ampiezza del segnale è 1,004 volte quella di ΛCDM (w0wa: 1,017).

| Misura | Valore osservato | Posizione del modello |
|---|---|---|
| Lente debole, KiDS-Legacy ([arXiv:2503.19441](https://arxiv.org/abs/2503.19441)) | S8 = 0,815 (+0,016, −0,021) | Come ΛCDM con i parametri di Planck (circa 0,81–0,82): in accordo |
| Lente debole e raggruppamento, DES anno 6 ([arXiv:2601.14559](https://arxiv.org/abs/2601.14559)) | S8 = 0,789 ± 0,012 | Come ΛCDM: circa 2 sigma sopra, la stessa lieve tensione del modello standard |
| Velocità di crescita, 7 misure di fσ8 |  | χ² = 5,8 contro 6,3 di ΛCDM |

Due conseguenze. La prima è che il controllo è superato: l'aumento di materia non è un problema. La seconda è che il modello risulta quasi indistinguibile da w0wa anche in queste misure: lente, raggruppamento delle galassie e velocità di crescita vedono tutte la stessa combinazione di materia e contrasto. Il conteggio degli ammassi di galassie non è stato calcolato; ci si aspetta che cambi poco, perché la materia aggiunta è diffusa e recente, ma resta da verificare.

### 4.8 Una previsione verificabile: la linea nel piano (w0, wa)

In parole semplici: chi misura l'energia oscura pubblica i risultati con due numeri, w0 (com'è oggi) e wa (quanto era diversa nel passato). Qui il modello non viene sostituito da w0wa, ma tradotto in quei due numeri per confrontarlo con le misure pubblicate. Tradotto, dice che l'energia oscura è cresciuta nel passato ed è in calo oggi.

Il modello ha un solo parametro libero in più di ΛCDM, l'intensità del pozzo. Variandolo, il punto (w0, wa) si muove lungo una linea, come un treno sul binario, mentre w0wa può occupare tutto il piano. Se le misure future cadranno fuori da quella linea, il modello sarà smentito.

Il modello ha un parametro in meno di w0wa, quindi non può occupare tutto il piano (w0, wa): i suoi punti stanno su una linea. Per ogni intensità del pozzo si è cercata la coppia (w0, wa) che riproduce la storia di espansione del modello (entro lo 0,5% fino a redshift 2,5). Con l'orologio stellare e il pozzo sulla massa accumulata la linea è:

```latex
w_a \simeq -0{,}41 \;-\; 1{,}61\,(1+w_0)
```

Con le altre due forme del pozzo la pendenza è −1,46 (pozzo costante) e −1,23 (pozzo sullo spazio creato); il punto di partenza è lo stesso. La linea non passa per ΛCDM (w0 = −1, wa = 0): senza pozzo il modello dà w0 = −0,97 e wa = −0,46, ed è escluso dai dati.

| Supernovae | Modello: w0, wa | w0wa libero sugli stessi dati: w0, wa |
|---|---|---|
| Pantheon+ | −0,82, −0,70 | −0,86, −0,51 |
| DES-Dovekie | −0,81, −0,72 | −0,83, −0,63 |
| DES anno 5 (versione precedente) | −0,76, −0,79 | −0,77, −0,77 |
| Union3 | −0,72, −0,85 | −0,70, −0,98 |

Confronto rapido (Planck compresso). In tutti e quattro i campioni il punto preferito da w0wa libero cade entro 0,2 in wa dalla linea del modello, cioè dentro l'incertezza attuale su wa (circa ±0,2). Con Pantheon+ il modello fissa wa fra −0,75 e −0,65, un intervallo circa quattro volte più stretto di quello di w0wa libero.

Questa è la previsione che i prossimi dati possono smentire: quando l'incertezza su wa scenderà sotto 0,1, il punto misurato dovrà trovarsi sulla linea. La direzione lungo cui i dati attuali sono incerti ha pendenza diversa (circa −2,9), quindi le due linee si separano allontanandosi dal punto di incrocio, vicino a w0 = −0,78.

## 5. Il futuro dell'espansione

Con il pozzo proporzionale alla massa accumulata il modello prevede che l'espansione si fermi fra circa 30 miliardi di anni e che l'universo ricollassi. Due orologi indipendenti danno le stesse date.

|  | Orologio stellare (catena completa) | Orologio emergente (simulazione) |
|---|---|---|
| Spazio creato esaurito | fra 17 miliardi di anni (13–26) | fra 16 |
| Espansione ferma | fra 31 (26–41) | fra 29 |
| Dimensione massima | 3,2 volte l'attuale (2,6–4,7) | 3,0 |
| Collasso | fra 49 (42–62) | fra 45, fine della simulazione a un decimo della dimensione attuale |

Fra parentesi l'intervallo al 68%. Questi valori sostituiscono quelli più approssimati del confronto rapido.

**La previsione è condizionata.** L'espansione si ferma perché l'energia dello spazio creato scende sotto zero: le masse continuano ad assorbire dopo aver esaurito lo spazio creato dagli eventi contati dal modello. Se l'assorbimento si fermasse a zero, resterebbe un universo di sola materia che si espande per sempre. I dati di oggi non distinguono i due casi. Nell'ipotesi di partenza tutto lo spazio, anche quello dell'espansione primordiale, nasce da eventi, e l'assorbimento oltre lo zero è naturale; ma il modello non conta gli eventi primordiali, e questa parte resta un'assunzione.

**Il collasso simulato.** Nella simulazione a N corpi portata oltre l'inversione, alla fine il 99,9% della materia è collassato. La materia cresce di circa 80 volte durante il ciclo, perché il pozzo continua a restituire energia alla materia mentre il debito di spazio aumenta. Non compare alcun rimbalzo: il modello descrive il collasso, non la ripartenza.

## 6. Esperimenti esplorativi su grafi

Le equazioni del modello non sono derivate dal grafo causale. Otto esperimenti numerici esplorano se le idee qualitative dell'ipotesi reggono su grafi senza coordinate. Sono giocattoli, quasi tutti a due dimensioni, con regole scelte a mano: mostrano che un comportamento è possibile, non che sia quello vero.

| Esperimento | Domanda | Risultato |
|---|---|---|
| 1. Crescita di un insieme causale con sorgente e pozzo | Che cosa separa gli eventi? | La larghezza del grafo segue la probabilità di legame; il grafo collassa a una catena solo se c'è il pozzo |
| 2. Gas su un anello | La regola "gli eventi creano spazio, le masse lo assorbono" produce un'espansione? | Espansione, arresto e collasso senza equazioni imposte; espansione non uniforme, che si media su zone grandi |
| 3. Rete di triangoli senza coordinate | Lo spazio creato in un punto resta lì? | Senza rilassamento i vertici nuovi non creano distanza; con poco rilassamento lo spazio resta locale; con rilassamento forte si sparge ovunque (93%) |
| 4. Materia ed eventi sulla rete, energia conservata | Le due singolarità sono uguali? | Da 4 vertici a circa 11.000 e ritorno a 4; massa ed energia identiche all'inizio e alla fine, ma 300 particelle all'inizio e un solo oggetto alla fine |
| 5. Decadimento | Che cosa fa ripartire un ciclo? | Singolarità instabile: cicli con lo stesso massimo entro il 2%. Massa critica: cicli senza singolarità, massimi variabili del 15%. Decadimento di ogni oggetto caldo: equilibrio irregolare, nessun collasso |
| 6. Segnali sulla rete | La velocità massima è uguale in ogni direzione? Che cosa fa la creazione di spazio? | Rete regolare: 10% di differenza fra direzioni; rete casuale: 2,4%, compatibile con il rumore. La creazione di spazio ritarda i segnali e fa comparire un orizzonte |
| 7. Orologio in moto | Un oggetto in moto compie meno eventi interni? | Sì: il rallentamento si avvicina al fattore della relatività al crescere della dimensione dell'oggetto (2,03 contro 2,29 a nove decimi della velocità massima) |
| 8. Oggetto tenuto insieme da segnali | Un oggetto in moto si accorcia? | Sì, se la sua lunghezza lungo il moto è fissata da segnali invece che imposta. A otto decimi della velocità massima il fattore è 0,631 per un oggetto di 97 legami, contro 0,600 della relatività; lo scarto cala con la dimensione (12,5%, 9,0%, 5,1% per 25, 49, 97 legami) |

Tre indicazioni utili per l'ipotesi.

- **Località della regola, diffusione dell'effetto.** I dati escludono un'energia oscura concentrata vicino agli eventi. Sul grafo la creazione può essere locale e l'effetto diffuso, se la geometria si rilassa molto più in fretta di quanto lo spazio venga creato. L'esperimento non dice perché dovrebbe essere così.
- **Cicli.** Con un bilancio chiuso dell'energia l'universo torna allo stato minimo con la stessa massa e la stessa energia. Per ripartire serve un'instabilità che scatti solo nello stato estremo, o al superamento di una massa critica.
- **Tempo e moto.** Anche su un grafo con un passo unico per tutti, il tempo di un oggetto fatto di segnali interni rallenta quando l'oggetto si muove. Anche la contrazione delle lunghezze compare, se la lunghezza dell'oggetto lungo il moto è fissata da segnali invece che imposta (esperimento 8). Entrambi gli effetti si avvicinano a quelli della relatività quanto più l'oggetto è grande rispetto alla grana del grafo. Il grafo però ha un sistema di riferimento privilegiato, il proprio: un teorema di Bombelli, Henson e Sorkin mostra che una rete con un numero finito di vicini per nodo non può rappresentare uno spaziotempo discreto senza sceglierne uno, mentre gli insiemi causali costruiti a caso nello spaziotempo lo evitano.

## 7. Discussione

### 7.1 Punti di forza

- **Pochi parametri con un significato fisico.** Spazio creato per evento e spazio assorbito per massa, contro i due parametri descrittivi di w0wa.
- **Due orologi indipendenti.** La formazione stellare osservata e il collasso gravitazionale calcolato danno lo stesso accordo, su tre campioni di supernovae.
- **Una previsione che distingue.** Un'energia oscura che diventa negativa e un'espansione che si ferma, a differenza di ΛCDM.

### 7.2 Limiti

- **Il segnale non è nuovo.** La preferenza statistica è in buona parte la stessa che DESI riporta per w0wa. Il modello ne offre un'interpretazione fisica; non la scopre e non la rafforza oltre il risparmio di un parametro.
- **La significatività è sotto la soglia di una scoperta**, e dipende dal campione di supernovae.
- **Manca il meccanismo.** Non c'è una teoria che dica perché un evento dovrebbe creare spazio, né che giustifichi un'energia creata circa 3000 volte superiore a quella dissipata.
- **L'orologio emergente ha un parametro scelto sui dati** (la soglia di collasso o la massa minima), anche se il valore cade nell'intervallo fisicamente atteso.
- **Le perturbazioni sono trattate in modo approssimato.** Il fluido efficace e la correzione di lente valgono per lo scambio diffuso; effetto Sachs–Wolfe integrato, ricostruzione della lente e conteggio degli ammassi non sono stati confrontati con i dati in modo completo.
- **Il legame con il grafo causale è qualitativo.** Gli esperimenti della sezione 6 non derivano le equazioni della sezione 2.
- **Poche previsioni che lo distinguano.** Espansione, lente e crescita coincidono con quelle di w0wa entro l'1–2%. Ciò che distingue il modello è la linea obbligata nel piano (w0, wa) della sezione 4.8, il legame con la storia delle strutture e il futuro, che non è osservabile.
- **Nessuna verifica indipendente.** Analisi e codice sono stati prodotti in sessioni di lavoro assistite da un modello di intelligenza artificiale; durante il lavoro sono stati trovati e corretti diversi errori, e altri possono essere rimasti.

### 7.3 Relazione con la letteratura

| Lavoro | Relazione con questo modello |
|---|---|
| DESI DR2 (arXiv:2503.14738) | Fonte della preferenza per un'energia oscura variabile; i nostri valori per w0wa riproducono i loro |
| Vuoto interagente (letteratura ampia) | Il modello ne è un caso particolare, con lo scambio legato a un orologio astrofisico |
| Croker e altri, buchi neri accoppiati alla cosmologia (arXiv:2405.12282) | Stessa idea di una sorgente legata alla formazione stellare, senza pozzo; la nostra variante con sola sorgente riproduce il loro H0 vicino a 70 |
| Li e Zhang, inversione di segno cosmica (arXiv:2506.18477) | Ricostruiscono dai dati DESI DR2, senza assumere una forma, lo scambio di energia fra materia oscura ed energia oscura: lo scambio cambia verso nel tempo, dalla materia all'energia oscura nel passato e al contrario oggi. È il comportamento del nostro modello, in cui prima prevale la sorgente e poi il pozzo |
| Luu, Qiu e Tye, durata dell'universo (arXiv:2506.24011) | Collasso previsto da un modello diverso (un campo leggero più una costante cosmologica negativa), con un'età totale di circa 33 miliardi di anni |
| Sorkin, costante cosmologica sempre presente | Energia oscura dagli insiemi causali, da fluttuazioni casuali invece che da eventi della materia |
| Bombelli, Henson e Sorkin, discretezza senza rottura di simmetria (arXiv:gr-qc/0605006) | Spiega perché le reti degli esperimenti della sezione 6, con pochi vicini per nodo, hanno per forza un sistema di riferimento privilegiato, e perché gli insiemi causali costruiti a caso nello spaziotempo no |

Le proposte di una proporzione variabile fra spazio e tempo nell'universo primordiale (velocità della luce variabile, gravità di Hořava–Lifshitz, triangolazioni dinamiche causali) sono vicine alle idee qualitative discusse nella sezione 6, ma non sono state usate nei calcoli.

## 8. Conclusioni e lavoro futuro

Spazio creato dagli eventi della materia e riassorbito dalle masse, al posto di un'energia oscura come sostanza a sé, descrive i dati cosmologici attuali quanto la migliore parametrizzazione in uso, con un parametro in meno e una lettura fisica. Il risultato regge sostituendo l'orologio osservato con uno calcolato. La preferenza rispetto al modello standard è fra 2,9 e 3,6 sigma nominali: un indizio, non una scoperta.

L'ipotesi di partenza ne esce confermata in una parte e smentita in un'altra. La parte confermata è l'accumulo: lo spazio creato resta, e i dati accettano questa forma mentre rifiutano le alternative provate. La parte smentita è la località dell'effetto: l'espansione non è più veloce dove ci sono più eventi, e la tensione di Hubble resta inspiegata.

Lavoro futuro, in ordine di priorità:

1. Verifica indipendente dei calcoli principali da parte di un cosmologo.
1. Catena completa delle incertezze con l'orologio emergente e lo spettro completo di Planck.
1. Confronto con l'effetto Sachs–Wolfe integrato e con la ricostruzione della lente del fondo cosmico.
1. Orologio emergente con un vero cercatore di aloni, per fissare la massa minima senza guardare i dati.
1. Un insieme causale che cresce evento per evento in più dimensioni, con la materia fatta di catene di eventi: è il passo che legherebbe le equazioni al grafo, ed è un problema aperto anche nella ricerca.

## Glossario

| Termine | Significato |
|---|---|
| Energia oscura | Ciò che fa accelerare l'espansione dell'universo. Non si sa che cosa sia; se ne misura solo l'effetto. |
| ΛCDM | Il modello standard della cosmologia: energia oscura costante nel tempo (Λ) più materia oscura fredda (CDM). |
| w0wa | Non una teoria ma una descrizione: lascia che l'energia oscura cambi nel tempo e usa due numeri per dire come. w0 dice com'è oggi, wa quanto cambia andando nel passato. Con w0 = −1 e wa = 0 si ritrova ΛCDM. |
| Parametro libero | Un numero che la teoria non fissa e che si ricava dai dati. A parità di accordo con i dati, meno parametri è meglio. |
| χ² (chi quadro) | Misura di quanto un modello si discosta dai dati: più è basso, migliore è l'accordo. Qui conta la differenza rispetto a ΛCDM. |
| Sigma | Quanto è improbabile che una differenza sia dovuta al caso. A 3 sigma si parla di indizio, a 5 di scoperta. |
| Redshift (z) | Misura di quanto è lontano nel tempo ciò che si osserva: z = 0 è oggi, z = 1 circa 8 miliardi di anni fa, z = 3 circa 11,5. |
| Sorgente | Nel modello, la creazione di spazio da parte degli eventi. |
| Pozzo | Nel modello, il riassorbimento di spazio da parte delle masse. |
| Spazio creato | Lo spazio prodotto dagli eventi e non ancora riassorbito. Nel modello prende il posto dell'energia oscura: produce lo stesso effetto senza essere una sostanza a sé. |
| Evento | Un fatto fisico che dissipa energia. Nei calcoli: la nascita di stelle, oppure materia che cade in una regione già addensata. |
| Orologio degli eventi | La curva che dice quanti eventi avvengono in ogni epoca. Stellare: dalle osservazioni delle stelle. Emergente: calcolato dall'addensarsi della materia. |
| Sorgente diluita, non diluita | Due modi di contare gli eventi mentre l'universo si espande. Diluita: conta gli eventi per volume reale, che cresce con l'espansione. Non diluita: li conta per volume che segue l'espansione. I dati accettano solo la prima. |
| Vuoto che scambia energia | La forma matematica del modello: lo spazio creato prende energia dalla materia quando nasce e la restituisce quando è riassorbito. |
| Scambio diffuso, locale | Diffuso: l'energia è tolta e restituita alla materia in modo uniforme ovunque. Locale: solo vicino agli eventi o alle masse. I dati escludono la versione locale. |
| Materia collassata, alone | Materia che la gravità ha raccolto in un grumo stabile. Un alone è uno di questi grumi; quelli grandi ospitano le galassie. |
| Soglia di collasso | Nelle simulazioni, la densità oltre la quale una regione è considerata collassata. |
| Simulazione a N corpi | Calcolo che segue milioni di particelle che si attraggono per gravità, per vedere come si formano le strutture. |
| Costante di Hubble (H0) | La velocità di espansione dell'universo oggi. |
| Tensione di Hubble | Il disaccordo fra H0 misurato vicino a noi (circa 73) e quello dedotto dall'universo primordiale (circa 67,5). |
| Planck | Satellite che ha misurato la luce più antica dell'universo (fondo cosmico a microonde). |
| DESI | Strumento che misura la distanza di milioni di galassie usando un'impronta lasciata dal suono nell'universo primordiale (oscillazioni acustiche dei barioni). |
| Supernovae | Esplosioni stellari di luminosità nota, usate per misurare le distanze. Pantheon+, Union3 e DES sono tre raccolte diverse. |
| Lente gravitazionale debole, S8 | La materia deforma leggermente le immagini delle galassie lontane. S8 è il numero che riassume quanta materia c'è e quanto è raggruppata. |
| Contrasto | Quanto una zona è più densa della media. |
| Mpc (megaparsec) | Unità di distanza: circa 3,3 milioni di anni luce. |
| Inversione, collasso | Inversione: il momento in cui l'espansione si ferma. Collasso: la contrazione che segue. |
| Grafo causale, insieme causale | Rappresentazione dell'universo come rete: i nodi sono eventi, i legami dicono quale evento può influenzare quale. |
| Vertice, legame | Sui grafi degli esperimenti: un vertice è un atomo di spazio, un legame collega due vertici vicini. |
| Rilassamento | Negli esperimenti sui grafi, il riassestarsi della rete verso una forma regolare dopo che sono stati aggiunti vertici. |

## 9. Riproducibilità e riferimenti

Tutti i risultati si rifanno con il codice dell'archivio modello_sorgente_pozzo_v18.zip (cartella `code/`); gli esperimenti sui grafi dal 4 al 7 sono stati rieseguiti su una seconda macchina con numeri identici. Il codice è nella cartella `code/` di questo repository; il diario completo del lavoro, con le correzioni fatte in corso d'opera, è disponibile su richiesta. Versione inglese: [article_en.md](article_en.md).

| Parte | Script |
|---|---|
| Modello e confronto rapido | fit3.py |
| Planck completo, orologio stellare | cmb_full.py, mcmc_cmb.py |
| Planck completo, orologio emergente | cmb_orologio.py |
| Simulazioni a N corpi | universo3d.py, adatta3d.py, corsa3d.sh, collasso.sh, scala_fine.sh |
| Esperimenti sui grafi | grafo/percolazione.py, grafo/esperimento1.py … esperimento8.py |
| Catena delle incertezze, orologio emergente | mcmc_orologio.py, analizza_orologio.py, catena_orologio.sh |

**Riferimenti**

- Ahmed, Dodelson, Greene, Sorkin (2004), Everpresent Λ, Physical Review D 69, 103523, [arXiv:astro-ph/0209274](https://arxiv.org/abs/astro-ph/0209274).
- Alam e altri, eBOSS Collaboration (2021), Completed SDSS-IV extended Baryon Oscillation Spectroscopic Survey: cosmological implications from two decades of spectroscopic surveys at the Apache Point Observatory, Physical Review D 103, 083533, [arXiv:2007.08991](https://arxiv.org/abs/2007.08991). Fonte della compilazione delle misure di fσ8.
- Beutler e altri (2012), The 6dF Galaxy Survey: z ≈ 0 measurements of the growth rate and σ8, Monthly Notices of the Royal Astronomical Society 423, 3430, [arXiv:1204.4725](https://arxiv.org/abs/1204.4725).
- Bombelli, Henson, Sorkin (2009), Discreteness without symmetry breaking: a theorem, Modern Physics Letters A 24, 2579, [arXiv:gr-qc/0605006](https://arxiv.org/abs/gr-qc/0605006).
- Bombelli, Lee, Meyer, Sorkin (1987), Space-time as a causal set, Physical Review Letters 59, 521.
- Carrick, Turnbull, Lavaux, Hudson (2015), Cosmological parameters from the comparison of peculiar velocities with predictions from the 2M++ density field, Monthly Notices of the Royal Astronomical Society 450, 317, [arXiv:1504.04627](https://arxiv.org/abs/1504.04627). Mappa usata nella versione distribuita da [pvhub](https://github.com/KSaid-1/pvhub).
- Croker e altri (2024), DESI dark energy time evolution is recovered by cosmologically coupled black holes, JCAP 10 (2024) 094, [arXiv:2405.12282](https://arxiv.org/abs/2405.12282).
- DES Collaboration (2026), Dark Energy Survey Year 6 results: cosmological constraints from galaxy clustering and weak lensing, accettato da Physical Review D, [arXiv:2601.14559](https://arxiv.org/abs/2601.14559).
- DESI Collaboration (2025), DESI DR2 results II: measurements of baryon acoustic oscillations and cosmological constraints, Physical Review D 112, [arXiv:2503.14738](https://arxiv.org/abs/2503.14738).
- Eisenstein, Hu (1998), Baryonic features in the matter transfer function, Astrophysical Journal 496, 605, [arXiv:astro-ph/9709112](https://arxiv.org/abs/astro-ph/9709112).
- Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: the MCMC hammer, PASP 125, 306, [arXiv:1202.3665](https://arxiv.org/abs/1202.3665).
- Lewis, Challinor, Lasenby (2000), Efficient computation of CMB anisotropies in closed FRW models (CAMB), Astrophysical Journal 538, 473, [arXiv:astro-ph/9911177](https://arxiv.org/abs/astro-ph/9911177).
- Li, Zhang (2025), Cosmic sign-reversal: non-parametric reconstruction of interacting dark energy with DESI DR2, [arXiv:2506.18477](https://arxiv.org/abs/2506.18477).
- Luu, Qiu, Tye (2025), The lifespan of our universe, JCAP (2025), [arXiv:2506.24011](https://arxiv.org/abs/2506.24011).
- Madau, Dickinson (2014), Cosmic star-formation history, Annual Review of Astronomy and Astrophysics 52, 415, [arXiv:1403.0007](https://arxiv.org/abs/1403.0007).
- Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Astronomy & Astrophysics 641, A6, [arXiv:1807.06209](https://arxiv.org/abs/1807.06209).
- Popovic e altri (2025), The Dark Energy Survey supernova program: a reanalysis of cosmology results and evidence for evolving dark energy with an updated Type Ia supernova calibration (DES-Dovekie), [arXiv:2511.07517](https://arxiv.org/abs/2511.07517).
- Press, Schechter (1974), Formation of galaxies and clusters of galaxies by self-similar gravitational condensation, Astrophysical Journal 187, 425.
- Prince, Dunkley (2019), Data compression in cosmology: a compressed likelihood for Planck data (planck-lite-py), Physical Review D 100, 083502, [arXiv:1909.05869](https://arxiv.org/abs/1909.05869).
- Rideout, Sorkin (2000), A classical sequential growth dynamics for causal sets, Physical Review D 61, 024002, [arXiv:gr-qc/9904062](https://arxiv.org/abs/gr-qc/9904062).
- Rubin e altri (2025), Union through UNITY: cosmology with 2,000 SNe using a unified Bayesian framework (Union3), Astrophysical Journal, [arXiv:2311.12098](https://arxiv.org/abs/2311.12098).
- Scolnic e altri (2022), The Pantheon+ analysis: the full data set and light-curve release, Astrophysical Journal 938, 113, [arXiv:2112.03863](https://arxiv.org/abs/2112.03863).
- Wright e altri (2025), KiDS-Legacy: cosmological constraints from cosmic shear with the complete Kilo-Degree Survey, Astronomy & Astrophysics 703, A158, [arXiv:2503.19441](https://arxiv.org/abs/2503.19441).

Riferimenti verificati il 15 novembre 2026 sulle pagine di arXiv, NASA ADS o delle riviste. Le 7 misure di fσ8 sono quelle riassunte da Alam e altri (2021): 6dFGS (Beutler e altri 2012), SDSS MGS, BOSS DR12 e i tre campioni eBOSS DR16 (galassie rosse, galassie a righe di emissione, quasar).
