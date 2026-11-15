"""Esperimento 8: l'accorciamento di un oggetto in moto, quando la sua forma e' tenuta da segnali.
Stessa rete casuale dell'esperimento 7: un lampo si allarga di un legame per passo in tutte le direzioni,
senza usare coordinate.
Oggetto: tre marcatori A, B, C.
- C sta di traverso rispetto al moto, a distanza fissa Lp da A: e' l'unita' di misura dell'oggetto.
  (Una scala di riferimento serve comunque: negli oggetti reali la fissano gli atomi. In relativita' le
  lunghezze di traverso al moto non cambiano.)
- B sta lungo il moto e la sua distanza NON e' imposta: B si mette dove l'andata e ritorno di un lampo A->B->A
  dura quanto quella A->C->A. E' la "regola del radar": l'oggetto resta rigido se i suoi tempi interni
  restano uguali in tutte le direzioni (un interferometro di Michelson che si regola da solo).
Lo sperimentatore sposta l'oggetto a velocita' v; le coordinate servono solo a spostarlo e a misurarlo.
Domanda: la distanza di equilibrio A-B si accorcia col moto? La relativita' prevede il fattore sqrt(1 - v^2/c^2).
Uso: python esperimento8.py [N punti] [tic per misura] [seme] [Lp] [velocita' v/c separate da virgole]
Fino a 400.000 punti la rete e' costruita come nell'esperimento 7; oltre, con un metodo che usa meno memoria
(copie dei punti solo in una fascia attorno ai bordi), cosi' si possono usare oggetti piu' grandi rispetto alla grana."""
import sys, math, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree, Delaunay
from esperimento7 import rete_casuale


def rete_casuale_margine(N, rng):
    """rete casuale su un toro: triangolazione di Delaunay con copie dei punti solo in una fascia di bordo"""
    P = rng.random((N, 2))
    m = 12 / math.sqrt(N)
    pts, idx = [P], [np.arange(N)]
    for ox in (-1, 0, 1):
        for oy in (-1, 0, 1):
            if ox == 0 and oy == 0:
                continue
            Q = P + np.array([ox, oy])
            sel = np.nonzero((Q[:, 0] > -m) & (Q[:, 0] < 1 + m) & (Q[:, 1] > -m) & (Q[:, 1] < 1 + m))[0]
            pts.append(Q[sel])
            idx.append(sel)
    X = np.concatenate(pts)
    I = np.concatenate(idx)
    T = Delaunay(X).simplices
    e = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    e = e[(e[:, 0] < N) | (e[:, 1] < N)]
    e = np.sort(I[e], axis=1)
    e = np.unique(e, axis=0)
    return P, e[e[:, 0] != e[:, 1]]


class Banco:
    """la rete casuale, il cercatore del punto piu' vicino e le misure di andata e ritorno"""

    def __init__(self, N, seme):
        self.rng = np.random.default_rng(seme)
        if N <= 400000:
            P, nbr = rete_casuale(N, self.rng)
            lun = np.array([len(x) for x in nbr])
            righe = np.repeat(np.arange(N), lun)
            colonne = np.fromiter((w for x in nbr for w in x), dtype=np.int64, count=int(lun.sum()))
        else:
            P, e = rete_casuale_margine(N, self.rng)
            righe = np.concatenate([e[:, 0], e[:, 1]])
            colonne = np.concatenate([e[:, 1], e[:, 0]])
        self.G = csr_matrix((np.ones(len(colonne), np.float32), (righe, colonne)), shape=(N, N))
        self.tree = cKDTree(P, boxsize=1.0)
        k = np.arange(0, len(righe), max(1, len(righe) // 20000))
        self.ell = float(np.mean(np.linalg.norm(((P[colonne[k]] - P[righe[k]] + 0.5) % 1) - 0.5, axis=1)))

    def andata_ritorno(self, direzione, L, v, ntic, cap):
        """tempo medio (in passi) di un lampo A->X->A con X a distanza L da A, lungo il moto ('lungo') o di traverso;
        l'oggetto si muove a velocita' v (coordinate per passo) in una direzione a caso. Ritorna media ed errore."""
        T = np.empty(ntic)
        passi = np.arange(1, cap + 1)
        for k in range(ntic):
            th = self.rng.random() * 2 * math.pi
            d = np.array([math.cos(th), math.sin(th)])
            perp = np.array([-d[1], d[0]])
            off = (d if direzione == 'lungo' else perp) * L
            xA = self.rng.random(2)
            t = 0
            for tratto in (0, 1):                 # tratto 0: da A a X ; tratto 1: da X ad A
                sorgente = xA + d * v * t + (off if tratto == 1 else 0.0)
                s0 = int(self.tree.query(sorgente % 1.0)[1])
                dist = dijkstra(self.G, directed=True, indices=s0, unweighted=True, limit=cap)
                arrivo = (xA[None, :] + d[None, :] * v * (t + passi[:, None]) + (off[None, :] if tratto == 0 else 0.0)) % 1.0
                nodi = self.tree.query(arrivo)[1]
                ok = np.nonzero(dist[nodi] <= passi)[0]
                if len(ok) == 0:
                    return float('nan'), float('nan')
                t += int(passi[ok[0]])
            T[k] = t
        return float(T.mean()), float(T.std(ddof=1) / math.sqrt(ntic))


def equilibrio(B, b, c, Lg, TAC, eAC, ntic):
    """distanza A-B alla quale l'andata e ritorno lungo il moto dura quanto quella di traverso (TAC).
    Metodo delle secanti partendo da Lg (stima iniziale). Ritorna distanza, errore, misure fatte."""
    v = b * c
    misure = []

    def misura(L):
        cap = int(1.5 * L / (c * (1 - b))) + 40
        T, e = B.andata_ritorno('lungo', L, v, ntic, cap)
        misure.append((L, T, e))
        return T, e

    L1, L2 = 0.95 * Lg, 1.05 * Lg
    T1, e1 = misura(L1)
    T2, e2 = misura(L2)
    for _ in range(2):
        L = L1 + (TAC - T1) * (L2 - L1) / (T2 - T1)
        if min(L1, L2) - 0.01 * Lg <= L <= max(L1, L2) + 0.01 * Lg:
            break
        T3, e3 = misura(L)                       # stima fuori dall'intervallo: nuova misura, tengo le due piu' vicine
        if abs(L - L1) > abs(L - L2):
            L1, T1, e1 = L, T3, e3
        else:
            L2, T2, e2 = L, T3, e3
    k = (T2 - T1) / (L2 - L1)
    L = L1 + (TAC - T1) / k
    w = (L - L1) / (L2 - L1)
    eT = math.sqrt(eAC ** 2 + ((1 - w) * e1) ** 2 + (w * e2) ** 2)
    return L, eT / abs(k), misure


if __name__ == '__main__':
    a = sys.argv[1:]
    N = int(a[0]) if a else 300000
    ntic = int(a[1]) if len(a) > 1 else 150
    seme = int(a[2]) if len(a) > 2 else 1
    Lp = float(a[3]) if len(a) > 3 else 0.12
    velocita = [float(x) for x in a[4].split(',')] if len(a) > 4 else [0.2, 0.4, 0.6, 0.8]
    t0 = time.time()
    B = Banco(N, seme)
    print('rete casuale di %d punti pronta in %.0f s' % (N, time.time() - t0), flush=True)

    # a riposo: unita' di misura (C di traverso) e velocita' massima c
    TAC0, eAC0 = B.andata_ritorno('traverso', Lp, 0.0, 2 * ntic, int(1.6 * Lp / (1.12 * B.ell)) + 40)
    c = 2 * Lp / TAC0
    L0, eL0, _ = equilibrio(B, 0.0, c, Lp, TAC0, eAC0, ntic)
    print('A riposo: C di traverso a %.1f legami, andata e ritorno %.2f passi, c = %.5f per passo' % (TAC0 / 2, TAC0, c))
    print('A riposo: B si mette a %.4f contro %.4f di C (rapporto %.4f +- %.4f: la rete non ha direzioni preferite)'
          % (L0, Lp, L0 / Lp, eL0 / Lp), flush=True)
    print(' v/c | tempo interno (andata e ritorno A-C), rallentamento | distanza di equilibrio A-B rispetto a riposo | relativita\' | differenza')
    for b in velocita:
        g = 1 / math.sqrt(1 - b * b)
        if (Lp / g) * 1.1 / (1 - b) >= 0.45:
            print(' %.2f | oltre il limite della scatola (il lampo farebbe il giro dell\'universo chiuso)' % b, flush=True)
            continue
        TAC, eAC = B.andata_ritorno('traverso', Lp, b * c, ntic, int(1.6 * Lp * g / c) + 40)
        L, eL, mis = equilibrio(B, b, c, Lp * TAC0 / TAC, TAC, eAC, ntic)
        r = L / L0
        er = r * math.sqrt((eL / L) ** 2 + (eL0 / L0) ** 2)
        print(' %.2f |                 %.3f (relativita\' %.3f)          |              %.4f +- %.4f              |    %.4f    |  %+.1f%%'
              % (b, TAC / TAC0, g, r, er, 1 / g, 100 * (r * g - 1)), flush=True)
    print('tempo %.0f s' % (time.time() - t0))
