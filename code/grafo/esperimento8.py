"""Experiment 8: the contraction of a moving object, when its shape is held by signals.
Same random network as experiment 7: a flash spreads by one link per step in all directions,
without using coordinates.
Object: three markers A, B, C.
- C sits across the motion, at fixed distance Lp from A: it is the object's unit of length.
  (A reference scale is needed anyway: in real objects atoms set it. In relativity
  lengths across the motion do not change.)
- B sits along the motion and its distance is NOT imposed: B places itself where the round trip of a flash A->B->A
  lasts as long as A->C->A. This is the "radar rule": the object stays rigid if its internal times
  stay equal in all directions (a self-adjusting Michelson interferometer).
The experimenter moves the object at speed v; coordinates are used only to move and measure it.
Question: does the equilibrium distance A-B shrink with motion? Relativity predicts the factor sqrt(1 - v^2/c^2).
Usage: python esperimento8.py [N points] [ticks per measurement] [seed] [Lp] [speeds v/c separated by commas]
Up to 400,000 points the network is built as in experiment 7; beyond, with a method that uses less memory
(copies of the points only in a band around the edges), so larger objects relative to the grain can be used."""
import sys, math, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree, Delaunay
from esperimento7 import rete_casuale


def rete_casuale_margine(N, rng):
    """random network on a torus: Delaunay triangulation with copies of the points only in an edge band"""
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
    """the random network, the nearest-point finder and the round-trip measurements"""

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
        """mean time (in steps) of a flash A->X->A with X at distance L from A, along the motion ('lungo') or across;
        the object moves at speed v (coordinates per step) in a random direction. Returns mean and error."""
        T = np.empty(ntic)
        passi = np.arange(1, cap + 1)
        for k in range(ntic):
            th = self.rng.random() * 2 * math.pi
            d = np.array([math.cos(th), math.sin(th)])
            perp = np.array([-d[1], d[0]])
            off = (d if direzione == 'lungo' else perp) * L
            xA = self.rng.random(2)
            t = 0
            for tratto in (0, 1):                 # leg 0: from A to X ; leg 1: from X to A
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
    """A-B distance at which the round trip along the motion lasts as long as the one across (TAC).
    Secant method starting from Lg (initial guess). Returns distance, error, measurements made."""
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
        T3, e3 = misura(L)                       # estimate outside the interval: new measurement, keep the two closest
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
    print('random network of %d points ready in %.0f s' % (N, time.time() - t0), flush=True)

    # at rest: unit of length (C across) and maximum speed c
    TAC0, eAC0 = B.andata_ritorno('traverso', Lp, 0.0, 2 * ntic, int(1.6 * Lp / (1.12 * B.ell)) + 40)
    c = 2 * Lp / TAC0
    L0, eL0, _ = equilibrio(B, 0.0, c, Lp, TAC0, eAC0, ntic)
    print('At rest: C across at %.1f links, round trip %.2f steps, c = %.5f per step' % (TAC0 / 2, TAC0, c))
    print('At rest: B settles at %.4f against %.4f for C (ratio %.4f +- %.4f: the network has no preferred directions)'
          % (L0, Lp, L0 / Lp, eL0 / Lp), flush=True)
    print(' v/c | internal time (round trip A-C), slowdown | equilibrium distance A-B relative to rest | relativity | difference')
    for b in velocita:
        g = 1 / math.sqrt(1 - b * b)
        if (Lp / g) * 1.1 / (1 - b) >= 0.45:
            print(' %.2f | beyond the box limit (the flash would go around the closed universe)' % b, flush=True)
            continue
        TAC, eAC = B.andata_ritorno('traverso', Lp, b * c, ntic, int(1.6 * Lp * g / c) + 40)
        L, eL, mis = equilibrio(B, b, c, Lp * TAC0 / TAC, TAC, eAC, ntic)
        r = L / L0
        er = r * math.sqrt((eL / L) ** 2 + (eL0 / L0) ** 2)
        print(' %.2f |                 %.3f (relativity %.3f)          |              %.4f +- %.4f              |    %.4f    |  %+.1f%%'
              % (b, TAC / TAC0, g, r, er, 1 / g, 100 * (r * g - 1)), flush=True)
    print('time %.0f s' % (time.time() - t0))
