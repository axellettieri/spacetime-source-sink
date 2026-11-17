"""Experiment 9: entropy and an internal arrow of time along the cycles of experiment 5.
The dynamics is exactly that of esperimento5.py (same rules, same random numbers): this script only adds measurements.

ENTROPY (counts of configurations compatible with what we see, in nats):
  S_geo  = V * H_deg * ln 2   geometry: V vertices, H_deg = Shannon entropy (bits) of the distribution of vertex degrees.
                              A closed triangulated sphere always has E = 3V - 6 edges, so edges add no information;
                              the degrees say how irregular the network is. The 4-vertex start (all degrees 3) has S_geo = 0.
  S_pos  = N_obj * ln V       matter: number of ways to place the N_obj objects on the V vertices.
  S      = S_geo + S_pos      total entropy of the graph used in the summary.
  S_bh   = sum of m^2         comparison only: what a black-hole-like entropy (Bekenstein-Hawking grows as mass^2) would give.
Both are cheap: one pass over the vertices and the objects, no random numbers (so the dynamics is untouched).

INTERNAL ARROW (does not use the step counter):
  Every vertex created by an event keeps a RECORD: the vertex where the event happened (its parent).
  This bookkeeping is added from outside and does not change the dynamics.
  Between two snapshots X and Y the script counts
    rec(Y<-X) = vertices alive in Y with an ancestor that is alive in X and gone in Y   (Y remembers something X had and lost)
    rec(X<-Y) = the same with X and Y exchanged.
  The snapshot holding records of the other is the later one. The script compares this verdict with the step counter.
  For contrast it also computes a naive SIZE arrow: which snapshot has more new descendants of the vertices they share.
MEMORY: fraction of the vertices of a reference graph that the current graph still records (alive, or ancestor of a
  living vertex). References: the maximum of the previous cycle (across each bounce) and the largest graph of the first
  10% of the run (to see whether the past is forgotten also without a bounce).

Usage: python esperimento9.py A|B [steps] [seed]
   A = singularity decays only at the minimum (TC=0)  -> cycles
   B = hot objects decay earlier (TC=0.9)              -> irregular equilibrium
Writes esperimento9_<A|B>.csv (one row per sample) and prints a summary."""
import sys, math, time
from collections import Counter
import esperimento4, esperimento5

PAR = {}                       # record: vertex -> vertex where the event that created it happened
KIDS = {}                      # inverse of PAR
_split = esperimento4.Sfera.split
def split_rec(R, a, b):
    v = _split(R, a, b); PAR[v] = a; KIDS.setdefault(a, []).append(v); return v
esperimento4.Sfera.split = split_rec

def ancestors_closure(alive):
    """alive vertices plus all their ancestors (records), alive or not"""
    seen = set()
    for v in alive:
        while v not in seen:
            seen.add(v)
            if v not in PAR: break
            v = PAR[v]
    return seen

def descendants_alive(roots, alive):
    """number of alive vertices descending (strictly) from the roots"""
    n = 0; st = [k for r in roots for k in KIDS.get(r, ())]; seen = set()
    while st:
        v = st.pop()
        if v in seen: continue
        seen.add(v)
        if v in alive: n += 1
        st.extend(KIDS.get(v, ()))
    return n

ROWS = []
STATE = dict(prev=None, cur_max=None, prev_max=None, vmax_cycle=0, cycle=1, k=0, prevN=None, first=None, vfirst=0, tfirst=0)
BOUNCES = []      # (index of the last sample before the bounce, index of the first after)
def osserva(t, R, objs):
    V = len(R.nbr); alive = set(R.nbr)
    deg = Counter(len(n) for n in R.nbr.values()); H = -sum(c / V * math.log2(c / V) for c in deg.values())
    Sgeo = V * H * math.log(2); Spos = len(objs) * math.log(V); Sbh = sum(m * m for _, m, _ in objs)
    mass = Counter(m for _, m, _ in objs); N = len(objs); Hm = -sum(c / N * math.log2(c / N) for c in mass.values())
    row = dict(t=t, V=V, E=3 * V - 6, Hdeg=H, Sgeo=Sgeo, Spos=Spos, S=Sgeo + Spos, Sbh=Sbh, Nobj=N, Hmass=Hm,
               rec_fw='', rec_bw='', size_fw='', size_bw='', mem_prev='', mem_first='', cycle=STATE['cycle'])
    # records arrow with the previous sample (X = previous by counter, Y = now): uses only the records
    if STATE['prev'] is not None:
        X = STATE['prev']; gone_X = X - alive; gone_Y = alive - X; common = X & alive
        rec_fw = descendants_alive(gone_X, alive)      # Y remembers vertices that X had and Y lost
        rec_bw = descendants_alive(gone_Y, X)          # X remembers vertices that Y has and X lacks (impossible if Y is later)
        size_fw = descendants_alive(common, gone_Y)    # new descendants of shared vertices that only Y has
        size_bw = descendants_alive(common, gone_X)    # descendants of shared vertices that only X has
        row.update(rec_fw=rec_fw, rec_bw=rec_bw, size_fw=size_fw, size_bw=size_bw)
    # BIG BANG (for the report only): the singularity decays into single particles, so the number of objects jumps
    if STATE['prevN'] is not None and STATE['prevN'] <= 3 and N >= 50:
        STATE['prev_max'] = STATE['cur_max']; STATE['cycle'] += 1; STATE['vmax_cycle'] = 0
        BOUNCES.append((len(ROWS) - 1, len(ROWS)))
    row['cycle'] = STATE['cycle']
    if V > STATE['vmax_cycle']: STATE['vmax_cycle'] = V; STATE['cur_max'] = alive
    # reference for the memory of the distant past: the largest graph in the first 10% of the run
    if t <= STATE['tfirst'] and V > STATE['vfirst']: STATE['vfirst'] = V; STATE['first'] = alive; STATE['tfirst_seen'] = t
    if STATE['k'] % 5 == 0 or (BOUNCES and BOUNCES[-1][1] == len(ROWS)):
        A = ancestors_closure(alive)
        if STATE['prev_max'] is not None: P = STATE['prev_max']; row['mem_prev'] = len(P & A) / len(P)
        if STATE['first'] is not None and t > STATE['tfirst']: F = STATE['first']; row['mem_first'] = len(F & A) / len(F)
    STATE['prevN'] = N
    STATE['prev'] = alive; STATE['k'] += 1
    ROWS.append(row)

def report(name):
    keys = list(ROWS[0])
    with open('esperimento9_%s.csv' % name, 'w') as f:
        f.write(';'.join(keys) + '\n')
        for r in ROWS: f.write(';'.join(('%.4g' % r[k]) if isinstance(r[k], float) else str(r[k]) for k in keys) + '\n')
    print('\n=== ENTROPY ALONG THE CYCLES ===')
    Vs = [r['V'] for r in ROWS]
    print(' (S = geometry + matter, in nats; S_bh = what a black-hole-like entropy, growing as mass^2, would give)')
    print(' cycle | big bang: step     V       S | maximum: step     V       S | last sample before the next bounce: step   V     S    objects  S_bh')
    cyc = sorted(set(r['cycle'] for r in ROWS))
    for c in cyc:
        idx = [i for i, r in enumerate(ROWS) if r['cycle'] == c]; i0 = idx[0]; imax = max(idx, key=lambda i: ROWS[i]['V']); il = idx[-1]
        r0, rm, rl = ROWS[i0], ROWS[imax], ROWS[il]
        last = ' %7d %5d %7.1f   %4d  %6d' % (rl['t'], rl['V'], rl['S'], rl['Nobj'], rl['Sbh']) if c != cyc[-1] else '   (run ends inside this cycle)'
        print(' %5d | %7d %6d %8.1f | %7d %6d %8.1f |%s' % (c, r0['t'], r0['V'], r0['S'], rm['t'], rm['V'], rm['S'], last))
    if BOUNCES:
        print(' entropy just before / just after each bounce (the bounce itself lasts a few steps, between two samples):')
        for i, j in BOUNCES:
            print('   step %7d: V %4d, S %7.1f, objects %3d, S_bh %6d  ->  step %7d: V %5d, S %7.1f, objects %3d, S_bh %4d' % (
                ROWS[i]['t'], ROWS[i]['V'], ROWS[i]['S'], ROWS[i]['Nobj'], ROWS[i]['Sbh'], ROWS[j]['t'], ROWS[j]['V'], ROWS[j]['S'], ROWS[j]['Nobj'], ROWS[j]['Sbh']))
    else:
        print(' no bounce (no decay of a singularity) in this run')
    # correlation between entropy and size
    S = [r['S'] for r in ROWS]; mV = sum(Vs) / len(Vs); mS = sum(S) / len(S)
    cov = sum((a - mV) * (b - mS) for a, b in zip(Vs, S)); cV = sum((a - mV) ** 2 for a in Vs); cS = sum((b - mS) ** 2 for b in S)
    print(' correlation between entropy S and number of vertices V: %.3f' % (cov / math.sqrt(cV * cS) if cV and cS else float('nan')))
    print('\n=== INTERNAL ARROW (records) COMPARED WITH THE STEP COUNTER ===')
    ph = {'expansion': [0, 0, 0, 0], 'collapse': [0, 0, 0, 0], 'still': [0, 0, 0, 0]}   # agree, disagree, undetermined, size-arrow agree
    bw_tot = 0
    for a in range(1, len(ROWS)):
        r = ROWS[a]; dV = ROWS[a]['V'] - ROWS[a - 1]['V']
        p = 'expansion' if dV > 0.002 * max(1, ROWS[a - 1]['V']) else 'collapse' if dV < -0.002 * max(1, ROWS[a - 1]['V']) else 'still'
        bw_tot += r['rec_bw']
        if r['rec_fw'] > r['rec_bw']: ph[p][0] += 1
        elif r['rec_fw'] < r['rec_bw']: ph[p][1] += 1
        else: ph[p][2] += 1
        if r['size_fw'] > r['size_bw']: ph[p][3] += 1
    print(' phase (by counter) | pairs | records arrow agrees | disagrees | undetermined | naive size arrow agrees')
    for p, (ag, dis, und, sz) in ph.items():
        n = ag + dis + und
        if n: print(' %-18s | %5d | %6.1f%%              | %5.1f%%    | %5.1f%%       | %5.1f%%' % (p, n, 100 * ag / n, 100 * dis / n, 100 * und / n, 100 * sz / n))
    print(' total records pointing the "wrong" way (must be 0 by causality): %d' % bw_tot)
    print('\n=== MEMORY OF THE PAST (fraction of the vertices of a reference graph that the current graph still records) ===')
    if BOUNCES:
        for i, j in BOUNCES:
            r = ROWS[j]; ra = ROWS[j]['rec_fw']
            print(' bounce at step ~%d: memory of the previous cycle maximum %.4f; records linking the two sides of the bounce: %s vertices' % (r['t'], r['mem_prev'] if r['mem_prev'] != '' else float('nan'), ra))
    mf = [(r['t'], r['mem_first']) for r in ROWS if r['mem_first'] != '']
    if mf:
        print(' memory of the largest graph of the first 10%% of the run (step %d, %d vertices):' % (STATE['tfirst_seen'], STATE['vfirst']))
        n = len(mf)
        for q in (0, n // 8, n // 4, n // 2, 3 * n // 4, n - 1):
            print('   at step %7d: %.4f' % mf[q])
    print('\nwritten esperimento9_%s.csv' % name)

if __name__ == '__main__':
    a = sys.argv[1:]; ver = (a[0] if a else 'A').upper(); passi = int(a[1]) if len(a) > 1 else 450000; seme = int(a[2]) if len(a) > 2 else 1
    TC = 0.0 if ver == 'A' else 0.9
    t0 = time.time(); STATE['tfirst'] = passi // 10
    esperimento5.esegui(300, 0.02, 0.9, 0.05, 0.005, passi, TC, 0.01, seme, osserva=osserva)
    report(ver); print('time %.0f s' % (time.time() - t0))
