"""Summary of the 3D runs: python riassunto3d.py   (reads all .json files produced by universo3d.py in the folder)"""
import json,glob,numpy as np
for f in sorted(glob.glob('*.json')):
    try: J=json.load(open(f)); h=np.array(J['storia']); P=J['par']
    except Exception: continue
    esp=h[:,1]>0                                                    # 'today' = crossing of the present size DURING EXPANSION (not on the way back)
    i=int(np.argmin(np.where(esp,np.abs(h[:,0]-1),np.inf)))
    print('\n%s | spreading %s Mpc | alpha %.4g gamma %.4g threshold %.1f | end: %s (a=%.2f)'%(f,P['lam'],P['alpha'],P['gamma'],P['dc'],J['fine'],h[-1,0]))
    print('   today: H0 %.1f km/s/Mpc, created space %.3f, matter %.3f, collapsed fraction %.2f, non-uniformity of created space %.0f%%'%(100*h[i,1],h[i,2],h[i,3],h[i,4],100*h[i,7]/max(abs(h[i,2]),1e-9)))
    for k,v in J['hubble_locale'].items():
        print('   local expansion ('+k+'): '+' | '.join('within %s Mpc %+.2f%% (spread between observers %.2f%%)'%(R,d['media_pct'],d['scarto_pct']) for R,d in v.items()))
    if J.get('tempo_Gyr'):
        t=np.array(J['tempo_Gyr']); j=int(np.argmax(h[:,0])); print('   age today %.1f billion years; maximum size %.2f times the present one in %.1f billion years; end of run in %.1f billion years (a=%.2f)'%(t[i],h[j,0],t[j]-t[i],t[-1]-t[i],h[-1,0]))
