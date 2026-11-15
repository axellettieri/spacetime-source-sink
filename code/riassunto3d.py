"""Riassunto delle corse 3D: python riassunto3d.py   (legge tutti i .json prodotti da universo3d.py nella cartella)"""
import json,glob,numpy as np
for f in sorted(glob.glob('*.json')):
    try: J=json.load(open(f)); h=np.array(J['storia']); P=J['par']
    except Exception: continue
    esp=h[:,1]>0                                                    # 'oggi' = passaggio per la dimensione attuale IN ESPANSIONE (non al ritorno)
    i=int(np.argmin(np.where(esp,np.abs(h[:,0]-1),np.inf)))
    print('\n%s | sparpagliamento %s Mpc | alpha %.4g gamma %.4g soglia %.1f | fine: %s (a=%.2f)'%(f,P['lam'],P['alpha'],P['gamma'],P['dc'],J['fine'],h[-1,0]))
    print('   oggi: H0 %.1f km/s/Mpc, spazio creato %.3f, materia %.3f, frazione collassata %.2f, disuniformita\' dello spazio creato %.0f%%'%(100*h[i,1],h[i,2],h[i,3],h[i,4],100*h[i,7]/max(abs(h[i,2]),1e-9)))
    for k,v in J['hubble_locale'].items():
        print('   espansione locale ('+k+'): '+' | '.join('entro %s Mpc %+.2f%% (scarto fra osservatori %.2f%%)'%(R,d['media_pct'],d['scarto_pct']) for R,d in v.items()))
    if J.get('tempo_Gyr'):
        t=np.array(J['tempo_Gyr']); j=int(np.argmax(h[:,0])); print('   eta\' oggi %.1f miliardi di anni; dimensione massima %.2f volte l\'attuale fra %.1f miliardi di anni; fine corsa fra %.1f miliardi di anni (a=%.2f)'%(t[i],h[j,0],t[j]-t[i],t[-1]-t[i],h[-1,0]))
