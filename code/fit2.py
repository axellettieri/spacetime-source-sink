import os
import numpy as np, json, math
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'fit.py')).read().split("models={")[0]
exec(src)
from scipy.integrate import solve_ivp
# Variante "crunch": lo scarico non dipende dal livello della vasca.
#  kind='cum'  : scarico proporzionale alla massa compatta accumulata M*(t)=integrale di psi
#  kind='const': scarico a ritmo costante
def make(n,mm,kind,zmax=30,N=240):
    u=np.linspace(np.log(1+zmax),0,N+1)
    def rhs(u_,y,Om,al,ga):
        x,M=y; z=math.expm1(u_); E2=Om*(1+z)**3+Or*(1+z)**4+x
        E=math.sqrt(max(E2,1e-8)); ps=(1+z)**2.7/(1+((1+z)/2.9)**5.6)
        sink=ga*(M if kind=='cum' else 1.0)*(1+z)**mm
        return (-(al*ps*(1+z)**n-sink)/E, -ps/E)
    def solve(Om,ga,al):
        x=0.0;M=0.0; ys=[(x,M)]; ul=u.tolist()
        for i in range(N):
            hh=ul[i+1]-ul[i]
            a1,b1=rhs(ul[i],(x,M),Om,al,ga);a2,b2=rhs(ul[i]+hh/2,(x+hh*a1/2,M+hh*b1/2),Om,al,ga)
            a3,b3=rhs(ul[i]+hh/2,(x+hh*a2/2,M+hh*b2/2),Om,al,ga);a4,b4=rhs(ul[i+1],(x+hh*a3,M+hh*b3),Om,al,ga)
            x+=hh*(a1+2*a2+2*a3+a4)/6; M+=hh*(b1+2*b2+2*b3+b4)/6; ys.append((x,M))
        return np.array(ys)
    def model(p):
        Om,ga=p; tgt=1-Om-Or
        g=lambda la: solve(Om,ga,np.exp(la))[-1,0]-tgt
        la=brentq(g,-40,15,xtol=1e-10); ys=solve(Om,ga,np.exp(la))
        it=interp1d(u[::-1],ys[::-1,0],kind='cubic')
        fn=lambda z: it(np.log1p(z)); fn.alpha=np.exp(la); fn.M0=ys[-1,1]; fn.xmin=ys[:,0].min()
        return fn
    def future(Om,ga,al,M0,tmax=400):
        # tempo in unita' 1/H0; y=(ln a, x, M)
        def f(t,y):
            la,x,M=y; a=np.exp(la); z=1/a-1
            E=np.sqrt(max(Om/a**3+x,0))
            sink=ga*(M if kind=='cum' else 1.0)*(1+z)**mm
            return [E, al*psi(z)*(1+z)**n-sink, psi(z)]
        ev=lambda t,y: Om/np.exp(3*y[0])+y[1]-1e-7
        ev.terminal=True; ev.direction=-1
        s=solve_ivp(f,[0,tmax],[0,1-Om-Or,M0],events=ev,rtol=1e-8,atol=1e-10,max_step=0.05)
        tz=None
        # istante in cui l'energia oscura passa per zero
        xx=s.y[1]; i=np.where(xx<0)[0]
        if len(i): tz=float(s.t[i[0]])
        return (float(s.t_events[0][0]) if len(s.t_events[0]) else None, tz, float(np.exp(s.y[0,-1])))
    model.future=future
    return model
res={}
for kind in ('cum','const'):
  for n in (0,3):
    for mm in (0,3):
        mk=make(n,mm,kind); best=None
        for g0 in (0.03,0.5):
            def F(v):
                if v[0]<0.05 or v[0]>0.6 or v[1]<=0: return 1e9
                try: return chi2(mk(v[:2]),v[0],v[2])
                except Exception: return 1e9
            r=minimize(F,[0.3,g0,30.0],method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-3,maxiter=3000))
            if best is None or r.fun<best.fun: best=r
        v=best.x; xf=mk(v[:2]); cc=chi2(xf,v[0],v[2],True)
        if not isinstance(cc,tuple): print(kind,n,mm,'FIT FALLITO'); continue
        cb,cs=cc
        tt,tz,af=mk.future(v[0],v[1],xf.alpha,xf.M0)
        zz=np.array([0,0.25,0.5,1,1.5,2,2.33])
        res[f'{kind}_n{n}_m{mm}']=dict(Om=float(v[0]),gamma=float(v[1]),alpha=float(xf.alpha),chi2=float(best.fun),bao=float(cb),sn=float(cs),xmin_past=float(xf.xmin),
             t_turn_H0=tt,t_x0_H0=tz,a_end=af,x_ratio=[float(a) for a in xf(zz)/xf(0.0)])
        print(kind,n,mm,res[f'{kind}_n{n}_m{mm}'],flush=True)
json.dump(res,open('results2.json','w'),indent=1)
