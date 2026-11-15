import os
"""Shell test: does the expansion rate measured by supernovae change with distance?"""
import numpy as np, math
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'fit.py')).read().split("models={")[0]; exec(src)
from scipy.optimize import minimize_scalar
def base(Om):
    E=np.sqrt(Om*(1+zg)**3+(1-Om)); I=np.concatenate([[0],np.cumsum(np.diff(zg)*(1/E[1:]+1/E[:-1])/2)])
    return mbs-5*np.log10((1+zsh)*interp1d(zg,I,kind='cubic')(zs))
def gls(A,d):
    F=A.T@Ci@A; cov=np.linalg.inv(F); b=cov@(A.T@Ci@d); r=d-A@b; return b,cov,float(r@Ci@r)
K=math.log(10)/5
edges=[0.01,0.023,0.04,0.06,0.1,0.15,0.3,0.5,2.3]
for Om in (0.30,0.334):
    d=base(Om); A=np.array([(zs>=a)&(zs<b) for a,b in zip(edges[:-1],edges[1:])]).T.astype(float)
    b,cov,c2=gls(A,d)
    # reference: weighted mean of the shells beyond z=0.15
    ref=[5,6,7]; w=np.linalg.inv(cov[np.ix_(ref,ref)]); wv=w.sum(0)/w.sum(); 
    print(f'--- Om={Om}  chi2={c2:.1f}')
    for i,(a,bb) in enumerate(zip(edges[:-1],edges[1:])):
        L=np.zeros(len(b)); L[i]=1; L[ref]-=wv
        dm=L@b; e=math.sqrt(L@cov@L)
        print(f'z {a:.3f}-{bb:.3f}  N={int(A[:,i].sum()):4d}  dH/H = {-K*dm*100:+.2f} +- {K*e*100:.2f} %' + ('  (reference)' if i in ref else ''))
# step test: a single jump at z_edge, Om free
print('--- single step (Om free)')
for ze in (0.04,0.07,0.1,0.15,0.2,0.3):
    def f(Om):
        A=np.array([np.ones_like(zs),(zs<ze).astype(float)]).T; return gls(A,base(Om))[2]
    r=minimize_scalar(f,bounds=(0.1,0.6),method='bounded'); A=np.array([np.ones_like(zs),(zs<ze).astype(float)]).T
    b,cov,c2=gls(A,base(r.x)); c0=minimize_scalar(lambda Om:gls(np.ones((len(zs),1)),base(Om))[2],bounds=(0.1,0.6),method='bounded').fun
    print(f'edge z={ze:.2f}: H inside - H outside = {-K*b[1]*100:+.2f} +- {K*math.sqrt(cov[1,1])*100:.2f} %   Om={r.x:.3f}  dchi2={c2-c0:+.2f}')
print('expected to explain the tension: +8 %  (0.167 mag)')
