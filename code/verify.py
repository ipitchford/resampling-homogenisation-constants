#!/usr/bin/env python3
"""Independent finite checks. These checks do not prove the unrestricted theorems."""


class CheckFailure(Exception):
    """Raised when a verification check fails (works under python -O)."""


def check(condition, label="check"):
    if not condition:
        raise CheckFailure(label)

from fractions import Fraction as Q
from itertools import product
from math import factorial, log
from pathlib import Path
import json, platform, time
import numpy as np
import scipy
from scipy.optimize import brentq, minimize
import sympy as sp
BASE=Path(__file__).resolve().parents[1]

def compositions(n,k):
    if k==1:
        yield (n,); return
    for j in range(n+1):
        for rest in compositions(n-j,k-1): yield (j,)+rest

def tv_product(P,Qv):
    ans=Q(0)
    for state in product(*[range(len(p)) for p in P]):
        a=b=Q(1)
        for i,j in enumerate(state): a*=P[i][j];b*=Qv[i][j]
        ans+=abs(a-b)
    return ans/2

def occupancy(n):
    return [(k,l,Q(factorial(n),factorial(k)*factorial(l)*factorial(n-k-l))*Q(1,n)**(k+l)*Q(n-2,n)**(n-k-l)) for k in range(n+1) for l in range(n-k+1)]

def Fn(n,d):
    return sum(w*(1-(1-d)**max(k,l)) for k,l,w in occupancy(n))

def extremiser(n,d):
    # Five distinguishable outcomes a0,a1,b0,b1,c; inactive labels coalesce safely.
    zero=Q(0);one=Q(1)
    p=[[one,zero,zero,zero,zero],[zero,zero,1-d,d,zero]]
    q=[[1-d,d,zero,zero,zero],[zero,zero,one,zero,zero]]
    for _ in range(n-2):p.append([zero,zero,zero,zero,one]);q.append([zero,zero,zero,zero,one])
    return p,q

def poly_modulus(d):
    r=1-d
    if d==0:return (0.,1.,0.)
    if d==1:return (1.,.6,1/3)
    r0=brentq(lambda z:z**3+3*z*z-z-1,.5,.8)
    r1=(13**.5-3)/2
    if r>=r0-1e-13: x=1.
    else:
        u=brentq(lambda u:u**6+2*u**5+r*u**4-r*r*u*u-r**3,r,r**.5,xtol=1e-14)
        x=log(u)/log(r)
    if r>=r1-1e-13:y=0.
    else:
        v=brentq(lambda v:v**4-r*v*v-2*r*v-r*r,r**.5,1.,xtol=1e-14)
        y=log(v)/log(r)
    f=r**(3*x)+r**(3-3*x)+3*r**(1+x)+3*r**(2*x)+3*r**(2-x)
    g=r**(3-3*y)+3*r**(1+y)+3*r**(2-y)+3*r**(2-2*y)
    return 1-(f+g+6*r)/27,x,y

def direct_overlap3(a,b):
    ans=0.
    for k in compositions(3,3):
        w=factorial(3)/np.prod([factorial(i) for i in k])/27
        ans+=w*np.exp(-max(np.dot(k,a),np.dot(k,b)))
    return ans

def main():
    start=time.perf_counter(); report={"scope":"Exact finite arithmetic, symbolic identities and independent floating-point optimisation; no proof-assistant verification.","python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__}
    exact=[]
    for n in range(2,6):
        for d in [Q(1,100),Q(1,10),Q(1,2),Q(9,10),Q(1)]:
            p,q=extremiser(n,d)
            dp=tv_product(p,q)
            bp=[sum(col)/n for col in zip(*p)];bq=[sum(col)/n for col in zip(*q)]
            h=tv_product([bp]*n,[bq]*n)
            check(dp==d and h==Fn(n,d))
            exact.append({"n":n,"D":str(d),"H_direct":str(h),"ratio":str(h/d)})
    report['direct_rational_product_checks']=exact
    cns=[]
    for n in range(2,31):
        c=sum(w*max(k,l) for k,l,w in occupancy(n));check(c<=Q(14,9))
        cns.append({"n":n,"C":str(c),"decimal":float(c)})
    check(all(Q(cns[i]['C'])>Q(cns[i+1]['C']) for i in range(1,len(cns)-1)))
    report['occupancy_constants']=cns
    z,x,y=sp.symbols('r x y', positive=True)
    d=sp.symbols('d')
    f2=sp.Rational(3,2)*d-sp.Rational(1,2)*d*d
    f3=sp.Rational(14,9)*d-sp.Rational(2,3)*d*d+sp.Rational(2,27)*d**3
    check(sp.expand(f3-f2-d*(4*d*d-9*d+3)/54)==0)
    n3=sum(sp.Rational(w.numerator,w.denominator)*(1-(1-d)**max(k,l)) for k,l,w in occupancy(3))
    check(sp.expand(n3-f3)==0)
    f=z**(3*x)+z**(3-3*x)+3*z**(1+x)+3*z**(2*x)+3*z**(2-x)
    g=z**(3-3*y)+3*z**(1+y)+3*z**(2-y)+3*z**(2-2*y)
    low=sp.expand((1-(f.subs(x,1)+g.subs(y,0)+6*z)/27).subs(z,1-d))
    check(sp.simplify(low-f3)==0)
    report['symbolic_checks']=['F3 exact multinomial polynomial','F3-F2 crossover factorisation','full BAEC3 sparse regime equals F3']
    numeric=[]; rng=np.random.default_rng(20261115)
    for dd in [.001,.1,.3248691294333441,.4,.6,.6972243622680054,.8,.95,.999]:
        h,xx,yy=poly_modulus(dd);m=-np.log1p(-dd)
        a=m*np.array([xx,1-xx,0]);b=m*np.array([0,yy,1-yy])
        hd=1-direct_overlap3(a,b)
        check(abs(h-hd)<1e-11 and xx+yy<=1+1e-10 and 2*xx+yy>=1-1e-10)
        # Independently minimise the convex chamber in four coordinates.
        def objective(v):
            aa=m*np.array([v[0],v[1],1-v[0]-v[1]])
            bb=m*np.array([v[2],v[3],1-v[2]-v[3]])
            return direct_overlap3(aa,bb)
        constraints=[{'type':'ineq','fun':lambda v:1-v[0]-v[1]}, {'type':'ineq','fun':lambda v:1-v[2]-v[3]}, {'type':'ineq','fun':lambda v:v[1]-v[3]}, {'type':'ineq','fun':lambda v:v[0]-v[2]-v[1]+v[3]}]
        opt=minimize(objective,[.7,.2,.1,.1],method='SLSQP',bounds=[(0,1)]*4,constraints=constraints,options={'ftol':1e-13,'maxiter':2000})
        check(opt.success and abs((1-opt.fun)-h)<2e-8)
        # Broad general BAEC samples: pointwise finite check, not a global optimisation certificate.
        worst=-1.
        for _ in range(200):
            aa=m*rng.dirichlet(np.ones(3));bb=m*rng.dirichlet(np.ones(3))
            hr=1-direct_overlap3(aa,bb);worst=max(worst,hr)
            check(hr<=h+1e-10)
        numeric.append({'D':dd,'H_exact_formula_numeric':h,'x':xx,'y':yy,'H_independent_chamber_optimisation':1-opt.fun,'largest_H_in_200_random_pairs':worst})
    report['BAEC3_numerical_checks']=numeric
    # Scope check: removing source labels from the opposite-erasure construction makes the averages identical.
    dd=Q(1,10)
    pp,qq=extremiser(3,dd)
    def collapse(rows):
        return [[row[0]+row[2]+row[4],row[1]+row[3]] for row in rows]
    pp=collapse(pp);qq=collapse(qq)
    bp=[sum(col)/3 for col in zip(*pp)];bq=[sum(col)/3 for col in zip(*qq)]
    check(bp==bq and bp==[1-dd/3,dd/3])
    check(tv_product([bp]*3,[bq]*3)==0)
    report['scope_check']='Collapsing all common symbols and both rare symbols gives identical homogeneous Bernoulli laws (H=0); tagged construction does not refute the binary-alphabet conjecture.'
    report['elapsed_seconds']=time.perf_counter()-start;report['status']='All listed checks passed'
    (BASE/'logs'/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['direct_rational_product_checks','occupancy_constants','BAEC3_numerical_checks']},indent=2))
    print('Exact direct-product cases:',len(exact),'occupancy sizes:',len(cns),'BAEC3 comparisons:',len(numeric),'random BAEC pairs:',len(numeric)*200)
if __name__=='__main__':main()
