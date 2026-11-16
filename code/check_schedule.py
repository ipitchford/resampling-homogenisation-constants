"""Finite exact checks of schedule theorem, not universal proof verification."""
from fractions import Fraction as F
from itertools import product
from math import prod
import json

def moments(law,n):
    p=sum(w for N,w in law if any(N))
    c=max(sum(w*max(N[i],N[j]) for N,w in law) for i in range(n) for j in range(n))
    return p,c,max(p,c)

def check(law,n,grid):
    p,c,k=moments(law,n)
    count=0
    for values in product(grid,repeat=2*n):
        r=values[:n];s=values[n:]
        d=1-min(prod(r),prod(s))
        h=sum(w*(1-min(prod(r[i]**N[i] for i in range(n)),prod(s[i]**N[i] for i in range(n)))) for N,w in law)
        if h>k*d: raise RuntimeError(('violated',law,r,s,d,h,k))
        count+=1
    return {'p':str(p),'C':str(c),'K':str(k),'cases':count}

laws={
    'empty': (2,[((0,0),F(1))]),
    'deterministic_baseline':(2,[((1,1),F(1))]),
    'bursty':(2,[((0,0),F(1,2)),((4,0),F(1,4)),((0,2),F(1,4))]),
    'single_uniform_n3':(3,[((1,0,0),F(1,3)),((0,1,0),F(1,3)),((0,0,1),F(1,3))]),
    'two_opposing_bursts':(2,[((3,0),F(1,2)),((0,3),F(1,2))]),
}
out={name:check(law,n,[F(1,4),F(1,2),F(1)]) for name,(n,law) in laws.items()}
# Necessity of the p term when C<p.
r=F(1,100)
d=1-r**3;h=1-r
if not h>F(2,3)*d: raise RuntimeError('p-term witness failed')
out['p_term_witness']={'D':str(d),'H':str(h),'H_over_D':str(h/d),'C_only':str(F(2,3))}
print(json.dumps({'status':'passed','tests':out},indent=2))
