"""Independent finite exact-arithmetic checks; not a proof of the universal theorem."""
from fractions import Fraction as F
from itertools import product
from math import factorial, prod
import json

def occupancy(n):
    p=[F(0) for _ in range(n+1)]
    for k in range(n+1):
        for l in range(n-k+1):
            j=n-k-l
            w=F(factorial(n),factorial(k)*factorial(l)*factorial(j))
            w*=F(1,n)**(k+l)*F(n-2,n)**j
            p[max(k,l)]+=w
    if sum(p)!=1: raise RuntimeError('occupancy does not sum to1')
    return p

def curve(n,d):
    return sum(w*(1-(1-d)**m) for m,w in enumerate(occupancy(n)))

def tv_iid(p,q,n):
    return sum(abs(prod(p[j] for j in x)-prod(q[j] for j in x)) for x in product(range(len(p)),repeat=n))/2

rows=[]
for n in range(2,9):
    p=occupancy(n)
    C=sum(j*w for j,w in enumerate(p))
    rows.append({'n':n,'C_exact':str(C),'C_decimal':float(C),'M_pmf':[str(x) for x in p]})
for n in [2,3,4,5]:
    for d in [F(0),F(1,5),F(1,2),F(1)]:
        # Atom order: a0,a1,b0,b1,c.
        p=[F(1,n),F(0),F(1-d,n),F(d,n),F(n-2,n)]
        q=[F(1-d,n),F(d,n),F(1,n),F(0),F(n-2,n)]
        v=tv_iid(p,q,n)
        if v!=curve(n,d): raise RuntimeError(('direct TV mismatch',n,d,v,curve(n,d)))
for d in [F(k,100) for k in range(101)]:
    if curve(2,d)!=F(3,2)*d-F(1,2)*d*d: raise RuntimeError('F2 polynomial')
    if curve(3,d)!=F(14,9)*d-F(2,3)*d*d+F(2,27)*d**3: raise RuntimeError('F3 polynomial')
print(json.dumps({'status':'passed','scope':'16 direct finite-product TV checks,202 polynomial checks,occupancy probabilities for n2..8','rows':rows},indent=2))
