#!/usr/bin/env python3
"""Exact finite checks for the unrestricted asymptotic witness; not a universal proof."""


class CheckFailure(Exception):
    """Raised when a verification check fails (works under python -O)."""


def check(condition, label="check"):
    if not condition:
        raise CheckFailure(label)

from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path
import json

BASE=Path(__file__).resolve().parents[1]

def compositions(total,dim):
    if dim==1:
        yield (total,)
        return
    for a in range(total+1):
        for rest in compositions(total-a,dim-1):
            yield (a,)+rest

def overlap_tensor(ps,qs):
    result=F(0)
    for state in product(*[range(len(p)) for p in ps]):
        a=b=F(1)
        for i,j in enumerate(state):
            a*=ps[i][j]
            b*=qs[i][j]
        result+=min(a,b)
    return result

def rational_string(x):
    return str(x)

rows=[]
q_zero_checks=[]
# (n,k,B): B=kn+1 is the radix of the first version; any integer B>=k+1 is admissible (Theorem 7.1).
for n,k,radix in [(1,2,3),(2,1,3),(2,2,5),(3,1,4),(2,2,3),(3,1,2),(3,1,3),(4,1,2),(3,2,3)]:
    total=n*k
    weights=[radix**i for i in range(n)]
    S=sum(weights); T=n*n*S
    A=[(n-1)*S+n*w for w in weights]
    B=[(n+1)*S-n*w for w in weights]
    check(all(e>0 for e in A+B) and sum(A)==sum(B)==T)
    r=F(1,2)**T
    ps=[];qs=[]
    for i in range(n):
        p=[F(0)]*(n+2);q=[F(0)]*(n+2)
        p[i]=F(1,2)**A[i];p[n]=1-p[i]
        q[i]=F(1,2)**B[i];q[n+1]=1-q[i]
        ps.append(p);qs.append(q)
    check(overlap_tensor(ps,qs)==r)
    count_rows=[]
    overlap_count=F(0);weight_sum=F(0);balanced_weight=F(0)
    off_gaps=[]
    for counts in compositions(total,n):
        coeff=F(factorial(total),n**total)
        for c in counts:coeff/=factorial(c)
        ea=sum(c*a for c,a in zip(counts,A))
        eb=sum(c*b for c,b in zip(counts,B))
        exponent=max(ea,eb)
        integer_difference=sum(c*w for c,w in zip(counts,weights))-k*S
        check(exponent==k*T+n*abs(integer_difference))
        balanced=counts==(k,)*n
        check((integer_difference==0)==balanced)
        check(integer_difference % (radix-1)==0)
        if balanced:
            check(exponent==k*T)
            balanced_weight+=coeff
        else:
            check(exponent>=k*T+n*(radix-1))
            off_gaps.append(exponent-k*T)
        overlap_count+=coeff*F(1,2)**exponent
        weight_sum+=coeff
        count_rows.append({'counts':counts,'weight':str(coeff),'exponent':exponent})
    check(weight_sum==1)
    q_expected=F(factorial(total),(factorial(k)**n)*(n**total))
    check(balanced_weight==q_expected)
    lower=q_expected*r**k
    upper=lower+(1-q_expected)*(r**k)*F(1,2)**(n*(radix-1))
    check(lower<=overlap_count<=upper)
    pbar=[sum(col)/n for col in zip(*ps)]
    qbar=[sum(col)/n for col in zip(*qs)]
    overlap_direct=overlap_tensor([pbar]*total,[qbar]*total)
    check(overlap_direct==overlap_count)
    # Every deterministic nonbalanced schedule has q=0 and the same pointwise upper bound.
    check(all(row['exponent']>=k*T+n*(radix-1) for row in count_rows if tuple(row['counts'])!=(k,)*n))
    if (n,k,radix)==(2,2,5):
        # One explicit q=0 revealed schedule: repeat coordinate 1 four times.
        counts=(total,0)
        q_zero_overlap=overlap_tensor([ps[0]]*total,[qs[0]]*total)
        q_zero_upper=(r**k)*F(1,2)**(n*(radix-1))
        check(F(0)<=q_zero_overlap<=q_zero_upper)
        q_zero_checks.append({'n':n,'k':k,'counts':counts,'q':'0',
                              'overlap_direct':str(q_zero_overlap),'upper_bound':str(q_zero_upper),
                              'normalised_overlap':str(q_zero_overlap/r**k),'status':'passed'})
    rows.append({'n':n,'k':k,'radix':radix,'weights':weights,'S':S,'r':str(r),'epsilon':str(F(1,n*S)),'delta':str(F((radix-1)**2,n*(radix**n-1))),
                 'alpha_mass_exponents':A,'beta_mass_exponents':B,
                 'count_vectors_checked':len(count_rows),'balanced_weight':str(q_expected),
                 'minimum_nonbalanced_exponent_gap':min(off_gaps) if off_gaps else None,
                 'overlap_direct':str(overlap_direct),'lower_bound':str(lower),'upper_bound':str(upper),
                 'count_formula_equals_direct_tensor':True,'sandwich_passed':True})
result={'scope':'Exact rational checks of finite constructions; no proof-assistant verification.',
        'cases':rows,'q_zero_schedule_checks':q_zero_checks,'status':'All listed checks passed'}
output=BASE/'logs'/'unrestricted_asymptotic.json'
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'cases':len(rows),'q_zero_schedule_checks':len(q_zero_checks),'count_vectors':sum(row['count_vectors_checked'] for row in rows),'output':str(output)},indent=2))
