"""Exact direct check of the example obstructing a common-part reduction.

This checks a finite rational example; the general obstruction is proved in
evidence/unrestricted_boundary.md and the manuscript.
"""


class CheckFailure(Exception):
    """Raised when a verification check fails (works under python -O)."""


def check(condition, label="check"):
    if not condition:
        raise CheckFailure(label)

from fractions import Fraction as Q
from itertools import product
from math import prod
from pathlib import Path
import json


def tensor_tv(ps, qs):
    return sum(abs(prod(ps[i][x] for i, x in enumerate(xs))
                   - prod(qs[i][x] for i, x in enumerate(xs)))
               for xs in product(*(range(len(p)) for p in ps))) / 2


def average(rows):
    return [sum(row[j] for row in rows) / len(rows)
            for j in range(len(rows[0]))]


def main():
    ps = [[Q(1, 5), Q(4, 5)], [Q(1, 2), Q(1, 2)]]
    qs = [[Q(0), Q(1)], [Q(2, 5), Q(3, 5)]]
    d = tensor_tv(ps, qs)
    pb, qb = average(ps), average(qs)
    h = tensor_tv([pb, pb], [qb, qb])
    check((d, h, h/d) == (Q(1, 5), Q(87, 400), Q(87, 80)))

    # Each lifted row has its own common, P-private, and Q-private symbol.
    common = [1-tensor_tv([p], [q]) for p, q in zip(ps, qs)]
    lifted_p, lifted_q = [], []
    for i, r in enumerate(common):
        p, q = [Q(0)] * 6, [Q(0)] * 6
        p[3*i] = q[3*i] = r
        p[3*i+1] = 1-r
        q[3*i+2] = 1-r
        lifted_p.append(p)
        lifted_q.append(q)
    dp = tensor_tv(lifted_p, lifted_q)
    lpb, lqb = average(lifted_p), average(lifted_q)
    hp = tensor_tv([lpb, lpb], [lqb, lqb])
    check((dp, hp, hp/dp) == (Q(7, 25), Q(111, 400), Q(111, 112)))
    check(dp > d and hp > h and hp/dp < h/d)
    result = {
        'status': 'passed',
        'scope': 'One exact rational finite-product example; no universal bound inferred.',
        'original': dict(D=str(d), H=str(h), H_over_D=str(h/d)),
        'common_part_refinement': dict(D=str(dp), H=str(hp), H_over_D=str(hp/dp))
    }
    (Path(__file__).resolve().parents[1] / 'logs' / 'boundary.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
