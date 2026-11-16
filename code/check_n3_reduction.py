"""Independent finite and symbolic audit of the n=3 reduction.

This checks algebraic identities and a finite lattice; it does not certify
the unrestricted optimisation theorem or its analytic inequalities.
"""


class CheckFailure(Exception):
    """Raised when a verification check fails (works under python -O)."""


def check(condition, label="check"):
    if not condition:
        raise CheckFailure(label)

from fractions import Fraction as Q
from itertools import product
from math import factorial, cos, pi, sqrt
from pathlib import Path
import json
import sympy as sp
from scipy.integrate import quad


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def overlap(a, b):
    answer = Q(0)
    for counts in compositions(3, 3):
        coefficient = Q(factorial(3), 27)
        for k in counts:
            coefficient /= factorial(k)
        exponent = max(sum(x*k for x, k in zip(a, counts)),
                       sum(x*k for x, k in zip(b, counts)))
        answer += coefficient * Q(1, 2)**exponent
    return answer


def chamber(a, b):
    a, b = list(a), list(b)
    d = [x-y for x, y in zip(a, b)]
    if sum(v < 0 for v in d) > sum(v > 0 for v in d):
        a, b = b, a
        d = [-v for v in d]
    order = sorted(range(3), key=lambda i: d[i], reverse=True)
    a = [a[i] for i in order]
    b = [b[i] for i in order]
    d = [x-y for x, y in zip(a, b)]
    check(d[0] >= d[1] >= 0 >= d[2])
    return a, b


def chamber_expression(a, b):
    m = sum(a)
    terms = [(1, 3*a[0]), (1, 3*a[1]),
             (3, 2*a[0]+a[1]), (3, 2*a[0]+a[2]),
             (3, a[0]+2*a[1]), (1, 3*b[2]),
             (3, 2*b[1]+b[2]), (3, b[0]+2*b[2]),
             (3, b[1]+2*b[2]), (6, m)]
    return sum(Q(weight, 27)*Q(1, 2)**power for weight, power in terms)


def generating_function(n, d):
    z = 1-d
    result = Q(0)
    for k, l, neutral in compositions(n, 3):
        weight = Q(factorial(n), factorial(k)*factorial(l)*factorial(neutral))
        weight *= Q(1, n)**(k+l)*Q(n-2, n)**neutral
        result += weight*z**max(k, l)
    return result


def main():
    lattice_cases = 0
    strict_transfer_cases = 0
    for m in range(1, 9):
        for aa, bb in product(list(compositions(m, 3)), repeat=2):
            a, b = chamber(aa, bb)
            original = overlap(a, b)
            check(original == overlap(aa, bb) == chamber_expression(a, b))
            a2 = [a[0]+a[2], a[1], 0]
            b2 = [0, b[1], b[2]+b[0]]
            reduced = overlap(a2, b2)
            check(reduced <= original)
            if a[2] or b[0]:
                check(reduced < original)
                strict_transfer_cases += 1
            d = [x-y for x, y in zip(a2, b2)]
            check(d[0] >= d[1] >= 0 >= d[2])
            check(chamber_expression(a2, b2) == reduced)
            lattice_cases += 1

    r, u, v = sp.symbols('r u v', positive=True)
    fprime = -u**3 + r**3/u**3-r*u-2*u**2+r**2/u
    gprime = r**3/v**3-r*v+r**2/v+2*r**2/v**2
    sextic = u**6+2*u**5+r*u**4-r**2*u**2-r**3
    quartic = v**4-r*v**2-2*r*v-r**2
    check(sp.expand(-u**3*fprime-sextic) == 0)
    check(sp.expand(-v**3*gprime/r-quartic) == 0)
    g_at_one_minus_x = u**3-r**2/u+r*u+2*u**2
    check(sp.simplify(fprime+g_at_one_minus_x-r**3/u**3) == 0)
    cubic = r**3+3*r**2-r-1
    quadratic = r**2+3*r-1
    check(sp.rem(cubic, quadratic, r) == -1)

    numerical = []
    for n in range(3, 13):
        for d in [Q(1, 20), Q(1, 2), Q(19, 20), Q(1)]:
            q = sqrt(float(1-d))
            def integrand(t):
                kernel = (1-q*q)/(1-2*q*cos(t)+q*q)
                return kernel*(1-2*(1-q*cos(t))/n)**n/(2*pi)
            integral, error = quad(integrand, -pi, pi, epsabs=1e-12)
            direct = float(generating_function(n, d))
            check(abs(integral-direct) < 1e-11)
            numerical.append({'n': n, 'D': str(d),
                              'absolute_error': abs(integral-direct),
                              'quadrature_error_estimate': error})
    report = {
        'status': 'passed',
        'scope': 'Finite lattice, four exact symbolic identities, and finite quadrature checks; no proof-assistant verification.',
        'lattice': {'sum_ranges': [1, 8], 'base_common_mass': '1/2',
                    'pairs_checked': lattice_cases,
                    'strict_boundary_transfers': strict_transfer_cases},
        'symbolic_identities': ['sextic derivative equation',
                                'quartic derivative equation',
                                'coupled-constraint derivative identity',
                                'threshold cubic remainder -1 modulo quadratic'],
        'Poisson_kernel_checks': numerical,
    }
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
