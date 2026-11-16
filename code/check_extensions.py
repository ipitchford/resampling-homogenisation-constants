#!/usr/bin/env python3
"""Checks for the results added in revision (Propositions 4.3 and 7.3, Corollary 5.2, Remark 7.4).

Finite exact or high-precision checks only; the general statements rest on the proofs in the manuscript.
Every check raises CheckFailure on failure, so the script also fails under ``python -O``.
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import json
import sys

import mpmath as mp
import sympy as sp


class CheckFailure(Exception):
    """Raised when a verification check fails."""


def check(condition, label):
    if not condition:
        raise CheckFailure(label)


BASE = Path(__file__).resolve().parents[1]
mp.mp.dps = 40
report = {}

# 1. C_n = E max(K,L) for (K,L,n-K-L) ~ Mult(n; 1/n,1/n,1-2/n): exact values and the Poisson limit.
def C_exact(n):
    total = F(0)
    for k in range(n + 1):
        for l in range(n - k + 1):
            w = F(factorial(n), factorial(k) * factorial(l) * factorial(n - k - l))
            w *= F(1, n) ** (k + l) * F(n - 2, n) ** (n - k - l)
            total += w * max(k, l)
    return total


def C_float(n):
    # log-space evaluation for large n
    s = mp.mpf(0)
    ln = mp.log
    for k in range(0, min(n, 80) + 1):
        for l in range(0, min(n - k, 80) + 1):
            lw = mp.loggamma(n + 1) - mp.loggamma(k + 1) - mp.loggamma(l + 1) - mp.loggamma(n - k - l + 1)
            lw += (k + l) * ln(mp.mpf(1) / n) + (n - k - l) * ln(mp.mpf(n - 2) / n)
            s += mp.e ** lw * max(k, l)
    return s


C_inf_bessel = 1 + mp.e ** -2 * (mp.besseli(0, 2) + mp.besseli(1, 2))
# independent series: E max(X,Y), X,Y iid Poisson(1)
pois = [mp.e ** -1 / mp.factorial(j) for j in range(120)]
C_inf_series = mp.fsum(pois[a] * pois[b] * max(a, b) for a in range(120) for b in range(120))
check(abs(C_inf_bessel - C_inf_series) < mp.mpf(10) ** -30, "C_inf Bessel formula")
# telescoping identity sum_{k>=1} k I_k(2) = I_0(2) + I_1(2)
tel = mp.nsum(lambda k: k * mp.besseli(int(k), 2), [1, mp.inf])
check(abs(tel - (mp.besseli(0, 2) + mp.besseli(1, 2))) < mp.mpf(10) ** -25, "Bessel telescoping identity")

exact = {n: C_exact(n) for n in range(2, 41)}
check(exact[2] == F(3, 2) and exact[3] == F(14, 9), "C_2, C_3")
check(all(exact[n] > exact[n + 1] for n in range(3, 40)), "strict decrease n>=3 (to 40)")
check(all(mp.mpf(exact[n].numerator) / exact[n].denominator > C_inf_bessel for n in range(3, 41)), "C_n > C_inf")
large = {n: C_float(n) for n in (100, 400, 1600)}
check(all(large[a] > large[b] > C_inf_bessel for a, b in ((100, 400), (400, 1600))), "large-n monotone")
check(large[1600] - C_inf_bessel < mp.mpf("1e-3"), "convergence to C_inf")
report["C_infinity"] = {
    "bessel_formula": mp.nstr(C_inf_bessel, 20),
    "poisson_series": mp.nstr(C_inf_series, 20),
    "C_n_exact_selected": {str(n): str(exact[n]) for n in (2, 3, 4, 5, 10, 20, 40)},
    "C_n_decimal_large_n": {str(n): mp.nstr(v, 12) for n, v in large.items()},
    "one_over_C_infinity": mp.nstr(1 / C_inf_bessel, 12),
}

# 2. Corollary 5.2 (n=2, arbitrary measures): 1 - F_2(D) = r(1+r)/2 with r = 1-D.
D = sp.symbols("D")
F2 = sp.Rational(3, 2) * D - D ** 2 / 2
r = 1 - D
check(sp.expand(1 - F2 - r * (1 + r) / 2) == 0, "n=2 overlap identity")
report["n2_identity"] = "1 - (3D/2 - D^2/2) = r(1+r)/2, r = 1-D"


# 3. Proposition 7.3: base-B separation for every integer B >= k+1.
def compositions(total, dim):
    if dim == 1:
        yield (total,)
        return
    for a in range(total + 1):
        for rest in compositions(total - a, dim - 1):
            yield (a,) + rest


sep_rows = []
for n in range(2, 6):
    for k in range(1, 4):
        vecs = list(compositions(k * n, n))
        for B in range(k + 1, k * n + 2):
            w = [B ** i for i in range(n)]
            target = k * sum(w)
            gaps = []
            for N in vecs:
                delta = sum(c * x for c, x in zip(N, w)) - target
                bal = N == (k,) * n
                check((delta == 0) == bal, f"uniqueness n={n} k={k} B={B} N={N}")
                check(delta % (B - 1) == 0, f"divisibility n={n} k={k} B={B}")
                if not bal:
                    gaps.append(abs(delta))
            check(min(gaps) >= B - 1, f"gap n={n} k={k} B={B}")
            sep_rows.append({"n": n, "k": k, "B": B, "vectors": len(vecs), "min_abs_Delta": min(gaps),
                             "delta": str(F((B - 1) ** 2, n * (B ** n - 1)))})
report["separation_checks"] = len(sep_rows)


def delta_B(n, k, B):
    return F((B - 1) ** 2, n * (B ** n - 1))


table = []
for n, k in [(2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (3, 2), (4, 2), (5, 3)]:
    first = delta_B(n, k, k * n + 1)
    best_B = max(range(k + 1, k * n + 2), key=lambda B: delta_B(n, k, B))
    small = delta_B(n, k, k + 1)
    M = max(sum(1 for N in compositions(s, n) if max(N) <= k) for s in range(k * n + 1))
    pigeon = F(2 * k, M - 1)
    check(delta_B(n, k, best_B) <= pigeon, "pigeonhole upper bound for the symmetric family")
    table.append({"n": n, "k": k, "delta_first_version_B=kn+1": str(first),
                  "delta_B=k+1": str(small), "best_B_in_range": best_B,
                  "delta_best": str(delta_B(n, k, best_B)),
                  "symmetric_family_upper_bound_2k/(M-1)": str(pigeon),
                  "decimal": [float(first), float(small), float(delta_B(n, k, best_B)), float(pigeon)]})
report["delta_table"] = table

# 3b. Remark after Corollary 6.2: 1 - H_3(D) = (2/9) r + c r^{6/5} + o(r^{6/5}), c = 5 * 2^{-2/5} / 27.
def bisect(fun, lo, hi, iters=400):
    flo = fun(lo)
    check(flo * fun(hi) < 0, "bracket")
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = fun(mid)
        if fm == 0:
            return mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def overlap3(rr):
    rr = mp.mpf(rr)
    u = bisect(lambda t: t ** 6 + 2 * t ** 5 + rr * t ** 4 - rr ** 2 * t ** 2 - rr ** 3, rr, mp.sqrt(rr))
    v = bisect(lambda t: t ** 4 - rr * t ** 2 - 2 * rr * t - rr ** 2, mp.sqrt(rr), mp.mpf(1))
    f = u ** 3 + rr ** 3 / u ** 3 + 3 * rr * u + 3 * u ** 2 + 3 * rr ** 2 / u
    g = rr ** 3 / v ** 3 + 3 * rr * v + 3 * rr ** 2 / v + 3 * rr ** 2 / v ** 2
    return (f + g + 6 * rr) / 27


mp.mp.dps = 80
c65 = 5 * mp.mpf(2) ** (-mp.mpf(2) / 5) / 27
ratios = {}
for e in (10, 20, 40):
    rr = mp.mpf(10) ** -e
    ratios[str(e)] = (overlap3(rr) - 2 * rr / 9) / rr ** (mp.mpf(6) / 5)
check(abs(ratios["40"] - c65) < abs(ratios["10"] - c65) and abs(ratios["40"] - c65) < mp.mpf("1e-3"), "next-order term")
report["three_component_next_order"] = {"predicted_coefficient": mp.nstr(c65, 12),
                                        "observed_at_r=1e-10,1e-20,1e-40": [mp.nstr(ratios[k], 12) for k in ("10", "20", "40")]}
mp.mp.dps = 40

# 4. Lower bound valid at every distance (Theorem 7.1, uniform resampling, k=1): ratio of Bayes errors >= n!/n^n.
report["balanced_coefficients_k1"] = {str(n): str(F(factorial(n), n ** n)) for n in range(2, 7)}

report["status"] = "All listed checks passed"
out = BASE / "logs" / "extensions.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({"status": report["status"], "C_infinity": report["C_infinity"]["bessel_formula"],
                  "separation_checks": len(sep_rows), "output": str(out)}, indent=2))
sys.exit(0)
