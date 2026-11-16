#!/usr/bin/env python3
"""Numerical search for large homogenisation ratios H/D over ARBITRARY finite experiments (Section 8.2).

For n pairs (P_i, Q_i) on an m-point alphabet, the script maximises
  * 'labelled': H_lab / D, where H_lab is the TV of the uniformly resampled experiment with component labels kept
    (H_lab >= H, so this is the harder target for the conjecture c* = 9/14), and
  * 'shared':   H / D with all components on one common m-point alphabet and labels forgotten
    (m = 2 is Kontorovich's Bernoulli setting, conjectured constant 8/9, i.e. ratio 9/8),
by Nelder-Mead from random softmax-parametrised starts with fixed seeds.

This is a bounded random-restart search: it can find counterexamples but proves nothing about the supremum.
Every check raises explicitly, so the script also fails under ``python -O``.
"""
from itertools import product
from math import factorial
from pathlib import Path
import json
import time

import numpy as np
from scipy.optimize import minimize

BASE = Path(__file__).resolve().parents[1]


class CheckFailure(Exception):
    pass


def check(condition, label):
    if not condition:
        raise CheckFailure(label)


def prod(dists):
    out = np.array([1.0])
    for d in dists:
        out = np.multiply.outer(out, d).ravel()
    return out


def tv(p, q):
    return 0.5 * float(np.abs(p - q).sum())


def counts(n):
    return [N for N in product(range(n + 1), repeat=n) if sum(N) == n]


def weight(N):
    n = len(N)
    w = factorial(n) / n ** n
    for c in N:
        w /= factorial(c)
    return w


def C_n(n):
    out = 0.0
    for k in range(n + 1):
        for l in range(n - k + 1):
            out += factorial(n) / (factorial(k) * factorial(l) * factorial(n - k - l)) * n ** -(k + l) * ((n - 2) / n) ** (n - k - l) * max(k, l)
    return out


def make_ratio(n, m, mode):
    Ns = counts(n)
    W = [weight(N) for N in Ns]

    def unpack(z):
        z = z.reshape(2 * n, m)
        e = np.exp(z - z.max(1, keepdims=True))
        e /= e.sum(1, keepdims=True)
        return e[:n], e[n:]

    def ratio(z):
        P, Q = unpack(z)
        D = tv(prod(P), prod(Q))
        if D < 1e-10:
            return 0.0
        if mode == "labelled":
            H = sum(w * tv(prod([P[i] for i in range(n) for _ in range(N[i])]),
                           prod([Q[i] for i in range(n) for _ in range(N[i])])) for w, N in zip(W, Ns))
        else:
            H = tv(prod([P.mean(0)] * n), prod([Q.mean(0)] * n))
        return H / D

    return ratio, unpack


CONFIGS = [("labelled", 2, 2), ("labelled", 2, 3), ("labelled", 3, 2), ("labelled", 3, 3), ("labelled", 3, 4),
           ("labelled", 4, 2), ("shared", 2, 2), ("shared", 3, 2), ("shared", 3, 3)]
RESTARTS = 40
MAXITER = 4000

rows = []
start = time.time()
for idx, (mode, n, m) in enumerate(CONFIGS):
    rng = np.random.default_rng(20261115 + idx)
    ratio, unpack = make_ratio(n, m, mode)
    best, bz = 0.0, None
    for _ in range(RESTARTS):
        z0 = rng.normal(size=2 * n * m) * 3
        res = minimize(lambda z: -ratio(z), z0, method="Nelder-Mead",
                       options={"maxiter": MAXITER, "xatol": 1e-9, "fatol": 1e-12})
        if -res.fun > best:
            best, bz = -res.fun, res.x
    P, Q = unpack(bz)
    D = tv(prod(P), prod(Q))
    ref = C_n(n) if mode == "labelled" else None
    if mode == "labelled":
        check(best <= ref + 1e-6, f"search exceeded C_n: {mode} n={n} m={m} best={best}")
    rows.append({"mode": mode, "n": n, "alphabet": m, "restarts": RESTARTS, "seed": 20261115 + idx,
                 "best_ratio_found": round(best, 8), "D_at_best": round(D, 8),
                 "C_n_reference": None if ref is None else round(ref, 8),
                 "P_at_best": np.round(P, 5).tolist(), "Q_at_best": np.round(Q, 5).tolist()})
    print(mode, n, m, round(best, 6), "C_n" if ref else "", None if ref is None else round(ref, 6), flush=True)

out = {"scope": "Bounded random-restart search over arbitrary finite experiments; evidence, not proof.",
       "configs": rows, "elapsed_seconds": round(time.time() - start, 1),
       "status": "No labelled ratio above C_n was found"}
(BASE / "logs").mkdir(exist_ok=True)
(BASE / "logs" / "search_unrestricted.json").write_text(json.dumps(out, indent=2) + "\n")
print(out["status"])
