# Resampling heterogeneous experiments: sharp amplification constants and high-signal error asymptotics

*In a natural class, resampling raises a test's advantage by at most 14/9; in general the homogenisation constant is at most 9/14*

Version 0.1.0-candidate. DOI [10.5281/zenodo.24267748](https://doi.org/10.5281/zenodo.24267748) (all versions:
[10.5281/zenodo.24267747](https://doi.org/10.5281/zenodo.24267747)). Creator: Anonymous. Release page:
https://evidencepress.org/releases/resampling-homogenisation-constants/

**Status: unrefereed research candidate.** The checking behind it is internal AI-assisted review and computation, not
peer review (see "How it was checked" below).

## What this is about

**The everyday problem.** You want to decide between two explanations using several different sources of evidence,
say n different sensors or n different sites. Design A looks at every source once. Design B makes n draws at random
with replacement, so some sources are seen twice and some not at all, as a bootstrap or a random-design study does. You
would expect B to be worse, since it may miss a source altogether. Surprisingly, it can be slightly better. The question
is by how much.

**How "better" is measured.** The quality of a design is the total-variation distance between what the data look like
under the two explanations. It is the advantage over guessing of the best possible test: with equal prior odds the best
test is right with probability (1 + distance)/2. Write D for the distance under design A and H for the distance under
design B.

**What was known.** Aryeh Kontorovich (arXiv preprints, 2026) proved that resampling can never help by more than a
constant factor: D >= cH for a universal constant c > 0.1489, whatever the sources and however many there are. In the
other direction, his yes/no (Bernoulli) example showed c <= 8/9, and 8/9 was conjectured to be optimal for Bernoulli data.

**What is new.**

1. *A better bound for the general constant* (Theorem 4.1): c <= 9/14, about 0.643. Two sources already give 2/3. So
   8/9 is not the answer for general data. The examples need at least three outcome symbols and say nothing about the
   Bernoulli conjecture.
2. *Exact worst cases in a natural class* (Theorems 3.1 and 4.1, Proposition 4.2). In an "erasure" experiment each
   observation either rules one explanation out or lands in a shared zone where both explanations behave alike, though
   they may reach it with different frequencies. For such experiments the exact worst-case factor sup H/D is found for
   every random sampling rule chosen independently of the data, and it depends only on how often pairs of sources are
   drawn. For uniform resampling the factor is 3/2 for two sources, 14/9 for three (its maximum), and decreases towards
   about 1.524 for many sources.
3. *Two sources, all data* (Theorem 5.1, Corollary 5.2). With two sources the problem is solved completely for
   arbitrary experiments: the largest possible H is 3D/2 - D^2/2. More generally the exact bound is given at every
   distance when only two sources differ, and the exact optimum for three erasure sources is found (Theorem 6.1).
4. *Strong evidence* (Theorem 7.1, Corollary 7.3). For arbitrary experiments, resampling with n draws always leaves at
   least n!/n^n of the fixed design's error probability, the chance that every source is drawn exactly once (2/9 for
   three sources), and this factor is attained when the evidence is already strong. A version with kn draws is proved
   too.

**A concrete example.** Three sources: two carry weak opposing "rules out" signals and one carries nothing. Looking
once at each gives advantage D. Three random draws give exactly 14D/9 - 2D^2/3 + 2D^3/27, about 14/9 times D when D is
small. That is 14/9 times the advantage over guessing, not 55.6% more accuracy.

**The idea of the proof.** In the erasure class, the information carried by any sampling outcome is a maximum of two
linear budgets, one for each explanation. The worst case puts each explanation's whole budget on a single source, so
only pairs of counts matter, which gives the exact factor. For strong evidence, a base-(k+1) "digit" construction makes
every unbalanced draw much less informative than the balanced one.

**Limits.** Whether 9/14 is the exact general constant is open: it is conjectured, and a numerical search found nothing
larger, but that is not a proof. The erasure-class results are not known to bound general data. The high-signal limit
can be approached slowly. The Bernoulli conjecture 8/9 is not addressed.

**Credit.** The homogenisation problem, the constant c > 0.1489, the example giving 8/9 and the labelled-lift argument
used here are Kontorovich's. The coefficient n!/n^n is the classical probability that n draws with replacement are all
distinct (Darroch, Jirina and Speed); privacy amplification by subsampling (Balle, Barthe and Gaboardi) studies the
one-differing-record counterpart. The novelty is in the sharp constants and the high-signal attainment theorem.

## What is in this package

| File | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper (16 pages). Section 1 is a plain introduction with a roadmap; Section 8 discusses interpretation, the conjecture and open problems. |
| [code/verify.py](code/verify.py) | Exact rational checks of the extremal families on finite product spaces, the occupancy constants C_n, symbolic identities, and the three-component formula against constrained optimisation. Output [logs/verification.json](logs/verification.json). |
| [code/check_schedule.py](code/check_schedule.py), [code/check_two_active.py](code/check_two_active.py), [code/check_n3_reduction.py](code/check_n3_reduction.py) | Checks of the schedule constant, the two-differing-coordinate bound and the three-component reduction (4916 rational chamber cases). Outputs `logs/*.stdout.json`. |
| [code/check_boundary.py](code/check_boundary.py) | The Bernoulli example of Section 8.3 (87/400, 111/400). Output [logs/boundary.json](logs/boundary.json). |
| [code/verify_unrestricted_asymptotic.py](code/verify_unrestricted_asymptotic.py) | The high-signal radix construction for nine (n, k, B) cases, 107 count vectors. Output [logs/unrestricted_asymptotic.json](logs/unrestricted_asymptotic.json). |
| [code/check_extensions.py](code/check_extensions.py) | The large-n limit, base-B separation (72 cases), the delta table and the next-order coefficient. Output [logs/extensions.json](logs/extensions.json). |
| [code/search_unrestricted.py](code/search_unrestricted.py) | Fixed-seed numerical search behind the conjecture (about 2.5 minutes). Output [logs/search_unrestricted.json](logs/search_unrestricted.json). |
| [code/make_figures.py](code/make_figures.py), `figures/sharp_homogenisation.{pdf,png,svg}` | Figure 1 and the script that draws it. |
| [RESEARCH_GATES.json](RESEARCH_GATES.json), [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | Contribution map, sources read and limitation triage; objections from internal review and how they were handled. |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits. |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [requirements.txt](requirements.txt) | Replay record and environment. |
| [refs.bib](refs.bib) | Bibliography built from Crossref and arXiv records. |
| [LICENSE.md](LICENSE.md), [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json), [MANIFEST.sha256](MANIFEST.sha256) | Licences, citation metadata and SHA-256 sums. |

## Reproduce the checks

From the package root, with [uv](https://docs.astral.sh/uv/) (the scripts write into `logs/` and `figures/`, so work on
a copy):

```
UV="uv run --no-project --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with sympy==1.14.0 --with mpmath==1.3.0 --with matplotlib==3.10.8"
for s in verify check_boundary check_schedule check_two_active check_n3_reduction verify_unrestricted_asymptotic check_extensions; do
  $UV python -I -B -O code/$s.py || echo "FAILED: $s"
done
$UV python -I -B code/make_figures.py
$UV python -I -B code/search_unrestricted.py     # about 2.5 minutes
```

Every script exits 1 if any check fails and uses no `assert`, so it behaves identically under `python -O`. The
recorded logs are reproduced exactly apart from wall-clock timings (and possibly the last digit of one optimiser float
in `verification.json` on other platforms). The exact checks test the stated formulas and constructions; the numerical
search supports the conjecture without proving it. The general statements rest on the written proofs.

## How it was checked

Internal, AI-assisted review only: five reviewer roles, an adversarial review and an audit of the new proofs by a second
AI model, clean replays of the scripts, one revision round, a re-review and a proof audit. This is not peer review,
independent reproduction or formal verification. The paper's closing statement describes the use of AI.

## Licence

Text, paper, figure and data: CC0-1.0. Code (`code/*.py`, `.github/workflows/*.yml`): MIT. See
[LICENSE.md](LICENSE.md).
