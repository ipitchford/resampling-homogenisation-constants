# AI index: Resampling heterogeneous experiments (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What this
is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is a plain
introduction with what was known, what is new and a roadmap). DOI 10.5281/zenodo.24267748 (concept
10.5281/zenodo.24267747). Status: unrefereed candidate; internal AI-assisted review only.

## Normalisation

n independent experiments; under hypothesis 1 experiment i has law P_i, under hypothesis 2 law Q_i, all on one measurable
space. TV(P,Q) = sup_A |P(A) - Q(A)|; overlap 1 - TV = int min(p,q), twice the equal-prior Bayes error. Fixed design:
D = TV(prod_i P_i, prod_i Q_i). Uniform resampling (homogenisation): Pbar = (1/n) sum_i P_i, Qbar likewise,
H = TV(Pbar^n, Qbar^n). c* = best constant in D >= c H over all n and all measurable spaces. Schedule: random counts
N = (N_1, ..., N_n) independent of the data, revealed (labelled) or not; H_nu = TV of the labelled experiment, H <= H_nu
(Kontorovich's labelled lift). Erasure experiment: likelihood ratio takes at most one value in (0, infinity). K, L: the
numbers of draws of two fixed experiments under Mult(n; 1/n, ..., 1/n); C_n = E max(K, L). KJan, KApr, KOct =
Kontorovich arXiv:2601.04079v3, 2604.03882v1, 2610.31164v1.

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A | c* <= 9/14 (about 0.643); n = 2 gives c* <= 2/3; every n >= 3 gives c* <= 1/C_n < 0.6563 | Theorem 4.1; five-outcome family of Section 4.1; Proposition 4.2 | written proof; improves the published bound c <= 8/9 (KApr p. 1); does not bear on the Bernoulli conjecture 8/9 (KJan) |
| B | Erasure experiments, any revealed schedule independent of the data: sup H_nu/D = max{P(N != 0), max_{i,j} E max(N_i, N_j)} | Theorem 3.1 | written proof; new in bounded search; BBG Thm 10 is the one-coordinate counterpart (hockey-stick at eps = 0) |
| C | Uniform resampling of erasure experiments: H <= C_n D, sharp as D -> 0; C_2 = 3/2, C_3 = 14/9, C_4 = 99/64, strictly decreasing for n >= 3 to C_inf = 1 + e^-2 (I_0(2) + I_1(2)) = 1.5237776118... | Theorem 4.1; Proposition 4.2 | written proof; exact values in logs/verification.json and logs/extensions.json |
| D | Arbitrary experiments differing in two coordinates: exact bound H <= F_n(D), attained at every D; n = 2 solved: sup H = 3D/2 - D^2/2; F_3(d) = 14d/9 - 2d^2/3 + 2d^3/27 | Theorem 5.1; Corollary 5.2 | written proof; the linear bound H <= 3D/2 for n = 2 is KOct Section 7.3 |
| E | Exact optimum for three erasure experiments at every distance, transitions D_0 = 0.3248691294..., D_1 = 0.6972243622...; crossover D_* = (9 - sqrt 33)/8; high-signal limit 1 - H_3 = (2/9) r + O(r^{6/5}) | Theorem 6.1; Corollary 6.2 | written proof; checked against constrained optimisation and 4916 rational chamber cases |
| F | Arbitrary experiments, any data-independent schedule with kn draws: 1 - H >= q (1 - D)^k at every distance, q = P(every experiment used exactly k times); leading term attained as D -> 1 with power-saving remainder; uniform k = 1: optimal Bayes-error ratio -> n!/n^n | Theorem 7.1; Remark 7.2; Corollary 7.3 | written proof; n!/n^n classical (Darroch-Jirina-Speed); its simultaneous attainment is new in bounded search |
| G | Conjecture c* = 9/14 | Conjecture 8.1 | conjecture only; fixed-seed search found no labelled ratio above C_n (tolerance 1e-6); shared binary n = 2 best 1.124951, shared ternary n = 3 best 1.436387 |

Not claimed: the exact general constant c*; any bound for general experiments from the erasure constants (Section 8.3
shows two natural transfer routes fail); the Bernoulli constant; finite-distance optima for four or more erasure
experiments or three or more arbitrary experiments; budgets not a multiple of n; adaptive schedules; hockey-stick
versions. The factor 14/9 multiplies total variation (advantage over chance), not accuracy.

## Evidence files

- [code/verify.py](code/verify.py) -> [logs/verification.json](logs/verification.json): 20 exact direct-product
  checks, occupancy constants for 29 sizes, symbolic identities (F_2, F_3, D_*), 9 comparisons of the three-component
  formula with constrained optimisation, 1800 random pairs.
- [code/check_schedule.py](code/check_schedule.py), [code/check_two_active.py](code/check_two_active.py),
  [code/check_n3_reduction.py](code/check_n3_reduction.py) -> [logs/check_schedule.stdout.json](logs/check_schedule.stdout.json),
  [logs/check_two_active.stdout.json](logs/check_two_active.stdout.json),
  [logs/check_n3_reduction.stdout.json](logs/check_n3_reduction.stdout.json).
- [code/check_boundary.py](code/check_boundary.py) -> [logs/boundary.json](logs/boundary.json): Bernoulli example
  D = 1/5, H = 87/400; refinement 111/400.
- [code/verify_unrestricted_asymptotic.py](code/verify_unrestricted_asymptotic.py) ->
  [logs/unrestricted_asymptotic.json](logs/unrestricted_asymptotic.json): radix family, 9 (n, k, B) cases, 107 count
  vectors, one q = 0 schedule.
- [code/check_extensions.py](code/check_extensions.py) -> [logs/extensions.json](logs/extensions.json): C_inf, 72
  separation checks, delta table, next-order coefficient 0.140344126529 against 0.140345102478 observed at r = 1e-40.
- [code/search_unrestricted.py](code/search_unrestricted.py) -> [logs/search_unrestricted.json](logs/search_unrestricted.json):
  40 Nelder-Mead restarts per configuration, one fixed seed per configuration.
- [code/make_figures.py](code/make_figures.py) -> [figures/sharp_homogenisation.pdf](figures/sharp_homogenisation.pdf).
- [requirements.txt](requirements.txt), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md):
  environment and the release replay (normal and `python -O`, negative controls).
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): internal review objections and dispositions.
  [RESEARCH_GATES.json](RESEARCH_GATES.json): contribution map, sources read in full, searches and limitation triage.
- [refs.bib](refs.bib): bibliography from Crossref and arXiv records. [LICENSE.md](LICENSE.md): CC0-1.0 text and data,
  MIT code. [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json): citation metadata. CI:
  [.github/workflows/replay.yml](.github/workflows/replay.yml) (blocking: manifest, exact checks under -O compared with
  the shipped logs, negative controls); [.github/workflows/slow.yml](.github/workflows/slow.yml) (non-blocking: figure
  and search). [MANIFEST.sha256](MANIFEST.sha256): SHA-256 of every file.

## Limits

No human peer review, independent human reproduction or formal verification. The exact checks test the formulas and
finite constructions; the general statements rest on the written proofs. The search is finite and numerical. Four of the
six cited works are arXiv preprints. Priority rests on a bounded search.
