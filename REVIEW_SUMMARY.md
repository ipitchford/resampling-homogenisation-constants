# Review summary (public, 0.1.0-candidate)

This is an internal, AI-assisted review record, not peer review. Five reviewer roles (editor, methodology, domain,
perspective, devil's advocate) read the first draft with its sources; a second AI model carried out an adversarial
review; the lead replayed the checks from a fresh copy. One revision round followed, with extensions. The revision was
then checked by a re-review, a proof audit of every new or rewritten proof, and an audit of the new proofs by the second
AI model. No reviewer or replay found a mathematical error in any theorem, and every load-bearing number reproduced.

## Objections and how they were handled

| Objection | Disposition |
|---|---|
| The relation to the published bounds was unclear: 9/14 improves the general upper bound c <= 8/9, and this was not said | The abstract, introduction, Theorem 4.1 and Proposition 4.2 now state that c* <= 9/14 (two sources: 2/3; every n >= 3: below 0.6563), that 8/9 is therefore not the general optimum, and that the Bernoulli conjecture is untouched |
| The class restriction (erasure experiments) was not visible in headline statements; the class might be read as a reduction target | Every headline result names its class and regime; Section 8.3 shows why the erasure class is not a reduction for general experiments |
| Credit: the labelled lift is Kontorovich's; the n = 2 linear bound appears in his local-to-global paper; the Bernoulli example of Section 8 is his | Attributed in Sections 1.4, 2 and 8.3 and Remark 5.3, with theorem-level locators; Remark 4.3 notes that any proof through the labelled lift cannot give c > 9/14 |
| Related work missing: privacy amplification by subsampling; the classical no-collision coefficient | Balle-Barthe-Gaboardi (Theorem 10, the one-differing-record counterpart) and Darroch-Jirina-Speed cited and compared in Section 1.4 |
| The erasure class was defined by reference to a cryptography paper that could not be checked | Replaced by an intrinsic definition (likelihood ratio with at most one value in (0, infinity)); that citation removed |
| Context and readability: the draft assumed the reader knew the problem; jargon in the abstract; notation clashes; process prose | New plain introduction (problem, what was known, what is new with class and regime for each, relation to other work, roadmap); each section opens with what it does; notation renamed where it clashed; author line "Anonymous" and an AI-use statement |
| Precision points (non-attainment wording, ordering of revealed labels, a binary-collapse sentence, compactness in Corollary 6.2, uniformity in Theorem 7.1, Bayes-error ratios only for D < 1) | Corrected |
| The checks were assert-based and so silent under `python -O` | Every script now raises an explicit `CheckFailure`; all rerun in normal and `-O` mode; negative controls fail as they should |
| Extension: the n = 2 problem for arbitrary experiments | Done: Corollary 5.2, sup H = 3D/2 - D^2/2 |
| Extension: the large-n behaviour of C_n | Done: Proposition 4.2, C_n decreases to 1 + e^-2 (I_0(2) + I_1(2)), about 1.5238 |
| Extension: a sharper high-signal construction | Done: any radix B >= k + 1 works, with an explicit remainder exponent (Theorem 7.1, Remark 7.2) |
| Extension: evidence on whether 9/14 is the general constant | Done as evidence only: Conjecture 8.1 with a documented fixed-seed numerical search |
| Points from the verification round (a vacuous sentence in Remark 4.3, an overstated sentence about extremal examples, "attained" versus "witnessed", the search tolerance, the erasure description in the abstract, bibliography year order, the figure caption) | Corrected before freezing |

## Not attempted (recorded as open problems in Section 8)

A proof of c* = 9/14, even for n = 3 or as D -> 0; the Bernoulli constant; the exact optimum for four or more erasure
experiments; budgets that are not a multiple of n; hockey-stick (differential privacy) versions; adaptive schedules.
Three further classical references suggested in review (on sampling with and without replacement, deficiency under
resampling, and comparison of experiments) were judged not load-bearing and are not cited.
