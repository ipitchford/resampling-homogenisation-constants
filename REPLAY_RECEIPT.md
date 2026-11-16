# Replay receipt (release 0.1.0-candidate, 16 November 2026)

The scripts were rerun by the publisher from fresh copies of the approved candidate before the release was frozen, with
the package copy as working directory, CPython 3.12.11 and the pinned requirements (see ENVIRONMENT.txt).

| Run | Mode | Result |
|---|---|---|
| `verify.py` | normal and `-O -I -B` | exit 0, "All listed checks passed": 20 direct products, 29 occupancy sizes, 9 BAEC3 comparisons, 1800 random pairs; log identical to `logs/verification.json` except the wall-clock `elapsed_seconds` |
| `check_boundary.py`, `verify_unrestricted_asymptotic.py`, `check_extensions.py` | normal and `-O -I -B` | exit 0; `logs/boundary.json`, `logs/unrestricted_asymptotic.json` (9 cases, 107 count vectors, 1 q = 0 schedule) and `logs/extensions.json` (C_inf 1.5237776118026086987, 72 separation checks, next-order coefficient 0.140344126529 against 0.140345102478 at r = 1e-40) byte-identical |
| `check_schedule.py`, `check_two_active.py`, `check_n3_reduction.py` | normal and `-O -I -B` | exit 0; standard output identical to `logs/*.stdout.json` |
| `search_unrestricted.py` | `-I -B` | exit 0, about 2.4 minutes; labelled best 1.5, 1.5, 1.555556, 1.555556, 1.555547, 1.546875 (never above C_n); shared 1.124951, 1.072201, 1.436387; log identical to `logs/search_unrestricted.json` except `elapsed_seconds` |
| `make_figures.py` | `-I -B` | exit 0; PNG byte-identical; PDF and SVG identical apart from the creation timestamp and random marker ids |
| negative control: `check_boundary.py` with 87/400 changed to 88/400 | `-O -I -B` | exit 1 (CheckFailure) |
| negative control: `verify_unrestricted_asymptotic.py` with the sandwich lower bound doubled | `-O -I -B` | exit 1 (CheckFailure) |

No script uses `assert`: every check raises `CheckFailure`, so the checks stay active under `-O`.

After this replay only the title-page version/DOI line of the manuscript (PDF rebuilt) and the release records (README,
metadata) changed; no script, figure or recorded log changed. GitHub CI reruns the exact checks under `-O`, compares the
logs with the shipped files and runs both negative controls (blocking); the figure and the search run in a separate
non-blocking workflow.
