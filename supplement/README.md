# Verification supplement

This directory contains the reproducible code and small data tables used to
check the finite computations for inversion-sequence Classes 152 and 207.
The written mathematical proofs are outside this code-only repository. The
Class 207 all-length equality follows from the summation argument; the
scripts below check its finite coefficients and algebraic identities.

## Files and scope

- `counts.csv`: lengths 0 through 12. The `B207` and `C207` columns are
  generated separately by the two tree implementations in `check_counts.py`.
  `A152` and `B152` use the central-binomial formula. Literal pattern
  enumeration independently checks all four columns through length 8.
- `class207_recurrence.py`: standard-library implementation of the Class 207
  C-site recurrence.
- `check_counts.py`: independent tree, literal pattern and fixed-maximum checks.
- `symbolic_check.py`: exact rational identities checked with SymPy.
- `opposite_catalysts.py`: exact signed-specialization check through length 60.
- `signed_check.py`: exact Gaussian-integer six-term specialization through
  length 60.
- `class207-proof/`: exact Lambert/product checks through order 128,
  symbolic certificates, and independent referee programs.
- `requirements.txt`: pinned symbolic-check dependencies.
- `expected-*.json`: concise expected outputs for the corresponding commands.

## Reproduce the checks

Run the standard-library checks from this directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 check_counts.py
PYTHONDONTWRITEBYTECODE=1 python3 class207_recurrence.py 12 --check counts.csv
PYTHONDONTWRITEBYTECODE=1 python3 opposite_catalysts.py
PYTHONDONTWRITEBYTECODE=1 python3 signed_check.py
```

Expected output is recorded in the `expected-*.json` files. The first command
prints a JSON object whose status is `PASS`; the recurrence command prints:

```text
Verified all 13 terms (0 <= n <= 12).
```

For the exact symbolic checks, use an environment outside the repository:

```sh
python3 -m venv /tmp/class152-207-symbolic
/tmp/class152-207-symbolic/bin/python -m pip install -r requirements.txt
PYTHONDONTWRITEBYTECODE=1 /tmp/class152-207-symbolic/bin/python symbolic_check.py
```

The Class 207 proof checks are run from `class207-proof/`:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 r3_exact.py --order 128 --counts 128
PYTHONDONTWRITEBYTECODE=1 python3 r3_symbolic.py
PYTHONDONTWRITEBYTECODE=1 python3 referee/c_kernel.py
PYTHONDONTWRITEBYTECODE=1 python3 referee/audit.py
PYTHONDONTWRITEBYTECODE=1 python3 referee/published_tree.py
PYTHONDONTWRITEBYTECODE=1 python3 referee/literal_formula.py
```

`r3_exact.py` writes the reproducible `r3-counts.csv`, which is ignored by
`.gitignore`; its expected checksum and the JSON reports are supplied in this
directory. No Lean build is claimed for Class 207. The finite checks use exact
integer or rational arithmetic, except for the explicitly supplementary
70-decimal audit in `referee/audit.py`.

## Citation

Please cite the associated article as:

> Toufik Mansour, Wenlin Zhang, Haobo Ma, *Counting formulas and refined
generating functions for inversion-sequence Classes 152 and 207*.

## Licenses

Python source code is MIT-licensed. The CSV and JSON data and expected-output
records are available under CC BY 4.0; see the repository `LICENSE` file.
