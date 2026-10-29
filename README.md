# Inversion-sequence Classes 152 and 207: verification data

This repository is the code-and-data companion for the counting results on
inversion-sequence Classes 152 and 207. It contains the finite enumerators,
recurrence implementations, exact symbolic certificates, and small expected
output tables. It deliberately contains no manuscript source or PDF.

## What is verified

- The two Class 152 avoidance classes agree with the displayed
  central-binomial partial-sum formula through the literal finite checks in
  `supplement/check_counts.py`.
- The independent Class 207 tree implementations agree through length 12,
  and the literal pattern enumerator checks the four recorded columns through
  length 8.
- The Class 207 C-site recurrence and product-orbit formula agree through
  order 128 (`supplement/class207-proof/r3_exact.py`).
- Exact symbolic identities, the kernel certificates, the published-tree
  boundary check, and the signed-specialization checks reproduce the JSON
  records shipped with the data.

The all-length Class 207 equality is a mathematical statement proved in the
associated article; these scripts provide finite and algebraic reproducibility
checks and do not replace that proof.

## Quick start

The standard-library checks need Python 3.10 or newer:

```sh
cd supplement
PYTHONDONTWRITEBYTECODE=1 python3 check_counts.py
PYTHONDONTWRITEBYTECODE=1 python3 class207_recurrence.py 12 --check counts.csv
PYTHONDONTWRITEBYTECODE=1 python3 opposite_catalysts.py
PYTHONDONTWRITEBYTECODE=1 python3 signed_check.py
```

The first command returns a JSON report with `"status": "PASS"`; the second
prints exactly:

```text
Verified all 13 terms (0 <= n <= 12).
```

For the full command list, dependencies, and exact expected JSON, read
[`supplement/README.md`](supplement/README.md). The symbolic checks use
SymPy and mpmath from `supplement/requirements.txt`; the environment should
be created outside this checkout. The order-128 exact check and its referee
checks are in `supplement/class207-proof/`.

## Expected outputs

The committed `expected-count-checks.json`, `expected-symbolic-checks.json`,
`expected-opposite-catalysts.json`, `expected-signed-check.json`, and
`class207-proof/*-check.json` files are the expected machine-readable outputs.
A clean-checkout run is recorded in [`checks/clean-checks.txt`](checks/clean-checks.txt)
with command lines, exit codes, and stdout.

## Citation

Toufik Mansour, Wenlin Zhang, Haobo Ma, *Counting formulas and refined
generating functions for inversion-sequence Classes 152 and 207*.

## License

The Python code is released under the MIT License. The CSV tables and JSON
expected-output data are released under CC BY 4.0. See [`LICENSE`](LICENSE).
