# Class 207: source and exact verification

Section 7 of the accompanying v2 manuscript proves equality of the three Class 207 totals by a cubic limit of
the nonterminating q-Saalschutz identity. The scripts provide exact algebraic
certificates and finite coefficient checks; the all-length conclusion uses
the proof, not extrapolation from these computations.

Run the following commands from this directory.

The main integer computation needs Python 3 and its standard library:

    python3 r3_exact.py --order 128 --counts 128

It expands each product orbit directly, without using the summation identity
to construct the orbit difference. Internal order is 136. It checks the
cleared residual through degree 129 and the total formula against the C
legal-site recurrence and the B four-type tree through length 128. The two
recurrences in r3_counts.py are copied unchanged from the earlier generator;
this reuse is not claimed as a new independent implementation.

This command regenerates r3-counts.csv. The full table is omitted
from the archive. Its expected SHA256 is:

    4ec9b7d8bb7ff14a9e87c9177562432e88aa8577727a4f3aeb9eb4fbd7f67cb8

Install SymPy and mpmath in an environment of your choice for the symbolic
and supplementary numerical checks:

    python3 -m venv .venv
    .venv/bin/python -m pip install sympy mpmath
    .venv/bin/python r3_symbolic.py
    .venv/bin/python referee/c_kernel.py
    .venv/bin/python referee/audit.py

The main symbolic program checks 19 exact rational identities and the cubic
logarithmic multipliers. Identities requiring the branch polynomial are
explicitly marked. Infinite-product constants are reduced by finite shifts.
The independent referee programs were written separately and import no
production checker. The two further checks need only the standard library:

    python3 referee/published_tree.py
    python3 referee/literal_formula.py

The first implements the published three-type B tree and checks its boundary
equation at x=2,3,-1 through length 14. The second constructs the formula
independently, inverts the Catalan substitution, and literally enumerates
all three avoidance classes through length 8. Concise JSON reports accompany
the programs. The 70-decimal off-branch numerical checks in audit.py are
supplementary; the theorem is justified by exact algebra and the written
convergence argument. No Lean verification is claimed for this result.

The B boundary f in the manuscript is the published-tree boundary, not the slack
series in the previous manuscript. Its value at x=1 is G_B-1. This distinction
is essential to the binding of the scalar equation to the counting problem.

Primary sources:

- Callan and Mansour, Integers 23 (2023), A78, Table 1 p. 28 and Section 4.10:
  https://math.colgate.edu/~integers/x78/x78.pdf
  https://doi.org/10.5281/zenodo.8399694
- NIST DLMF, Section 17.7(ii), equation 17.7.5:
  https://dlmf.nist.gov/17.7.E5

This archive contains source and concise reports. It omits third-party PDFs,
regenerable full tables, internal task records, and correspondence.
