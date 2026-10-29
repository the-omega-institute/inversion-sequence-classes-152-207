# Independent referee computations for Class 207

These programs were written separately from the production verification code.
They do not import the production generator or one another.
Run each program from this directory; it writes the correspondingly named JSON.

- `python3 literal_formula.py`: standard library only. Constructs the small
  root, Lambert series, and total formula with eight guard coefficients;
  compares the length coefficients through 8 with literal inversion-sequence
  enumeration for all three pattern sets A, B, and C.
- `python3 published_tree.py`: standard library only. Implements the three-type
  tree of Callan--Mansour, Table 1, p. 28. Checks its totals and the separate
  boundary `f(x)` scalar equation through degree 14 at x = 2, 3, -1.
- `python3 c_kernel.py`: requires SymPy. Starts from the manuscript's cleared
  C equation and common P equation, checks the diagonal identity modulo the
  kernel, and checks its transformation to the two-root input.
- `python3 audit.py`: requires SymPy and mpmath. Checks ten exact certificates,
  including both exceptional-term algebra and the all-index tail ratio.
  Infinite products are canceled with finite shift relations. It also checks
  the two-parameter summation numerically at three parameter pairs which do
  not satisfy the quadratic relation. The numerical computations use 70
  decimal working precision, 190 orbit terms, and 219 Lambert terms; they are
  supplementary error checks rather than a proof.

The exact checks require no network access. SymPy 1.14.0 and mpmath 1.3.0
were used in the referee environment. The mathematical proof is in the note.
The primary references verified during review are DLMF 17.7.5 (including
its condition ef = abcq), and Callan--Mansour, Integers 23 (2023), Paper A78,
Table 1, p. 28 and section 4.10. Primary-source PDFs are not part of this
source package.
