#!/usr/bin/env python3
"""C legal-site and B four-type recurrences, copied without changes from
our earlier generate_push2.py. These are separate from the formula builder;
the independent literal enumeration is in referee/literal_formula.py.
"""
from collections import Counter

def c_coefficients(nmax: int):
    F = [Counter()]
    P = [Counter()]
    Z = [Counter()]
    for n in range(1, nmax + 1):
        pn, zn = Counter(), Counter()
        if n == 1:
            pn[1] = 1
        for i, c in P[-1].items():
            pn[i + 1] += c
            for j in range(i):
                zn[i - j, j] += c
        for i, c in F[-1].items():
            for j in range(i):
                pn[j + 1] += c
            for a in range(i - 1):
                for b in range(i - 1 - a):
                    zn[a + 1, b] += c
        for (i, j), c in Z[-1].items():
            for b in range(j):
                zn[i + 1, b] += c
            for a in range(i):
                zn[i - a, j + a] += c
        fn = pn.copy()
        for (i, j), c in Z[-1].items():
            fn[i + 1] += c
        F.append(fn)
        P.append(pn)
        Z.append(zn)
    return F, P, Z


def b_totals(nmax: int):
    states = Counter({("R", 1, 0): 1})
    out = [1]
    for _n in range(1, nmax + 1):
        out.append(sum(states.values()))
        children = Counter()
        for (kind, slack, gap), multiplicity in states.items():
            if kind in ("R", "U", "W"):
                children[("R", slack + 1, 0)] += multiplicity
            if kind != "R":
                stop = gap if kind == "W" else gap + 1
                for k in range(1, stop):
                    children[("V", slack + 1, k)] += multiplicity
            for j in range(1, slack + 1):
                children[("W" if kind == "V" else "U", slack - j + 1, j)] += multiplicity
        states = children
    return out
