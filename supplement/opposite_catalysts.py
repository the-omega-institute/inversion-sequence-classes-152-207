#!/usr/bin/env python3
"""Exact C-system coefficients and the opposite-catalyst test (20 October 2026)."""
from collections import Counter
import json
from pathlib import Path

def coefficients(N):
    F, P, Z = [Counter()], [Counter()], [Counter()]
    for n in range(1, N+1):
        pn, zn = Counter(), Counter()
        if n == 1:
            pn[1] = 1
        for i, c in P[-1].items():
            pn[i+1] += c
            for j in range(i):
                zn[i-j, j] += c
        for i, c in F[-1].items():
            for j in range(i):
                pn[j+1] += c
            for a in range(i-1):
                for b in range(i-1-a):
                    zn[a+1, b] += c
        for (i, j), c in Z[-1].items():
            for b in range(j):
                zn[i+1, b] += c
            for a in range(i):
                zn[i-a, j+a] += c
        fn = pn.copy()
        for (i, j), c in Z[-1].items():
            fn[i+1] += c
        P.append(pn); Z.append(zn); F.append(fn)
    return F, P, Z

def anti(poly):
    ans=Counter()
    for (i,j), c in poly.items():
        ans[i+j] += (-1)**j*c
    return {i:c for i,c in sorted(ans.items()) if c}

def main():
    F,P,Z=coefficients(60)
    from class207_recurrence import counts
    assert [sum(f.values()) for f in F][1:] == counts(60)[1:]
    differences=[]
    for n, poly in enumerate(Z):
        a=anti(poly)
        d={i:2*c for i,c in a.items() if i%2==0}
        if d: differences.append({'n':n,'coefficient_of_sum':d})
    report={'order':60,'definition':'Z_C=H_C-P_C, Proposition 3',
            'identity':not differences,'first_failure':differences[0] if differences else None,
            'first_specializations':{n:anti(Z[n]) for n in range(2,9)},
            'nonzero_discrepancy_orders':[d['n'] for d in differences]}
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
