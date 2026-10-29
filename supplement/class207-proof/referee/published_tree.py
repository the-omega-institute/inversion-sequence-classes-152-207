from collections import Counter
from math import comb
import json
from pathlib import Path
N=14
states=Counter({('a',1,0):1})
A=[Counter() for _ in range(N+2)]
T=[1]
for n in range(1,N+2):
    T.append(sum(states.values()))
    for (kind,m,j),c in states.items():
        if kind=='a': A[n][m]+=c
    children=Counter()
    for (kind,m,j),c in states.items():
        if kind=='a':
            children['a',m+1,0]+=c
            for k in range(1,m+1): children['b',m,k]+=c
        else:
            mp=m+1-j
            children['a',mp+1,0]+=c
            for k in range(1,mp+1):children['b',mp,k]+=c
            for l in range(j): children['c',m+(kind=='c')-l,j-l]+=c
    states=children
f=[Counter() for _ in range(N+1)]
for n in range(N+1):
    for m,c in A[n+1].items():f[n][m-1]+=c
f[0][0]-=1
for n in range(1,N+1):assert sum(f[n].values())==T[n]
res=[]
for x in (2,3,-1):
    fv=[sum(c*x**m for m,c in row.items()) for row in f]
    phi=[0]*(N+1)
    for n,row in enumerate(f):
        for m,c in row.items():
            for k in range(N+1-n):
                if m:
                    phi[n+k]+=c*comb(m+k-1,k)*x**k
                elif k==0:phi[n]+=c
    lhs=[0]+phi[:-1]
    # RHS (1-zx)(1-x+zx+zx²) f + z(1-x+zx²)(f(1)-x)
    p=[1-x,2*x*x,-x*x*(1+x)]
    d=[1-x,x*x]
    rhs=[sum(p[k]*fv[n-k] for k in range(min(2,n)+1)) for n in range(N+1)]
    g=[sum(row.values()) for row in f];g[0]-=x
    for n in range(1,N+1):rhs[n]+=sum(d[k]*g[n-1-k] for k in range(min(1,n-1)+1))
    fail=next((n for n in range(N+1) if lhs[n]!=rhs[n]),None)
    res.append({'x':x,'first_failure':fail,'residual_prefix':[lhs[n]-rhs[n] for n in range(8)]})
print(json.dumps({'N':N,'totals':T[:13],'residuals':res},indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps({'N':N,'totals':T[:13],'residuals':res},indent=2)+'\n')
