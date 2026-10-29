"""Independent exact formula-to-definition check; standard library only."""
from math import comb
from pathlib import Path
import json
N=8
M=N+8
def mul(a,b):
    return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(M+1)]
def inv(a):
    assert a[0] in (-1,1)
    o=[a[0]]+[0]*M
    for n in range(1,M+1):o[n]=-sum(a[j]*o[n-j] for j in range(1,n+1))//a[0]
    return o
r=[0]*(M+1)
for n in range(1,M+1):
    r[n]=-(n==1)+2*r[n-1]+(2*r[n-2] if n>=2 else 0)-(r[n-3] if n>=3 else 0)
    if n>=2:r[n]-=sum(r[j]*r[n-2-j] for j in range(n-1))
a=[0]+r[:-1]
p=[1]+[0]*M
T=[0]*(M+1)
for k in range(1,M//2+1):
    p=mul(p,a)
    for n,c in enumerate(p):
        if c:
            for j in range(n,M+1,k):T[j]+=k*k*c
A=mul(T,[1,-2,0,2,-1]+[0]*(M-4))
H=[A[n+2] if n+2<=M else 0 for n in range(M+1)];H[0]+=1
assert H[0]==0 and H[1]==-1
h1=H[1:]+[0]
iH=inv(h1)
gq=[iH[n+1]+(2 if n==0 else 1 if n==1 else 0) for n in range(N+1)]
# Gq=q^-1*((1+q)^2 + inv(H/q)); the q^-1 constants cancel.
qz=[0]+[comb(2*n,n)//(n+1) for n in range(1,M+1)]
zpow=[1]+[0]*M
gz=[0]*(M+1)
for c in gq:
    gz=[a+c*b for a,b in zip(gz,zpow)];zpow=mul(zpow,qz)
def append_ok(w,v,forbidden):
    for j in range(1,len(w)):
        for i in range(j):
            a,b=w[i],w[j]
            # Encode each triple by its three pairwise comparisons.
            key=((a>b)-(a<b),(a>v)-(a<v),(b>v)-(b<v))
            if key in forbidden:return False
    return True
def code(p):
    a,b,c=map(int,p)
    return ((a>b)-(a<b),(a>c)-(a<c),(b>c)-(b<c))
classes={'A':['100','101','110','201'],'C':['100','101','110','210'],'B':['101','110','120','210']}
counts={}
for name,pats in classes.items():
    forbidden=set(map(code,pats));words=[()];totals=[1]
    for n in range(1,N+1):
        words=[w+(v,) for w in words for v in range(n) if append_ok(w,v,forbidden)]
        totals.append(len(words))
    counts[name]=totals
    assert totals==gz[:N+1],(name,totals,gz[:N+1])
output={'N':N,'formula_counts':gz[:N+1],'literal_counts':counts,'all_equal':True,'r_prefix':r[:8],'H_prefix':H[:8]}
Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
