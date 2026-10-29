#!/usr/bin/env python3
"""Exact checks of signed catalysts; independent literal site reconstruction."""
from collections import Counter
from itertools import product
import json
from opposite_catalysts import coefficients, anti

def rank3(a,b,c):
    vals=sorted({a,b,c})
    return tuple(vals.index(v) for v in (a,b,c))
BAD={(1,0,0),(1,0,1),(1,1,0),(2,1,0)}
def legal(w,v):
    return all(rank3(w[i],w[j],v) not in BAD for j in range(len(w)) for i in range(j))
def literal(N,Z):
    words=[()]; sizes=[1]
    for n in range(1,N+1):
        words=[w+(v,) for w in words for v in range(n) if legal(w,v)]
        sizes.append(len(words)); actual=Counter()
        for w in words:
            M=max(w); sites=[v for v in range(M+1) if legal(w,v)]
            for v in sites:
                if v==M:continue
                k=sum(t>v and t not in w for t in sites)
                actual[n-M,k]+=1
        assert actual==Z[n],(n,actual,Z[n])
    return sizes

def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[0],-a[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def power(a,n):
    p=(1,0)
    for _ in range(n):p=mul(p,a)
    return p

def evaluate(p,a,b):
    v=(0,0)
    for (i,j),c in p.items():v=add(v,mul((c,0),mul(power(a,i),power(b,j))))
    return v

def main():
    F,P,Z=coefficients(61)
    sizes=literal(8,Z)
    failures=[(n,{k:2*v for k,v in anti(p).items() if k%2==0}) for n,p in enumerate(Z)]
    failures=[(n,p) for n,p in failures if p]
    assert failures[0]==(5,{2:2})
    # x=1,-1 and y=1,-1,i; Gaussian integer arithmetic throughout.
    ys=[(1,0),(-1,0),(0,1)];betas=[(-1,-1),(-1,1),(2,0)]
    def W(n,a,b):
        K0=mul(add((1,0),neg(b)),add(a,neg(b)))
        # a is +1 or -1, so its inverse equals a.
        return add(mul(mul(K0,a),evaluate(Z[n+1],a,b)),mul(add(a,(-1,0)),evaluate(Z[n],a,b)))
    for n in range(61):
        v=(0,0)
        for b,beta in zip(ys,betas):v=add(v,mul(beta,add(W(n,(1,0),b),neg(W(n,(-1,0),b)))))
        assert v==(0,0),(n,v)
    print(json.dumps({'literal_C_counts_through_8':sizes,'whole_bivariate_Z_agrees_through':8,'signed_first_failure':{'z_degree':5,'x_degree':2,'coefficient':2},'six_term_fourth_root_zero_through_z':60,'method':'exact Gaussian integers; 61 C-system coefficients; independent literal forbidden-triple sites through length 8'},indent=2))
if __name__=='__main__':main()
