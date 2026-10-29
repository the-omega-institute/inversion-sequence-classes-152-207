#!/usr/bin/env python3
"""Exact integer verification of the Class 207 Lambert/product summation.

Only the Python standard library is needed. Products are expanded directly,
separately at every orbit point; the summation identity is not used to build D.
The C legal-site recurrence and B four-type tree are imported from the supplied
earlier generator (not claimed as independent implementations).
"""
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parent


class Ring:
    def __init__(self, n): self.n = n
    def c(self, a): return [a] + [0]*self.n
    def sh(self, a, k):
        if k >= 0: return ([0]*k+a)[:self.n+1]
        assert all(x == 0 for x in a[:-k])
        return a[-k:]+[0]*(-k)
    def add(self, a, b, scale=1): return [x+scale*y for x,y in zip(a,b)]
    def mul(self, a, b):
        out = self.c(0)
        aa = [(i,x) for i,x in enumerate(a) if x]
        bb = [(j,y) for j,y in enumerate(b) if y]
        for i,x in aa:
            for j,y in bb:
                if i+j > self.n: break
                out[i+j] += x*y
        return out
    def div(self, a, b):
        assert b[0] != 0
        out = self.c(0)
        for i in range(self.n+1):
            x = a[i]-sum(b[j]*out[i-j] for j in range(1,i+1))
            out[i], rem = divmod(x,b[0])
            assert rem == 0, (i, b[0])
        return out
    def om(self, a): return self.add(self.c(1),a,-1)
    def root(self):
        r = self.c(0)
        for n in range(1,self.n+1):
            r[n] = -(n == 1)
            for j,c in [(1,2),(2,2),(3,-1)]:
                if n >= j: r[n] += c*r[n-j]
            if n >= 2: r[n] -= sum(r[i]*r[n-2-i] for i in range(n-1))
        return r
    def U(self, v, r):
        # q/r=-L-q^2*r is regular; no inversion of a zero constant.
        L = self.c(1); L[1:4] = [-2,-2,1]
        qr = self.add([-x for x in L],self.sh(r,2),-1)
        a = self.mul(qr,v); b = self.sh(self.mul(r,v),2)
        c = self.sh(v,1); out = self.c(1)
        for j in range(self.n+1):
            out = self.mul(out,self.om(self.sh(a,j)))
            out = self.mul(out,self.om(self.sh(b,j)))
            den = self.om(self.sh(c,j))
            out = self.div(self.div(out,den),den)
        return out
    def W(self,v,r):
        fac = self.c(1); fac[1:3]=[-2,1]
        return self.div(self.mul(self.mul(fac,v),self.U(v,r)),self.om(self.sh(v,2)))
    def T(self,r):
        out=self.c(0); power=self.c(1); qr=self.sh(r,1)
        for k in range(1,self.n//2+1):
            power=self.mul(power,qr)
            for j in range(0,self.n+1,k):
                out=self.add(out,self.sh(power,j),k*k)
        return out
    def compose(self, coeffs, z):
        out=self.c(0); power=self.c(1)
        for c in coeffs:
            out=self.add(out,power,c); power=self.mul(power,z)
        return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--order',type=int,default=128)
    p.add_argument('--counts',type=int,default=128)
    args=p.parse_args(); start=perf_counter()
    N=args.order; A=Ring(N+8); r=A.root(); T=A.T(r)
    alpha=A.c(1); alpha[1:5]=[-2,0,2,-1]  # (1-q)^3(1+q)
    H=A.add(A.c(1),A.sh(A.mul(alpha,T),-2))
    M=[]
    for a in (A.c(1),r):
        m=A.c(0)
        for k in range(A.n+1):
            v=A.sh(a,k)
            if any(v): m=A.add(m,A.W(v,r))
        M.append(m)
    D=A.add(M[0],M[1],-1)
    U1=A.U(A.c(1),r)
    fac=A.c(1);fac[2]=-1
    qJ=A.mul(fac,U1)
    residual=A.add(A.sh(D,1),A.mul(qJ,H))
    assert residual[:N+2] == [0]*(N+2)
    Hunit=H[1:]+[0]
    inv=A.div(A.c(1),Hunit)
    assert inv[0] == -1
    # 1/H+q^-1+2+q; the q^-1 terms cancel.
    G=inv[1:]+[0]; G[0]+=2; G[1]+=1
    from r3_counts import c_coefficients,b_totals
    Cmax=args.counts
    F,_,_=c_coefficients(Cmax)
    C=[1]+[sum(x.values()) for x in F[1:]]
    B=b_totals(Cmax)
    assert C==B
    z=A.c(0)
    for n in range(1,A.n+1):z[n]=(-1)**(n-1)*n
    gq=A.compose(C,z)
    checked=min(N,Cmax)
    assert G[:checked+1]==gq[:checked+1]
    # A regenerable complete table is deliberately not attached.
    table='n,C,B\n'+''.join(f'{i},{c},{b}\n' for i,(c,b) in enumerate(zip(C,B)))
    (ROOT/'r3-counts.csv').write_text(table)
    data={'order':N,'internal_order':A.n,'product_residual_zero_through':N+1,
          'total_formula_matches_both_recurrences_through':checked,
          'C_B_totals_equal_through':Cmax,
          'r_prefix':r[:9],'H_prefix':H[:9],'D_prefix':D[:9],
          'total_prefix':C[:13], 'table_sha256':hashlib.sha256(table.encode()).hexdigest(),
          'exact_integer_divisions':True}
    (ROOT/'r3-exact-check.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2));print('elapsed_seconds',round(perf_counter()-start,2))


if __name__=='__main__': main()
