#!/usr/bin/env python3
"""Denominator-cleared symbolic certificates; requires SymPy.

Infinite products are canceled by finite shift relations, not evaluated.
Every residual below is an exact rational identity, reduced modulo E(r)
only when explicitly indicated.
"""
import json
from pathlib import Path
import sympy as s

q,r,u,v,h=s.symbols('q r u v h')
L=1-2*q-2*q*q+q**3; E=q+L*r+q*q*r*r
checks={}


def check(name,expr,branch=False):
    num,den=s.fraction(s.cancel(expr))
    if branch: num=s.rem(s.Poly(num,r),s.Poly(E,r)).as_expr()
    assert s.expand(num)==0,(name,s.factor(num))
    checks[name]={'residual':'0','quadratic_relation_used':branch}


# Canonicalize (q^j r^e;q)_infty to a formal base product times
# a finite rational expression. For e=0 use (q;q)_infty as base.
def poch(j,e):
    base=1 if e==0 else 0
    out=s.Symbol('P'+str(e))
    if j>=base:
        for k in range(base,j): out/=1-q**k*r**e
    else:
        for k in range(j,base): out*=1-q**k*r**e
    return out


D=q*u*u+L*u+q*q
jg=(u-1)*(u-q)*(u*u-q)*(1-q)**2/(u*D)
j0=(u-1)**2*(u-q)**2*(u*u-q)*(1+q)/(u*u*D)
check('Lambert_redundancy',j0+(1+q)**2/q*jg-(1+q)/q*(u-1)*(u-q)*(u*u-q)/(u*u))
partial=(1-q)**2/q*(u-q/u)+(1-q)**4*(1+q)/q**2*(1/(1-r*u)-1/(1-q*r/u))
check('partial_fraction',jg-partial,True)
beta=(1-q)*(1-r)*(1-q*r)/(q*r)
alpha=(1-q)**3*(1+q)/q**2
check('beta_specialization',beta+alpha,True)
check('telescoper',(1-q*v)**3-(1-q*q*v)*(1-v/r)*(1-q*r*v)+(1-q)**2/q*v*(1-q**3*v),True)
b1=-r*h*(1-q-h)**2/((1-q)*(1-q*r-h)*(1-r*(1-h)**2))
R=(1-r+r*h)*(1-q*r)*(1-h)/((1-q*r-h)*(1-r*(1-h)**2))
check('universal_cubic_identity',R-1-b1-r*q*h**3/((1-q)*(1-q*r-h)*(1-r*(1-h)**2)))

Wr=(1-q)**2*r/(1-q*q*r)*poch(1,0)*poch(2,2)/poch(1,1)**2
W1=(1-q)**2/(1-q*q)*poch(1,-1)*poch(2,1)/poch(1,0)**2
J1=poch(0,-1)*poch(1,1)/((1-q)*poch(1,0)**2)
Am=q*q/(1-q)**2*poch(-2,-1)*poch(1,1)**2*poch(0,1)/(poch(2,2)*poch(1,0)**3)
b2=(1-1/q)**2*q*q/((1-q*q)*(1-1/(q*r))*(1-1/r)*(1-r)*(1-q*r))
check('product_constant',Wr*Am/J1-beta)
check('tail_base_n2',Wr*Am*b2/W1+1)
wr_ratio=q*(1-q*v)**2*(1-q*q*v)/((1-q**3*v)*(1-q*v/r)*(1-q*q*r*v))
tail_ratio=q*(1-q*v)**2*(1-q*q*v)/((1-q**3*v)*(1-q*v/r)*(1-r*q*q*v))
check('tail_all_index_ratio',tail_ratio-wr_ratio)
check('J1_branch',J1/W1-(1+q)**2/q,True)
# Removing finite shifts in P(h)/R(h) leaves precisely the four products
# with exponents 3,1,-3,-1 in the confluence calculation.
t=1-h
finite_PoverS=(1-r*t)*(1-q*r)/((1-r*t*t)*(1-q*r/t))
check('P_over_R_product_shift',finite_PoverS-R)
check('cubic_constant',r*q/((1-q)*(1-q*r)*(1-r))-1/beta)
z,t,G,P1,P2=s.symbols('z t G P1 P2')
d=t*t-t+z; X=t*(t-1)/d; Xp=(z-t)/d; Y=(t-1)*(t-z)/d
kap=lambda x,y:(1-y)*(x-y)+z*x*(x-1)
check('C_kernel_first_root',kap(X,Y))
check('C_kernel_second_root',kap(Xp,Y))
tw=(2-Y)*((1-X)**2*P1-(1-Xp)**2*P2)-(1-Y)*(X-Xp)*(G+(1-Y)/z)
target=(d+t*z)*(z*z*P1-t**4*P2)-z*t*(t*t-z)*(d*G+t)
ratio=s.cancel(s.diff(tw,P1)/s.diff(target,P1))
check('C_two_root_binding',tw-ratio*target)
tq=(u-q)/((1+q)*(u-1));zq=q/(1+q)**2
dq=tq*tq-tq+zq
check('C_additive_substitution',tq*(tq*tq-zq)*(dq*G+tq)/(zq*(dq+tq*zq)*dq*dq)-(1+q)**4/(q*(1-q)**3)*(jg*G+j0))
xv=(1+q)*(1-v)/(1-q*v)
av=(1-q*q*v)*(q+L*v+q*q*v*v)/(q*(1-q*v)**3)
bv=(1-q)**2*v/(1-q*v)**2
check('B_tree_qdiff_multiplier',(1-zq*xv)*(1-xv+zq*xv+zq*xv*xv)/zq-av)
check('B_tree_qdiff_forcing',1-xv+zq*xv*xv-bv)
Jratio=(1-q*v)**3/((1-q*q*v)*(1-v/r)*(1-q*r*v))
check('B_J_normalization',Jratio*av-1,True)
check('B_W_normalization',Jratio*bv-(1-q)**2*v*(1-q*v)/((1-q*q*v)*(1-v/r)*(1-q*r*v)),True)
assert [3+(-1)**j-2**j for j in (1,2,3)]==[0,0,-6]
checks['logarithmic_derivative_multipliers']=[0,0,-6]
data={'certificates':checks,'all_passed':True,
      'branch_polynomial':'q+(1-2*q-2*q**2+q**3)*r+q**2*r**2'}
Path(__file__).with_name('r3-symbolic-check.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
