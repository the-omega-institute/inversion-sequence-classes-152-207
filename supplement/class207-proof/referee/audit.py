"""Independent exact and numerical referee audit; does not import production code."""
import json
from pathlib import Path
from collections import Counter
import sympy as s
import mpmath as mp
q,r,u,v,h,z,x,t=s.symbols('q r u v h z x t')
L=1-2*q-2*q*q+q**3
relation=q+L*r+q*q*r*r
checks={}
def zero(name,expr,relation=False):
    n=s.cancel(expr).as_numer_denom()[0]
    if relation:n=s.rem(n,s.Poly(q+L*r+q*q*r*r,r),r)
    val=s.factor(n)
    checks[name]=val==0
    assert val==0,(name,val)
D=q*u*u+L*u+q*q
jG=(u-1)*(u-q)*(u*u-q)*(1-q)**2/(u*D)
j0=(u-1)**2*(u-q)**2*(u*u-q)*(1+q)/(u*u*D)
zero('C_Lambert_redundancy',j0+(1+q)**2/q*jG-(1+q)/q*(u-1)*(u-q)*(u*u-q)/(u*u))
zero('C_partial_fractions',jG-(1-q)**2/q*(u-q/u)-(1-q)**4*(1+q)/q**2*(1/(1-r*u)-1/(1-q*r/u)),True)
beta=(1-q)*(1-r)*(1-q*r)/(q*r)
alpha=(1-q)**3*(1+q)/q**2
zero('beta_specialization',beta+alpha,True)
zero('B_telescoper',(1-q*v)**3-(1-q*q*v)*(1-v/r)*(1-q*r*v)+(1-q)**2/q*v*(1-q**3*v),True)
b1=-r*h*(1-q-h)**2/((1-q)*(1-q*r-h)*(1-r*(1-h)**2))
R=(1-r+r*h)*(1-q*r)*(1-h)/((1-q*r-h)*(1-r*(1-h)**2))
zero('universal_cubic_R',R-1-b1-r*q*h**3/((1-q)*(1-q*r-h)*(1-r*(1-h)**2)))
# Product shifts: I_{a,b}=(q^a r^b;q)_infty. Reduce a to 0 for b!=0,
# to 1 for b==0; retain a formal symbol for each independent base product.
def canonical(factors):
    scalar=s.Integer(1); powers=Counter()
    for a,b,e in factors:
        base=0 if b else 1
        powers[b]+=e
        if a>base:
            for k in range(base,a):scalar/=(1-q**k*r**b)**e
        elif a<base:
            for k in range(a,base):scalar*=(1-q**k*r**b)**e
    assert all(c==0 for c in powers.values()),powers
    return s.factor(scalar)
# Wr Aminus /J1: include scalar separately, then finite product quotient.
prod=canonical([(1,0,1),(2,2,1),(1,1,-2),(-2,-1,1),(1,1,2),(0,1,1),(2,2,-1),(1,0,-3),(0,-1,-1),(1,1,-1),(1,0,2)])
ratio=(1-q)**2*r/(1-q*q*r)*q*q/(1-q)**2*(1-q)*prod
zero('product_constant_23',ratio-beta)
# Aminus*b2/h^3=-W1/Wr. Multiply both sides by Wr and divide W1.
b2=(1-1/q)**2*(1-q)*q*q/((1-q)*(1-q*q)*(1-1/(q*r))*(1-1/r)*(1-r)*(1-q*r))
prod=canonical([(-2,-1,1),(1,1,2),(0,1,1),(2,2,-1),(1,0,-3),(1,0,1),(2,2,1),(1,1,-2),(1,-1,-1),(2,1,-1),(1,0,2)])
scalar=q*q/(1-q)**2*(1-q)**2*r/(1-q*q*r)/((1-q)**2/(1-q*q))
zero('tail_n2',scalar*prod*b2+1)
# General tail ratio, with Q=q^n; compare to W(q^{n-1})/W(q^{n-2}).
Q=s.symbols('Q')
tail=q*(1-Q/q)**2*(1-Q)/((1-q*Q)*(1-Q/(q*r))*(1-r*Q))
w_ratio=q*(1-q*(Q/q**2))**2*(1-q*q*(Q/q**2))/((1-q**3*(Q/q**2))*(1-q*(Q/q**2)/r)*(1-q*q*r*(Q/q**2)))
zero('tail_all_n_ratio',tail-w_ratio)
# cubic logarithmic multipliers: 3+(-1)^j-2^j, j=1,2,3
checks['logarithmic_multipliers']=[3+(-1)**j-2**j for j in (1,2,3)]==[0,0,-6]
# B published-tree kernel elimination, independent rational certificate.
f,g,F=s.symbols('f g F');y=1/(1-z*x)
expr=((1-z*x)*f-z*x-z*(f-g))*(x-y)-z*y*(f-F)
old=z*F-(1-z*x)*(1-x+z*x+z*x*x)*f-z*(1-x+z*x*x)*(g-x)
zero('published_tree_elimination',(1-z*x)*expr-old)
# Independent two-parameter numerical theorem (not on quadratic branch).
mp.mp.dps=70
numeric=[]
for qs,rs in [('0.12','-0.31'),('0.27','-0.63'),('0.4','-0.2')]:
    qq,rr=mp.mpf(qs),mp.mpf(rs)
    def U(vv):return mp.qp(qq*vv/rr,qq)*mp.qp(qq*qq*rr*vv,qq)/mp.qp(qq*vv,qq)**2
    def W(vv):return (1-qq)**2*vv*U(vv)/(1-qq*qq*vv)
    M1=mp.fsum(W(qq**k) for k in range(190))
    Mr=mp.fsum(W(rr*qq**k) for k in range(190))
    TT=mp.fsum(k*k*(qq*rr)**k/(1-qq**k) for k in range(1,220))
    JJ=mp.qp(1/rr,qq)*mp.qp(qq*rr,qq)/((1-qq)*mp.qp(qq,qq)**2)
    bb=(1-qq)*(1-rr)*(1-qq*rr)/(qq*rr)
    residual=M1-Mr+JJ*(1-bb*TT)
    numeric.append({'q':qs,'r':rs,'absolute_residual':mp.nstr(abs(residual),8)})
    assert abs(residual)<mp.mpf('1e-60')
output={'exact_certificates':checks,'numerical_independent_parameters':numeric,'numerical_scope':'70 decimal working precision, 190 orbit terms and 219 Lambert terms; supplementary to the written proof'}
Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
