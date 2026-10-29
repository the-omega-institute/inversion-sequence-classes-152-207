import json
import sympy as s
x,y,z,w,q,b,G,Px,Py,Fx,Fy,b0,t,L=s.symbols('x y z w q b G Px Py Fx Fy b0 t L')
checks=[]
def check(name,e):
 e=s.factor(e)
 if e!=0:raise AssertionError((name,e))
 checks.append(name)
D=1-x+z*x*x;A=D+z*x;M=z*x+z-1;phi=1/(1-z*x)
P=lambda v,F: z*v/(1-z*v)*(1+(F-b0)/(v-1))
V=(x-y)*(Fx-P(x,Fx))/(z*x)+(1-y)*(P(x,Fx)-P(y,Fy))+(x-y)*(Fx-b0)/(x-1)-Fx+Fy
scalar=M*A*Fx+z*x*(1-z*x)**2*Fy-x*z*z*D*(b0+1-1/z)
check('B scalar elimination',V.subs(y,phi)+scalar/(z*x*(1-z*x)*M))
X=lambda v:(1+q)*(1-v)/(1-q*v)
zz=q/(1+q)**2
check('B conjugacy',phi.subs({x:X(w),z:zz})-X(q*w))
Rraw=(x-1)*A/(z*z*x**3)
R=(1-w/q)*(1-w/(q*b))*(1-q*q*b*w)/(1-w)**3
bpoly=q**3*b*b+(1-2*q-2*q*q+q**3)*b+1
num=s.together(R-Rraw.subs({x:X(w),z:zz})).as_numer_denom()[0]
check('B multiplier factorization',s.rem(num,bpoly,b))
Braw=D*(G-x)/(z*x*x*(1-z*x))
B0=w*(q-1)**2*(q+1)*(G*q*w-G-q*w+q-w+1)/(q*(w-1)**2*(q*q*w-1))
check('B forcing',Braw.subs({x:X(w),z:zz})-B0)
c=q/(1-q)**2
C=lambda v:1-c*(v+1/v-2)
Y=(q-w)*(q*q*w-1)/(w*(q-1)**2*(q+1))
kappa=(1-y)*(x-y)+z*x*(x-1)
check('C first kernel root',kappa.subs({x:C(w),y:Y,z:zz}))
check('C second kernel root',kappa.subs({x:C(q*w),y:Y,z:zz}))
check('C root-sum parameter',Y-1-zz*(C(w)+C(q*w)-1))
# Original B cleared numerator, substituting F(t) affine in P(t).
F0,P0,fp=s.symbols('F0 P0 fp')
Ft=G-t-(1-t)*(1-z*t)*fp/(z*t)
cleared=F0-P0+z*x*(1-t)*(P0-fp)/(x-t)+z*x*((F0-G+1)/(x-1)-(F0-Ft)/(x-t))
AA=F0-P0+z*x*(1-t)*P0/(x-t)+z*x*(F0-G+1)/(x-1)+z*x*(G-t-F0)/(x-t)
check('full B transfer',cleared-AA+x*(1-t)*fp/(t*(x-t)))
# Class152 rational identities, independently in t with z=t(1-t).
cat=1/(1-t);zc=t*(1-t);q1=1+zc;q2=1/(1-zc)
Dq=lambda a:cat*t*t*(a-1)/(1-a*t)
d1=t**3/((1-t)*(1-t*t));d2=t**3/(1-t)**2
check('Class152 first weighted run',Dq(q1)-d1)
check('Class152 second weighted run',Dq(q2)-d2)
extra=d1*(t/(1-2*t)+t*t/(1-zc))
check('Class152 weighted left simplification',extra-cat*t**4/((1-zc)*(1-2*t)))
left=cat+extra
right=1+zc/(2*(1-zc))*(cat+1/(1-t/(1-2*zc)))
final=1+zc/(2*(1-zc))*(1+1/(1-2*t))
check('Class152 equality',left-right)
check('Class152 closed generating function',right-final)
# Corrected fixed point for the finite kernel recurrence.
alpha=z*x/(A*(x-1));beta=x*(x*x*z*z-2*x*z-z+1)/((x*z-1)*A*(x-1));gamma=x*z*(x*x*z-1)/((x*z-1)*A*(x-1))
zL=(L-1)/L**2
check('Catalan fixed-point solution',((beta*b0+gamma)/(1-alpha)).subs({x:L,z:zL})-(b0+1-L)/(zL*L))
# Bridge between the formal and analytic C parametrizations.
T=lambda v:v*(z*v-v+1)/(z*v*v-v+1)
check('C parameter bridge',C(w)/Y-X(q*w))
check('C first bridge root',T(X(q*w)).subs(z,zz)-C(w))
check('C second bridge root',T(phi.subs(x,X(q*w))).subs(z,zz)-C(q*w))
print(json.dumps({'exact_identities':checks},indent=2))
