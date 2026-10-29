"""Independent manuscript-to-two-root rational certificate."""
import sympy as s
import json
from pathlib import Path
z,x,y,t,P,R,Q,G,Z=s.symbols('z x y t P R Q G Z')
Fx=G-x-(1-x)*(1-z*x)/(z*x)*P
Fy=G-y-(1-y)*(1-z*y)/(z*y)*R
rhs=(x-y)*(Fx-P)/(z*x)-(1-y)*Z+(1-y)*(P-R)+(x-y)*(Fx-G+1)/(x-1)-Fx+Fy
Zsol=s.solve(rhs,Z)[0]
expr=(1-y)**2*(Zsol+R/(z*y))-(1-x)*(G-x)-(2-y)*(1-x)**2/(1-y)*P
num=s.cancel(expr).as_numer_denom()[0]
k=(1-y)*(x-y)+z*x*(x-1)
remainder=s.factor(s.rem(num,k,x));assert remainder==0
d=t*t-t+z;X=t*(t-1)/d;Xh=(z-t)/d;Y=(t-1)*(t-z)/d
Ctwo=(2-Y)*((1-X)**2*P-(1-Xh)**2*Q)-(1-Y)*(X-Xh)*(G+(1-Y)/z)
Input=(d+t*z)*(z*z*P-t**4*Q)-z*t*(t*t-z)*(d*G+t)
ratio=s.factor(s.diff(Ctwo,P)/s.diff(Input,P))
res=s.factor(Ctwo-ratio*Input);assert res==0
output={'C_diagonal_remainder_mod_kernel':str(remainder),'two_root_proportionality_factor':str(ratio),'two_root_residual':str(res)}
Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
