from sympy import *

alpha, beta, gamma, Ni, Nj = symbols('alpha beta gamma Ni Nj')
x= symbols('x', real=True)
S=Ni*Nj*exp(beta**2/(4*alpha)+gamma)*sqrt(pi/alpha)

x_expec=diff(S, beta)
x2_expec=diff(x_expec, beta)
x3_expec=diff(x2_expec, beta)
x4_expec=diff(x3_expec, beta)
x5_expec=diff(x4_expec, beta)
x6_expec=diff(x5_expec, beta)

x_expec/=S
print(expand(simplify(x_expec)))
x2_expec/=S
print(expand(simplify(x2_expec)))
x3_expec/=S
print(expand(simplify(x3_expec)))
x4_expec/=S
print(expand(simplify(x4_expec)))
x5_expec/=S
print(expand(simplify(x5_expec)))
x6_expec/=S
print(expand(simplify(x6_expec)))