import sympy as sp

a_i,b_i,mu_i,p_i= sp.symbols('a_i b_i mu_i p_i', real=True)
x_i=sp.symbols('x_i')
I=sp.I
pi=sp.pi
a_j,b_j,mu_j,p_j= sp.symbols('a_j b_j mu_j p_j', real=True)
a_k,b_k,mu_k,p_k= sp.symbols('a_k b_k mu_k p_k', real=True)
x=sp.symbols('x', real=True)

alpha_i=a_i**2+I*b_i
alpha_j=a_j**2+I*b_j
alpha_k=a_k**2+I*b_k
polynomial=sp.conjugate(-alpha_i*(x-mu_i)**2+I*p_i*(x-mu_i))
polynomial+=(-alpha_j*(x-mu_j)**2+I*p_j*(x-mu_j))
polynomial+=(-alpha_k*(x-mu_k)**2+I*p_k*(x-mu_k))

coeffs = sp.Poly(polynomial, x).all_coeffs()
A, B, C = coeffs  # quadratic polynomial: A*x^2 + B*x + C
A=sp.simplify(A)
B=sp.simplify(B)
C=sp.simplify(C)
print("A=",A)
print("B=",B)
print("C=",C)
mu_ijk=sp.simplify(-B/(2*A))
d=sp.simplify(C - B**2/(4*A))
print("mu_ijk=",sp.simplify(-B/(2*A)))
print("d=",sp.simplify(C - mu_ijk**2*A))
#Make sure that the two expressions are equivalent
polynomial_2=A*(x-mu_ijk)**2+d
print(sp.simplify(polynomial-polynomial_2))