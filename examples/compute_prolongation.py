import sympy as sym
import prolongations as pr


x = sym.symbols('x')
t = sym.symbols('t')
u = sym.Function('u')(x,t)



prolongate = pr.Prolongation([x,t], [u])

phit = prolongate.compute_vector_field_coefficient([1], 0)
phixx = prolongate.compute_vector_field_coefficient([0,0], 0)

d = sym.expand(phit-phixx)

d = d.subs(sym.Function("u_t")(x,t), sym.Function("u_xx")(x,t))

monomials = prolongate.get_monomials(d)



print(monomials)
print(prolongate.output_to_latex(monomials))