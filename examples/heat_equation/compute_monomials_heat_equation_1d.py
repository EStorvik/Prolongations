import sympy as sym
import prolongations as pr


# Create symbolic variables
x = sym.symbols('x')
t = sym.symbols('t')
u = sym.Function('u')(x,t)

# Create prolongation object
prolongate = pr.Prolongation([x,t], [u])

# Compute the active coefficients in the prolongation applied to the heat equation
phit = prolongate.compute_vector_field_coefficient([1], 0)
phixx = prolongate.compute_vector_field_coefficient([0,0], 0)

d = sym.expand(phit-phixx)

# Substitute the heat equation
d = d.subs(sym.Function("u_t")(x,t), sym.Function("u_xx")(x,t))

# Compute the monomials
monomials = prolongate.get_monomials(d)

# Print the monomials and the latex representation
print(monomials)
print(prolongate.output_to_latex(monomials))