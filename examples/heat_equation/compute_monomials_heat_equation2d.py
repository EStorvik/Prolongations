import sympy as sym
import prolongations as pr


# Create symbolic variables
x = sym.symbols('x')
y = sym.symbols('y')
t = sym.symbols('t')
u = sym.Function('u')(x, y, t)

# Create prolongation object
prolongate = pr.Prolongation([x, y, t], [u])

# Compute the active coefficients in the prolongation applied to the heat equation
phi_t = prolongate.compute_vector_field_coefficient([2], 0)
phi_xx = prolongate.compute_vector_field_coefficient([0,0], 0)
phi_yy = prolongate.compute_vector_field_coefficient([1,1], 0)

d = sym.expand(phi_t-phi_xx-phi_yy)

# Substitute the heat equation
d = d.subs(sym.Function("u_t")(x, y,t), sym.Function("u_xx")(x, y,t)+sym.Function("u_yy")(x, y,t))

# Compute the monomials
monomials = prolongate.get_monomials(d)

# Print the monomials and the latex representation
print(monomials)
print(prolongate.output_to_latex(monomials))