import sympy as sym
import prolongations as pr


# Create symbolic variables
x = sym.symbols('x')
t = sym.symbols('t')
u = sym.Function('u')(x,t)


# Create symbolic derivatives in correct style
u_t = sym.Function('u_t')(x,t)
u_x = sym.Function('u_x')(x,t)
u_xx = sym.Function('u_xx')(x,t)
u_xxxx = sym.Function('u_xxxx')(x,t)


# Create prolongation object
prolongate = pr.Prolongation([x,t], [u])

# Compute the active coefficients in the prolongation applied to the heat equation
phit = prolongate.compute_vector_field_coefficient([1], 0)
phi = prolongate.phi[0]
phix = prolongate.compute_vector_field_coefficient([0], 0)
phixx = prolongate.compute_vector_field_coefficient([0,0], 0)
phixxxx = prolongate.compute_vector_field_coefficient([0,0,0,0], 0)

d = sym.expand(phit-phi*(24*u_x**2+24*u*u_xx-12*u_xx)-phix*(48*u*u_x-24*u_x) - phixx*(12*u**2-12*u+2)+phixxxx)

# Substitute the heat equation
d = d.subs(u_t,24*u*u_x**2+12*u*2*u_xx-12*u_x**2-12*u*u_xx+2*u_xx-u_xxxx)

# Compute the monomials
monomials = prolongate.get_monomials(d)

# Print the monomials and the latex representation
print(monomials)
print(prolongate.output_to_latex(monomials))