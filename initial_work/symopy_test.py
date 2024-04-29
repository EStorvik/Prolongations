import sympy as sym

x = sym.symbols('x')
y = sym.symbols('y')

u = sym.Function('u')(x, y)
u_x = sym.Function('u_x')(x, y)
du_dx = sym.diff(u, x)

f = sym.Function('f')(x, y, u)
f_sym = sym.Function('f_sym')(x, y, u)

ksi = sym.Function('ksi')(x, y)
u_sym = u.subs(sym.Subs(u, (x, y), (ksi, y)))

df_du = sym.Derivative(f, u)
f_u = sym.Function('f_u')(x, y, u)

f_x = sym.Function('f_x')(x, y, u_sym)

Df_Dx = sym.diff(f, x)

Df_Dx = Df_Dx.subs(df_du, f_u)
Df_Dx = Df_Dx.subs(du_dx, u_x)
Df_Dx = Df_Dx.subs(sym.Subs(sym.Derivative(f_sym, x), (x,), (x,)))

print(Df_Dx)
