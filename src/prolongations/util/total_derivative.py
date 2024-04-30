import sympy as sym

from math import prod


def D_base(f: sym.Function, x: sym.Symbol):
    """Material derivative with simplification to notations"""
    name = f.name
    args = f.args
    function_args = [arg for arg in args if isinstance(arg, sym.Function)]

    x_name = x.name

    # Derivate f with respect to x
    Df_Dx = sym.diff(f, x)

    # Make substitutions
    if function_args == []:
        u_x = sym.Function(f'{name}_{x_name}')(*args)
        Df_Dx = Df_Dx.subs(sym.Derivative(f, x), u_x)

    else:
        for u in function_args:
            u_name = u.name
            u_args = u.args
            
            # u_x
            u_x = sym.Function(f'{u_name}_{x.name}')(*u_args)
            Df_Dx = Df_Dx.subs(sym.Derivative(u, x), u_x)

            # f_u
            f_u = sym.Function(f'{name}_{u_name}')(*args)
            Df_Dx = Df_Dx.subs(sym.Derivative(f, u), f_u)

            # f_x
            f_x = sym.Function(f'{name}_{x_name}')(*args)
            substitutions = {expr: f_x for expr in Df_Dx.atoms(sym.Subs)}
            Df_Dx = Df_Dx.subs(substitutions)


    return Df_Dx


def D(f, x):
    # Df_Dx: sym.Function = None
    if isinstance(f, sym.Function):
        return D_base(f, x)
    elif isinstance(f, sym.Add):
        args = f.args
        return sum(D(arg, x) for arg in args)
    elif isinstance(f, sym.Mul):
        args = f.args
        return sum(D(arg_d, x)*prod([arg_p for arg_p in args if arg_p != arg_d]) for arg_d in args)
    elif isinstance(f, sym.Pow):
        base, exp = f.args
        return exp*sym.Pow(D(base, x), exp-1)
    elif isinstance(f, sym.Integer):
        return 0
    else:
        raise NotImplementedError(f'Function type {type(f)} not implemented') 

