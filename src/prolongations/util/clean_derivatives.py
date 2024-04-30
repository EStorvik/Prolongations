import sympy as sym

from math import prod

def sort_derivatives(derivatives, args):
    order = [arg.name for arg in args]
    sorted_derivatives = []
    for arg in order:
        for derivative in derivatives:
            if arg == derivative:
                sorted_derivatives.append(derivative)
    return sorted_derivatives


def clean_base(f):
    name_list = f.name.split('_')
    derivatives = name_list[1:]
    # Sort list of derivatives to ensure same order as args
    sorted_derivatives = sort_derivatives(derivatives, f.args)
    # Create new name
    new_name = f'{name_list[0]}_{"".join(sorted_derivatives)}'
    return sym.Function(new_name)(*f.args)

def clean(f):
    if isinstance(f, sym.Function):
        return clean_base(f)
    elif isinstance(f, sym.Add):
        args = f.args
        return sum(clean(arg) for arg in args)
    elif isinstance(f, sym.Mul):
        args = f.args
        return prod([clean(arg) for arg in args])
    elif isinstance(f, sym.Pow):
        base, exp = f.args
        return sym.Pow(clean(base), exp)
    elif isinstance(f, sym.Integer):
        return f
    else:
        raise NotImplementedError(f'Function type {type(f)} not implemented')