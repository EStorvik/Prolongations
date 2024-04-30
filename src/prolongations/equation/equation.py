import prolongations as pr
import sympy as sym



class Equation():

    def __init__(self, equation: sym.Eq) -> None:
        self.equation = equation

    def __call__(self):
        pass