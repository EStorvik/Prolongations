import sympy as sym



class Prolongation():

    def __init__(self, x, u):
        self.x = x
        self.u = u

        self.monomials = {}
        self.ksi = sym.symbols("ksi")
        self.phi = sym.symbols("phi")

        self.v = [self.ksi, self.phi]

    def material_derivative(self, ksi):
        return 

        

    def prolongation(self, f, x, y, z):
        pass




def material_derivative(f,x):
    return f_x + u_x*f_u

x = sym.symbols('x')
y = sym.symbols('y')

u = sym.symbols('u')

k1 = sym.symbols('k1')
k2 = sym.symbols('k2')




