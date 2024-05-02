from future import __annotations__
import prolongations as pr
import sympy as sym

class Prolongation():
    """
    Prolongation class to compute vector field coefficients.

    Takes us to the prolongation of the given system of PDEs and computes the vector field coefficients. 
    There is still some manual work to be done in terms of translating the PDE, but this is WIP.
    """
    def __init__(self, independent_variables: list[sym.Symbol], dependent_variables: list[sym.Function]) -> None:
        """
        Initialize the Prolongation object

        args:
            independent_variables: list[sym.Symbol] - list of independent variables
            dependent_variables: list[sym.Function] - list of dependent variables

        attributes:
            independent_variables: list[sym.Symbol] - list of independent variables
            dependent_variables: list[sym.Function] - list of dependent variables
            xi: list[sym.Function] - list of xi functions
            phi: list[sym.Function] - list of phi functions

        """
        self.independent_variables = independent_variables
        self.dependent_variables = dependent_variables

        self.xi, self.phi = self._init_vector_field_coefficients()


    def get_monomials(self, expr) -> dict:
        """
        Computes a disctionary of monomials and their coefficients from a given expression.

        Usefull for computing the vector field coefficients. This is the main important method.
        """

        # Init output dictionary
        monomials = dict()

        # Simplify and expand the expression
        expr = sym.expand(expr)

        # get the terms of the expression
        terms = expr.args

        # Loop through all the terms and find the monomials
        for term in terms:
            # If it is a multiplication chech what is multiplied with a monomial
            if isinstance(term, sym.Mul):
                monomial_type =''
                monomial_coefficient = 1
                factors = term.args
                for factor in factors:
                    # Add integer factors to the coefficient
                    if isinstance(factor, sym.Integer):
                        monomial_coefficient *= factor
                    elif isinstance(factor, sym.Function):
                        # Add monomial factors to the monomial type
                        if self._is_monomial(factor.name):
                            monomial_type += factor.name
                        # Otherwise, add to the coefficient
                        else:
                            monomial_coefficient *= factor
                    # If power, check if monomial or coefficient and add to respective places
                    elif isinstance(factor, sym.Pow):
                        base, exp = factor.args
                        if self._is_monomial(base.name):
                            monomial_type += f'{base.name}^{exp}'
                        else:
                            monomial_coefficient *= factor
                    else:
                        raise NotImplementedError(f'Factor type {type(factor)} not implemented')
                # Add the data to the dictionary
                if monomial_type not in monomials:
                    monomials[monomial_type] = monomial_coefficient
                else:
                    monomials[monomial_type] += monomial_coefficient

            # Special case where it is just one term
            elif isinstance(term, sym.Function):
                if self._is_monomial(term.name):
                    monomial_type = term.name
                    monomial_coefficient = 1
                    if monomial_type not in monomials:
                        monomials[monomial_type] = monomial_coefficient
                    else:
                        monomials[monomial_type] += monomial_coefficient
                else:
                    if '' not in monomials:
                        monomials[''] = term
                    else:
                        monomials[''] += term
            else:
                raise NotImplementedError(f'Term type {type(term)} not implemented')
            
        return monomials


    # Check if function has the same name as (or is a derivative of) a monomial
    def _is_monomial(self, name):
        name_list = name.split('_')
        for var in self.dependent_variables:
            if name_list[0] == var.name:
                return True
            else:
                return False



    def compute_vector_field_coefficient(self, derivative: list[int], dependent_variable_index: int):
        """
        Computes the vector field coefficient for a given derivative list of derivatives and dependent variable index.

        args:
            derivative: list[int] - list of derivatives in index form. Uses index in the self.independent_variables list.
            dependent_variable_index: int - index of the dependent variable in the self.dependent_variables list.
        """
        assert len(self.dependent_variables)>dependent_variable_index; "The dependent variable index is out of range"

        # Get the dependent variable
        u: sym.Function = self.dependent_variables[dependent_variable_index]
        
        # Compute the vector field coefficient
        body = self.phi[dependent_variable_index]
        for i, x in enumerate(self.independent_variables):
            body -= self.xi[i]*pr.D(u, x)
        phi = pr.D(body,self.independent_variables[derivative[0]])
        for i in range(1, len(derivative)):
            phi = pr.D(phi, self.independent_variables[derivative[i]])
        
        u_dd = u
        for i in derivative:
            u_dd = pr.D(u_dd, self.independent_variables[i])

        for i, x in enumerate(self.independent_variables):
            phi+=self.xi[i]*pr.D(u_dd, x)

        # Clean the expression
        phi = sym.simplify(pr.clean(phi))
        return phi
        
        

    def _init_vector_field_coefficients(self) -> tuple:
        """
        Initialize the xi and phi functions.
        """
        xi = []
        phi = []
        for x in self.independent_variables:
            xi.append(sym.Function(f'xi^{x.name}')(*self.independent_variables, *self.dependent_variables))
        for u in self.dependent_variables:
            phi.append(sym.Function(f'phi^{u.name}')(*self.independent_variables, *self.dependent_variables))
        return xi, phi
    
    def output_to_latex(self, monomials: dict) -> str:
        """
        Output monomial dictionary to latex table.
        """

        # Create the table header
        output = "\\begin{table}[] \n \centering \n \\begin{tabular}{c|c} \n monomial & coefficient\\\\ \n \hline \n"

        # Loop through the monomials and add them to the table
        for monomial, coefficient in monomials.items():
            coefficient = self._prepare_coefficient(coefficient)
            if monomial == '':
                output += f'$1$ & ${coefficient}$ \\\\ \n'
            else:
                monomial = self._prepare_monomial(monomial)
                output += f'${monomial}$ & ${coefficient}$ \\\\ \n'
        
        # Close the table
        output += "\end{tabular}\n \caption{Caption}\n \label{tab:label}\n \end{table}"
        return output

    def _prepare_monomial(self, monomial: str) -> str:
        """
        Prepare the monomial for latex output.
        """
        monomial_list = monomial.split("_")
        out_monomial = f"{monomial_list[0]}_"+"{"
        for m in monomial_list[1:]:
            for t in m:
                if self._is_dependent(t):
                    out_monomial += "}"+f"{t}_"+"{"
                elif t == '^':
                    out_monomial += "}^{"
                else:
                    out_monomial += f"{t}"

        out_monomial += "}"
        return out_monomial

    def _is_dependent(self, t: str) -> bool:
        """
        Check if the given variable is dependent.
        """
        for var in self.dependent_variables:
            if t == var.name:
                return True
        return False
    
    def _prepare_coefficient(self, coefficient) -> str:
        """
        Prepare the coefficient for latex output.

        Currently not properly implemented. WIP.
        """
        # Create the out coefficient string
        out_coefficient = ""

        # Most likely the coefficient is an addition of terms
        if isinstance(coefficient, sym.Add):
            terms = coefficient.args
            # Split the terms and loop through them to prepare them for latex output case by case
            for term in terms:
                # If the term is a multiplication, split the factors and loop through them
                if isinstance(term, sym.Mul):
                    factors = term.args
                    for factor in factors:
                        # If it is an integer, just add it to the coefficient as a string
                        if isinstance(factor, sym.Integer):
                            out_coefficient += f"{factor}"
                        # If it is another function, prepare it for latex output and add as a string
                        else:
                            factor_name = self._prepare_function(factor)
                            out_coefficient += factor_name
#                        else:
#                           raise NotImplementedError(f'Factor type {type(factor)} not implemented')
                # If the term is a function, prepare it for latex output and add as a string
                else:
                    factor_name = self._prepare_function(factor)
                    out_coefficient += factor_name
                #else:
                #    raise NotImplementedError(f'Factor type {type(factor)} not implemented')
                out_coefficient += "+"
        # If the coefficient is a multiplication, split the factors and loop through them.
        elif isinstance(coefficient, sym.Mul):
            factors = coefficient.args
            for factor in factors:
                if isinstance(factor, sym.Integer):
                    out_coefficient += f"{factor}"
                else:
                    factor_name = self._prepare_function(factor)
                    out_coefficient += factor_name
                #else:
                #    raise NotImplementedError(f'Factor type {type(factor)} not implemented')
        else:
            factor_name = self._prepare_function(coefficient)
            out_coefficient += factor_name
        #else:
        #    raise NotImplementedError(f'Factor type {type(coefficient)} not implemented')
        out_coefficient = self._clean_double_plus(out_coefficient)
        return out_coefficient        

    def _clean_double_plus(self, string: str) -> str:
        """
        Clean the double plus signs from the string.
        """
        # Replace double plus
        out_string = string.replace("++", "+")
        # Replace plus minus
        out_string = out_string.replace("+-", "-")

        # Remove leading and trailing plus signs
        if out_string[0] == "+":
            out_string = out_string[1:]
        if out_string[-1] == "+":
            out_string = out_string[:-1]

        #replace multiplication by one
        out_string = out_string.replace("1", "")
        return out_string

    def _prepare_function(self, factor) -> str:
        factor_name = factor.name
        factor_name_list = factor_name.split('_')
    
        function = factor_name_list[0].split("^")[0]
        exponent = factor_name_list[0].split("^")[1]
        subscript = factor_name_list[1]

        return "\\"+function+"_"+"{"+subscript+"}^"+exponent