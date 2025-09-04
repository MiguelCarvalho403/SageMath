"""
This module contains pdsolve() and different helper functions that it
uses. It is heavily inspired by the ode module and hence the basic
infrastructure remains the same.
**Functions in this module**
    These are the user functions in this module:
    - pdsolve()     - Solves PDE's
    - classify_pde() - Classifies PDEs into possible hints for dsolve().
    - pde_separate() - Separate variables in partial differential equation either by
                       additive or multiplicative separation approach.
    These are the helper functions in this module:
    - pde_separate_add() - Helper function for searching additive separable solutions.
    - pde_separate_mul() - Helper function for searching multiplicative
                           separable solutions.
**Currently implemented solver methods**
The following methods are implemented for solving partial differential
equations.  See the docstrings of the various pde_hint() functions for
more information on each (run help(pde)):
  - 1st order linear homogeneous partial differential equations
    with constant coefficients.
  - 1st order linear general partial differential equations
    with constant coefficients.
  - 1st order linear partial differential equations with
    variable coefficients.
"""

#from sage.all import *
from sage.calculus.functional import diff
from sage.misc.functional import integral
from sage.symbolic.relation import solve
from sage.arith.misc import gcd
from sage.repl.rich_output.pretty_print import show
from sage.misc.html import html
from sage.misc.latex import latex
from sage.functions.generalized import sign
from sage.symbolic.ring import SR
from colorama import Fore


def sep_simplify(expression): # Simplify expressions 
    return expression.expand().factor().expand().combine()

def eq_operands(pde): # Simplify the expressions and split them into a list, in case they haven't been already
    if pde.is_relational():
        splitted_pde = list(map(lambda q: q.factor(), sep_simplify(pde.lhs()-pde.rhs()).operands()))
    else:
        splitted_pde = list(map(lambda q: q.factor(), sep_simplify(pde).operands()))
    return splitted_pde

def is_trivially_separable(splitted_pde, sep_variables, show_log):# Verify whether the PDE is separated into terms that contain only one variable, return True or False
    flag=True
    for k in splitted_pde:
        if k and (sep_variables[0] in k.variables() and sep_variables[1] in k.variables()):
            flag = False
    if show_log and flag:
        print(Fore.GREEN + 'is_trivially_separable' + Fore.RESET)
    return flag

def list_of_terms(splitted_pde, sep_variables): # Group the terms that contain the same variable
    return [sum([q for q in splitted_pde if v in q.variables()]) for v in sep_variables] # this function is ignoring constant independent terms

def pre_separate(pde, sep_variables, check, show_log): # Decompose the PDE symbolic expression into mixed-term and single-variable elements
    splitted_pde = pde.operands() # separate PDE into a list
    separated = is_trivially_separable(splitted_pde, sep_variables, show_log=False) # verify whether the PDE is separated into terms that contain only one variable, return True or False
    if separated: # only if splitted_pde is trivially separated
        pre_partially_separated_pde = list_of_terms(splitted_pde, sep_variables=sep_variables) # Group the terms that contain the same variable
    else:
        splitted_pde = eq_operands(pde) # Simplify the expressions and split them into a list, in case they haven't been already
        both_var = [] # store mixed terms
        for q in splitted_pde[:]: #Create a list of terms that contain more than one variable and remove them from splitted_pde 
            if sep_variables[0] in q.variables() and sep_variables[1] in q.variables(): # Check whether a given term in splitted_pde involves more than one variable
                both_var.append(q)
                splitted_pde.remove(q)
            # By the end, each element of the splitted_pde list involves only one variable
        if len(both_var): 
            pre_partially_separated_pde = [sum(both_var)] # group mixed terms in a symbolic expression and store them in a list
            pre_partially_separated_pde += list_of_terms(splitted_pde, sep_variables=sep_variables) # append symbolic expressions that only contain one variable
        else:
            pre_partially_separated_pde = list_of_terms(splitted_pde, sep_variables=sep_variables) # Construct a list where each symbolic expression involves a single variable
            # By the end, pre_partially_separated_pde will be separated into mixed-term elements and single-variable symbolic expression elements
    if check:
        check_result(pde, pre_partially_separated_pde, 'pre_separate test: ', show_log) # flag1
    return pre_partially_separated_pde



def partial_separate(pde, sep_variables, check, show_log): # Decompose the PDE into a list of terms, each involving only one variable
    partially_separated_pde = pre_separate(pde, sep_variables, check, show_log) # Decompose the PDE symbolic expression into mixed-term and single-variable elements
    separated = is_trivially_separable(partially_separated_pde, sep_variables, show_log=False)
    
    if separated: # Verify whether the PDE is separated into terms that contain only one variable
        return partially_separated_pde
    else: # if there are mixed terms in the list of symbolic expressions
        partial_pde = [0]*3
        while not separated  and (partial_pde[1] or partial_pde[2]): #Decompose the mixed-term symbolic expression into mixed-term and single-variable elements until no mixed terms remain
            partial_pde = pre_separate(partially_separated_pde[0], sep_variables, check, show_log)
            partially_separated_pde[0] = partial_pde[0]
            partially_separated_pde[1] += partial_pde[1]
            partially_separated_pde[2] += partial_pde[2]
            separated = is_trivially_separable(partially_separated_pde, sep_variables=sep_variables, show_log=False)
        if check:
            check_result(pde, partially_separated_pde, 'partial_separate test: ', show_log)
        return partially_separated_pde
                     
def list_to_pde(splitted_pde, factor): #path1f Convert a list of equation terms to a unique PDE
    return sum(map(lambda q: (q*factor).full_simplify(),splitted_pde))

def is_sum_separable(partial_pde, sep_variables): # Determine whether the PDE can be separated into a sum of single-variable terms
    if diff(partial_pde,sep_variables[0],sep_variables[1]).canonicalize_radical():
        # if partial_dpe = f(x) + g(y), then d(partial_dpe)/dxdy = d(f(x) + g(y))/dxdy = 0
         return False
    else:
        return True

def sum_separate(partial_pde, sep_variables, check, show_log): # Decompose the mixed terms(fun_var1 and fun_var_2) into a sum of single-variable components and return the sum of them
    #list_of_terms = partial_pde.expand().operands()
    #fun_var_1, fun_var_2 = [sum(filter(lambda term: vars in term.variables(), list_of_terms)) for vars in sep_variables]
    fun_var_1 = integral(diff(partial_pde,sep_variables[1]).canonicalize_radical(),sep_variables[1]) # edp = f(x) + g(y) => integral(diff(edp, y), y) = g(y)
    fun_var_2 = (partial_pde-fun_var_1).canonicalize_radical() # edp = f(x) + g(y) => f(x) = edp -g(y)
    if check:
        check_result(partial_pde, fun_var_1 + fun_var_2, 'sum_separate test: ', show_log)
    return fun_var_1 + fun_var_2 # return f(x) + g(y)
    
def separate_equations(partially_separated_pde, funcs, sep_variables, sep_constant, check, show_log):
    equations = []
    for k, eq in enumerate(partially_separated_pde):
        temp_eq = eq.numerator() + sign(2*k-1)*sep_constant*eq.denominator()
        equations.append(canonicalize_de(temp_eq==0, funcs[k], sep_variables[k], check, show_log)) # flag3
    if check:
        check_pos_separation(equations, sep_variables, show_log) 
        check_mix_diff(equations, sep_variables, show_log)
        check_separated(partially_separated_pde, equations, sep_constant, show_log)  
    return equations

def check_pos_separation(equations, sep_variables, show_log):
    step_message = 'separate_equations test 1: '
    if equations[0].has(sep_variables[1]) or equations[1].has(sep_variables[0]):
        message = Fore.RED + 'failed' + Fore.RESET
        print(step_message + message)
    else:
        if show_log:
            message = Fore.GREEN + 'passed!' + Fore.RESET
            print(step_message + message)
    pass

def check_mix_diff(equations, sep_variables, show_log):
    list_of_eqs = []
    for eq in equations:
        list_of_eqs.append(diff(eq.lhs(),sep_variables[0],sep_variables[1]))
    check_result(0, list_of_eqs, 'check_mix_diff test: ', show_log)
    pass

def check_separated(partially_separated_pde, equations, sep_constant, show_log):

    #print(f"partially_separated_pde: {partially_separated_pde}")
    #print(f"Equações: {equations}")
    #print(f"sep_constant: {sep_constant}")
    initial_eqs = [sign(1-2*k)*solve(equations[k], sep_constant)[0].rhs() for k in range(2)]
    check_result(sum(partially_separated_pde), initial_eqs, 'separate_equations test 2: ', show_log)
    pass

def partial_common_factor(partially_separated_pde, sep_variables, check, show_log):
    #common_factor, expr_list = eliminate_common_factor(eq_operands(partially_separated_pde[0]), return_factor=True, check=False, show_log=False)
    common_factor = get_common_factor(eq_operands(partially_separated_pde[0])) # Get commmon factor from mixed terms
    if common_factor.is_rational_expression():
        index_var_list = [] # list of index of variables that are not part of common factor
        for index, svar in enumerate(sep_variables):
            if svar not in common_factor.variables():
                index_var_list.append(index) 
        if len(index_var_list)==1: # if common factor has only one variable
            not_in_factor_index = index_var_list[0]+1 # contains the index of variable that is not part of the common factor
            if not partially_separated_pde[not_in_factor_index]: # checks for terms that do not contain the common factor  
                # If any term does not share the common factor, dividing the equation by that common factor will result in mixed variables.
                new_partially_separated_pde = [sep_simplify(SR(expr)/common_factor) for expr in partially_separated_pde] # eliminate common factor
                message = Fore.LIGHTBLUE_EX + f'factor {common_factor} eliminated' + Fore.RESET
                #return common_factor, new_partially_separated_pde
            else: # the commom factor is one
                message = Fore.YELLOW + f'common factor {common_factor} ignored' + Fore.RESET # add
                common_factor = SR(1)
                new_partially_separated_pde = partially_separated_pde
        else: # add
            new_partially_separated_pde = partially_separated_pde # add
            message = Fore.YELLOW + f'common factor {common_factor} ignored' + Fore.RESET # add
    else:
        message = Fore.YELLOW + f'common factor {common_factor} ignored' + Fore.RESET
        new_partially_separated_pde = partially_separated_pde
    if show_log: 
        print(message)
    if check:
        check_result(sum(partially_separated_pde), sum(new_partially_separated_pde)*common_factor, 'partial_common_factor test: ', show_log)
    return common_factor, new_partially_separated_pde


def get_common_factor(expr_list): # Get common factor from mixed terms, receive as input a list of mixed terms, returns common factor
    if not expr_list[1] and not expr_list[2]:
        expr_list = sum(expr_list).expand().operands()
    numerators, denominators = zip(*[expr.numerator_denominator() for expr in expr_list if expr])
    return gcd(numerators)/gcd(denominators)
    #return gcd(map(numerator,expr_list)) / gcd(map(denominator,expr_list))

def eliminate_common_factor(initial_expr_list, return_factor=False, check=True, show_log=True):
    common_factor = get_common_factor(initial_expr_list) # get common factor from EDP
    if common_factor.is_rational_expression():
        expr_list = [sep_simplify(SR(expr)/common_factor) for expr in initial_expr_list]
        message = Fore.LIGHTBLUE_EX + f'factor {common_factor} eliminated' + Fore.RESET
    else:
        expr_list = initial_expr_list
        message = Fore.YELLOW + f'common factor {common_factor} ignored' + Fore.RESET
        common_factor = SR(1)
    if show_log:
        print(message)
    if check:
        check_result(sum(initial_expr_list)/common_factor, expr_list, 'eliminate_common_factor test: ', show_log)
    if return_factor:
        return common_factor, expr_list
    else:
        return expr_list

def canonicalize_de(initial_de, fun, variable, check=True, show_log=True, expand=False):
    """[summary]

    Arguments:
        initial_de {[type]} -- [description]
        fun {[type]} -- [description]
        variable {[type]} -- [description]

    Returns:
        [type] -- [description]
    """    
#    from sage.symbolic.expression import is_SymbolicEquation
#    if is_SymbolicEquation(initial_de):
#        initial_de = initial_de.lhs()-initial_de.rhs()
    coefficients = get_de_coefficients(initial_de, fun, variable)
    coeff_2 = 1
    if coefficients[2]:
        coeff_2 = coefficients[2]
        for k in range(3):
            coefficients[k] = (coefficients[k] / coeff_2).canonicalize_radical().expand()
            if not expand: 
                coefficients[k] = coefficients[k].factor()
    new_de = sum([coefficients[k]*diff(fun, variable, k) for k in range(3)]) == 0
    if check:
        check_result(initial_de, new_de*coeff_2, step_message=f'canonicalize_de for {fun}: ', show_log=show_log)
        
        print(f"initial_de: {initial_de}")
        print(f"initial_de: {initial_de}")
        print(f"new_de*coeff_2: {new_de*coeff_2}")
    
    return new_de

def get_de_coefficients(initial_de, fun, variable):
    """[summary]

    Arguments:
        initial_de {[type]} -- [description]
        fun {[type]} -- [description]
        variable {[type]} -- [description]

    Returns:
        [type] -- [description]
    """  
    if initial_de.is_relational(): # verificar se é uma equação
        de = initial_de.lhs()-initial_de.rhs()
    else: de = initial_de  
    return {k:formal_diff(de,diff(fun, variable, k)).expand().factor() for k in range(3)}
    

def formal_diff(f, x): # formal derivative of `f` with respect to the symbolic function `x`
    """The formal derivative of `f` with respect to the symbolic function `x`.

    Arguments:
        f {[type]} -- [description]
        x {[type]} -- [description]

    Returns:
        [type] -- [description]
    """    
    tempX = SR.symbol()
    return f.subs({x: tempX}).diff(tempX).subs({tempX: x})


def change_function(de, old_fun, new_fun, new_fun_factor, variable, expand=False, show_log=True):
    """[summary]

    Arguments:
        de {[type]} -- [description]
        old_fun {[type]} -- [description]
        new_fun {[type]} -- [description]
        new_fun_factor {[type]} -- [description]
        variable {[type]} -- [description]

    Returns:
        [type] -- [description]
    """    
    new_de = de.subs({diff(old_fun, variable, k):diff(new_fun_factor, variable, k) for k in range(3)})
    return canonicalize_de(new_de, new_fun, variable, expand, show_log)

def change_variable(de, fun, old_var, new_var_as_function, new_var_as_var, diff_var_as_function=None, expand=False, show_log=True):
    """[summary]

    Arguments:
        de {[type]} -- [description]
        fun {[type]} -- [description]
        old_var {[type]} -- [description]
        new_var_as_function {[type]} -- [description]
        new_var_as_var {[type]} -- [description]

    Keyword Arguments:
        diff_var_as_function {[type]} -- [description] (default: {None})

    Returns:
        [type] -- [description]
    """    
    new_de = de.subs({diff(fun(old_var),old_var,k):diff(fun(new_var_as_function),old_var,k) for k in range(3)})
    if diff_var_as_function:
        new_de1 = new_de.subs({diff(new_var_as_function, old_var,k+1):diff(diff_var_as_function, old_var,k) for k in range(2)})
    else:
        new_de1 = new_de
    new_de2 = new_de1.subs({new_var_as_function:new_var_as_var})
    return canonicalize_de(new_de2, fun(new_var_as_var), new_var_as_var, expand, show_log)


def check_result(expression_1, expression_2, step_message='are equals? ', show_log=True): # flag2
    """[summary]

    Arguments:
        expression_1 {[type]} -- [description]
        expression_2 {[type]} -- [description]

    Keyword Arguments:
        step_message {str} -- [description] (default: {'are equals? '})
        show_log {bool} -- [description] (default: {True})
    """    
    from collections.abc import Iterable
    from sage.symbolic.expression import Expression
    if isinstance(expression_1, Iterable):
        expression_1 = sum(expression_1)
    if isinstance(expression_2, Iterable):
        expression_2 = sum(expression_2)
    if isinstance(expression_1, Expression) and expression_1.is_relational():
        expression_1 = expression_1.lhs()-expression_1.rhs()
    if isinstance(expression_2, Expression) and expression_2.is_relational():
        expression_2 = expression_2.lhs()-expression_2.rhs() 
    difference = (expression_1-expression_2).canonicalize_radical()
    if difference:
        message = Fore.RED + 'failed' + Fore.RESET
        print(step_message + message)
    else:
        if show_log:
            message = Fore.GREEN + 'passed!' + Fore.RESET
            print(step_message + message)
    pass

def show_equations(equations, from_separation=False):
    """[summary]

    Arguments:
        equations {[type]} -- [description]
    """    
    from sage.manifolds.utilities import ExpressionNice
    if from_separation:
        show(html(r'<br>'))
        show(html(r'$\large \color{green}{\text{Equations separated}}$'))
    show(html(' '))
    for eq in equations:  
        show(html(r'$\text{Parameters:} '+latex(eq.variables())+'$'))
        show(html(' '))
        show(html(r'$\large' + latex(ExpressionNice(eq))+'$'))
        show(html(r'<br>'))
    pass

def _separate_pde(pde,
                 funcs,
                 sep_variables,
                 sep_constant,
                 is_list = False,
                 factor=1,
                 check=True,
                 show_log=True,
                 show_eqs=True,
                 initial_pde=None):
    """Separate variables in partial differential equation either by additive
    or multiplicative separation approach. It tries to rewrite an equation so
    that one of the specified variables occurs on a different side of the
    equation than the others.

    Arguments:
        pde {[type]} -- [description]
        funcs {[type]} -- [description]
        sep_variables {[type]} -- [description]

    Keyword Arguments:
        sep_constant {[type]} -- [description] 
        is_list {bool} -- [description] (default: {False})
        factor {int} -- [description] (default: {1})
        check {bool} -- [description] (default: {True})
        show_eqs {bool} -- [description] (default: {True})
        initial_pde {[type]} -- [description] (default: {None})

    Returns:
        [type] -- [description]

    EXAMPLES:

    ***
    """
    if show_log:
        print(Fore.LIGHTMAGENTA_EX + 'initialized' + Fore.RESET)   
    if is_list:
        pde = list_to_pde(pde, factor) 
    if initial_pde is None:
        initial_pde = pde
    # substitute the original function to a sum
    # check if it is sum peparable, if not, proceed
    # substitute the original function by a product and divide the equation by the product
    # this should be in another function
    partially_separated_pde = partial_separate(pde, sep_variables, check, show_log) # Decompose the PDE into a list of terms, each involving only one variable
    if is_trivially_separable(partially_separated_pde, sep_variables, show_log): # Check whether EDP can be separated
        separated_equations = separate_equations(partially_separated_pde, funcs, sep_variables, sep_constant, check, show_log) # separate equations
        if show_eqs:
            show_equations(separated_equations, from_separation=True)
        return separated_equations 
    else:
        if is_sum_separable(partially_separated_pde[0], sep_variables): # Determine whether the PDE can be separated into a sum of single-variable terms
            partially_separated_pde[0] = sum_separate(partially_separated_pde[0], sep_variables, check, show_log) # separe PDE mix-terms
            return _separate_pde(sum(partially_separated_pde),
                                funcs,
                                sep_variables,
                                sep_constant,
                                check=check,
                                show_log=show_log,
                                show_eqs=show_eqs,
                                initial_pde=initial_pde)
        else:
            common_factor, partially_separated_pde = eliminate_common_factor(partially_separated_pde, return_factor = True, check=check, show_log=show_log)
            if not common_factor==1:
                return _separate_pde(sum(partially_separated_pde), funcs, sep_variables, sep_constant, check=check, show_log=show_log, show_eqs=show_eqs, initial_pde=initial_pde/common_factor)
            else:
                common_factor=1
                common_factor, partially_separated_pde = partial_common_factor(partially_separated_pde, sep_variables, check=check, show_log=show_log)
                if not common_factor==1:
                    return _separate_pde(sum(partially_separated_pde), funcs, sep_variables, sep_constant, check=check, show_log=show_log, show_eqs=show_eqs, initial_pde=initial_pde/common_factor)
                else:
                    return partially_separated_pde


def separate_pde(pde,
                 funcs,
                 sep_variables,
                 sep_constant,
                 is_list = False,
                 factor=1,
                 check=True,
                 show_log=True,
                 show_eqs=True,
                 initial_pde=None):
    separated_equations = _separate_pde(pde,
                                        funcs,
                                        sep_variables,
                                        sep_constant,
                                        is_list,
                                        factor,
                                        check,
                                        show_log,
                                        show_eqs,
                                        initial_pde)
    return separated_equations
#def separate_pde(initial_fun,sep_funcs):
#    pde = split_as_sum()
#    if is_sum_separable(pde, sep_variables):
#            _separate_pde(***)
#            return ***
#    pde = split_as_product_and_divide()
#


