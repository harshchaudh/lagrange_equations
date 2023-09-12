import sympy as sp

def input_value(value):
    correct = "False"
    inputs = ["True", "False"]
    while correct == "False":
        print(f'\nIs this correct: {value}? ')
        correct = input("Type 'True' or 'False': ")
        while correct not in inputs:
            correct = input("Type 'True' or 'False': ")
        if correct == "False":
            value = input("Retype the correct value: ")
    return value

def generalised_coordinates():
    no_generalised_coordinates = int(input("How many generalised co-ordinates exist for this mechanism: "))
    no_generalised_coordinates = input_value(no_generalised_coordinates)
    count = 1

    generalised_coordinates_list = []

    while count <= no_generalised_coordinates:
        temp = input("Enter generalised co-ordinate: ")
        generalised_coordinate = sp.Function(temp)(t)
        generalised_coordinates_list.append(generalised_coordinate)
        count += 1

    return generalised_coordinates_list

def sympified(equation):
    try:
        sp.sympify(equation)
    except:
        print("\nERROR - Try again")
        equation = input_value(equation)
        sympified(equation)

    return sp.sympify(equation)
    
print("\nEnter time derivates as 'j(t).diff(t, n); \nWhere, \n 'j' is the general co-ordinate entered, \n 'n' is the n-th time derivative.\n")
print("E.g. 'x(t).diff(t, 2) is the second derivate, namely the acceleration, of x with the respect to time.\n")

m, g, l, t = sp.symbols("m g l t")

# General co-ordinates
generalised_coordinates_list = generalised_coordinates()

# Kinetic Energy
T_KE = input("\nKinetic energy equation for the motion of the mechanism: ")
T_KE_str = input_value(T_KE)
T_KE_sp = sympified(T_KE_str)

# Potential Energy
V_PE = input("\nPotential energy equation for the motion of the mechanism: ")
V_PE_str = input_value(V_PE)
V_PE_sp = sympified(V_PE_str)

L = sp.simplify(T_KE_sp - V_PE_sp)

langragian_list = []

for generalised_coordinate in generalised_coordinates_list:
    term_1 = sp.diff(L, generalised_coordinate.diff(t))
    term_2 = sp.diff(term_1, t)

    term_3 = sp.diff(L, generalised_coordinate)

    langragian = sp.simplify(term_2 - term_3)
    langragian_list.append(langragian)

print(langragian_list)
    
