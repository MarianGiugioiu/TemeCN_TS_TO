import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(precision=2)

def canonical(list,type):
    if type == "max":
        return all(c[-1] == "<=" for c in list)
    else :
        return all(c[-1] == ">=" for c in list)

def equality(list):
    if list[-1] == "<=":
        return list.insert(-2,1) #adauga variabila cu coeficient 1
    else:
        return list.insert(-2,-1) #adauga variabila cu coeficient -1

def standardize(constraints, target):
    #transforma inegalitatile in egalitati adaugand variabile sistemului
    nr_free_coef = 0 #numar variabile adaugate
    nr_coef = len(constraints[0]) - 2 #numar variabile initiale
    for c in constraints:
        if c[-1] != "==":
            for d in constraints:
                if d != c:
                    d.insert(-2,0) #pentru toate inegalitatile diferite de cea curenta adaug variabila cu coeficient 0
            equality(c)
            nr_free_coef += 1

    #actualizez coeficientii functiei pe care o maximizez/minimizez
    coefs = [-x for x in target[:nr_coef]]
    free_coefs = [0] * nr_free_coef
    target = coefs + free_coefs + [1] + [0] + ['==']
    #adaug inca o variabila cu coeficient 0 la fiecare ecuatie
    for c in constraints:
        c.insert(-2,0)
    nr_free_coef += 1
    return (constraints,target,nr_free_coef)

def to_matrix(constraints, target):
    #matricea existinsa pentru sistemul creat
    mat = []
    for c in constraints:
        mat.append(np.hstack(c[:-1]))
    mat.append(target[:-1])
    return np.vstack(mat)

def pivot(A, px, py):
    #genereaza un pas din metoda Gauss cu pivotul pe pozitia (px,py)
    p = A[px,py]
    A[px] /= p
    for i in range(len(A)):
        if i != px:
            A[i] -= A[i,py] * A[px]

def primal_simplex(constraints,target):
    assert canonical(constraints, "max"), 'Linear program must be in canonical form!'
    nr_coef = len(target)
    X = [0] * nr_coef
    S = [-1] * len(constraints)
    (constraints, target, nr_free_coef) = standardize(constraints, target)
    A = to_matrix(constraints, target)
    while True:
        print(A)
        print()
        py = -1
        min_c = 0
        #caut cea mai mica valoare negativa din egalitatea creata din functia initiala
        for column, value in enumerate(A[-1]):
            if value < min_c:
                py = column
                min_c = value

        #am gasit optimul
        if py == -1:
            break
        px = -1
        min_r = np.inf
        #aleg elementul de pe coloana corespunzatoare minimului gasit mai sus care minimizeaza raportul dintre termenul drept al egalitatii de pe linia respectiva si elementul precizat
        for row, value in enumerate(A[:-1, -1]):
            elem = A[row, py]
            if elem > 0:
                r = value / elem
                if r < min_r:
                    px = row
                    min_r = r
        #optimul este infinit
        if px == -1:
            print('optim = inf')
            return

        pivot(A, px, py)
        S[px] = py

    #generez solutia
    for i in range(len(constraints)):
        X[S[i]] = A[i,-1]
    for i in range(nr_coef):
        print(f'x{i + 1} = {X[i]}')
    print(f'optim = {A[-1][-1]}')

def dual_simplex(constraints,target):
    assert canonical(constraints, "min"), 'Linear program must be in canonical form!'
    nr_coef = len(target)
    target.append(1.)
    target.append(('=='))
    B = to_matrix(constraints,target)
    C = B.transpose().tolist()
    for c in C:
        c.append('<=')
    del C[-1][-1]
    del C[-1][-1]

    #am construit matricea corespunzatoare programului liniar si am transpus-o
    (constraints, target, nr_free_coef) = standardize(C[:-1],C[-1])

    A = to_matrix(constraints, target)
    #se rezolva analog ca la cealalta metoda
    while True:
        py = -1
        min_c = 0
        for column, value in enumerate(A[-1]):
            if value < min_c:
                py = column
                min_c = value

        if py == -1:
            break
        px = -1
        min_r = np.inf
        for row, value in enumerate(A[:-1, -1]):
            elem = A[row, py]
            if elem > 0:
                r = value / elem
                if r < min_r:
                    px = row
                    min_r = r

        if px == -1:
            print('optim = inf')
            return

        pivot(A, px, py)

    for i in range(nr_coef):
        print(f'x{i + 1} = {A[-1,nr_coef+i]}')
    print(f'optim = {A[-1][-1]}')


def show(constraints,target,type):
    #face afisarea programului liniar
    string = type + " "

    string += str(target[0]) + "??x1"
    for i in range(1,len(target)):
        string += " + " + str(target[i]) + "??x" + str(i+1)
    string += "\n"

    for j in range(len(constraints)):
        string += str(constraints[j][0]) + "??x1"
        for i in range(1, len(target)):
            string += " + " + str(constraints[j][i]) + "??x" + str(i + 1)
        string += " " + str(constraints[j][-1]) + " " + str(constraints[j][-2])
        string += "\n"

    return string

def ex1():
    print("Simplex primal")
    '''constraints = [
        [2., 3., 2., 1000., '<='],
        [1., 1., 2., 800., '<=']
    ]
    type = "max"

    target = np.array([7., 8., 10.])'''
    constraints = [
        [1., -1., 2., 4., '<='],
        [2., 2., -1., 1., '<='],
        [2., 0., 4., 4., '<=']
    ]
    type = "max"

    target = np.array([2., 3., 1.])

    print("Programul liniar:")
    print(show(constraints, target, type))

    print("Rezultat:")
    primal_simplex(constraints, target)


def ex2():
    print("Simplex dual")
    constraints = [
        [3., 5., 2., 60., '>='],
        [4., 5., 4., 72., '>='],
        [2., 4., 5., 100., '>='],
    ]
    type = "min"

    target = [5., 10., 8.]
    print("Programul liniar:")
    print(show(constraints, target, type))

    print("Rezultat:")
    dual_simplex(constraints, target)

ex1()
print()
#ex2()