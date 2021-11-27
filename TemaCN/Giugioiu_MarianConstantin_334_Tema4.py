import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
eps = 1e-5 #epsilon pentru ex1

def f2_(x): #functia pentru ex2
    return (1 / (1.9 * np.sqrt(2.0 * np.pi))) * np.exp(-(x ** 2) / (2 * (1.9 ** 2)))

def expresie(f):
    x = sp.symbols('x')
    diff2 = f.diff().diff() #expresia asociata derivatei de ordin 2 a functiei
    f = sp.lambdify(x, f)
    diff2 = sp.lambdify(x, diff2)
    #folosesc lambdify pentru a putea evalua expresiile functiei si derivatei
    return f, diff2

def diferente_finite(X, Y):
    n = X.shape[0] #lungimea lui X
    diff2 = np.zeros(n) #derivata de ordin 2 aproximativa
    h = X[2] - X[1] #diferenta dintre 2 diviziuni
    # valoarea fiecarui element din X prin formula descrisa in ex1
    for i in range(1, n - 1):
        diff2[i] = (Y[i+1] - 2 * Y[i] + Y[i - 1]) / h ** 2
    return diff2

def intergrare(f, x, metoda):
    n = x.shape[0] #lungimea lui x
    sol = 0 #valoarea integralei aproximative
    #am folosit formulele date la curs
    #in loc sa scriu n ca 2 * m si sa fac un for de la 0 la m-1 folosind 2 * i in formula
    #am facut for din 2 in 2  de la 0 la n-2 folosind doar i in formula si
    #modificand numarul care se aduna corespunzator ca sa pastrez pozitiile
    if metoda == "dreptunghi":
        for i in range(0, n - 2, 2):
            sol += f(x[i + 1]) * (x[i + 2] - x[i])

    elif metoda == "trapez":
        for i in range(n - 1):
            sol += ((f(x[i]) + f(x[i + 1])) / 2) * (x[i + 1] - x[i])

    elif metoda == "simpson":
        for i in range(0, n - 2, 2):
            sol += ((f(x[i]) + 4 * f(x[i + 1]) + f(x[i + 2])) / 3) * ((x[i + 2] - x[i]) / 2)

    return sol


def ex1():
    #capetele intervalului
    a = -np.pi / 2
    b = np.pi
    N = 2 #nr diviziuni
    x = sp.symbols('x')
    f = sp.cos(0.6 * x) #expresia asociata functiei date
    f, diff2 = expresie(f)

    while True:
        X = np.linspace(a, b, N)
        h = X[1] - X[0]  # diferenta dintre 2 diviziuni(mereu aceeasi)
        #plecand de la formula din curs f''(x) = (f(x+h)-2f(x)+f(x-h))/h**2
        #si tinand cont ca X[i+1] = X[i] + h
        #putem obtine f''(X[i]) =  (f(X[i+1])-2f(X[i])+f(X[i-1]))/h**2
        #pentru a folosi aceasta formula mai trebuie adaugate la inceput si
        #sfarsit pentru X 2 valori cu diferenta h fata de vecinii lor
        #pentru acest lucru am creat X_aux dupa cum urmeaza pentru a inlocui X
        X_aux = np.zeros(N + 2)
        X_aux[1:-1] = X
        X_aux[0] = X[0] - h
        X_aux[-1] = X[-1] + h
        X = X_aux

        Y_f = f(X) #imaginea lui X prin functia data
        Y_diff2_aproximativ = diferente_finite(X, Y_f) #Imaginea lui X prin functia
        #data de formula descrisa mai sus
        Y_diff2 = diff2(X) #imaginea lui x prin derivata de ordin 2 exacta

        max = np.max(abs(Y_diff2[1:-1] - Y_diff2_aproximativ[1:-1]))
        #diferenta maxima dintre valorile derivatei exacte si celei aproximative
        #pentru fiecare valoare a lui x, va reprezenta eroarea de trunchiere

        if(max <= eps): #conditia de oprire
            print("Numarul de puncte al discretizarii intervalului :", N)
            plt.figure(0)  # Initializare figura
            #plotarea derivatei exacte
            plt.plot(X[1: -1], Y_diff2[1: -1], c='k', linewidth=2, label='derivata exacta')
            #plotarea derivatei aproximative cu linie punctata portocalie
            plt.plot(X[1: -1], Y_diff2_aproximativ[1: -1], c='orange', linewidth=2, linestyle='--',label='derivata aproximativa')
            plt.xlabel('x')  # Label pentru axa OX
            plt.ylabel("f''(x)")  # Label pentru axa OY
            plt.title('Derivata de ordin 2')
            plt.legend()
            plt.show()
            break
        N += 1


def ex2():
    x = np.linspace(-19, 19, 100)
    print("Valoare aproximativa a integralei cu formula de cuadratura sumata a dreptunghiului: ",intergrare(f2_, x, "dreptunghi"))
    print("Valoare aproximativa a integralei cu formula de cuadratura sumata a trapezului: ",intergrare(f2_, x, "trapez"))
    print("Valoare aproximativa a integralei cu formula de cuadratura sumata Simpson: ",intergrare(f2_, x, "simpson"))

ex1()
ex2()