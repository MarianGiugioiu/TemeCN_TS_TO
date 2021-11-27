import numpy as np
import matplotlib.pyplot as plt

#METODE
def metoda_bisectiei(f, a_old, b_old, eps):
    assert a_old < b_old, 'a trebuie sa fie mai mic strict ca b!'
    x_num = (a_old + b_old) / 2  # Prima aproximare
    N = int(np.floor(np.log2((b_old - a_old) / eps)))  # Criteriul de oprire
        # int pentru ca `range` de mai jos nu accepta float

    for _ in range(1, N):
        if f(x_num) == 0:
            break
        elif np.sign(f(a_old)) * np.sign(f(x_num)) < 0:
            b_old = x_num
        else:
            a_old = x_num

        x_num = (a_old + b_old) / 2

    return x_num, N

def metoda_newton_raphson(f, df, x_old, epsilon):

    x_new = x_old - f(x_old)/df(x_old)  # Prima aproximare
    counter = 1  # Numar de pasi
    while np.abs(f(x_new)) > epsilon:  # Primul criteriu de oprire
        x_old = x_new
        x_new = x_old - f(x_old) / df(x_old)
        counter += 1

    return x_new, counter

def pozitie_falsa(f, a_old, b_old, eps):
    counter = 0  #Numar de pasi
    x_old = (a_old * f(b_old) - b_old * f(a_old)) / (f(b_old) - f(a_old))  #Prima aproximare

    while True:
        counter += 1
        if(f(x_old) == 0):
            x_new = x_old
            break
        elif f(a_old) * f(x_old) < 0:
            #schimbarea capetelor intervalului in care caut solutia
            a_new  = a_old
            b_new = x_old
            #calcularea noii pozitii
            x_new = (a_new * f(b_new) - b_new * f(a_new)) / (f(b_new) - f(a_new))
        elif f(a_old) * f(x_old) > 0:
            #schimbarea capetelor intervalului in care caut solutia
            a_new = x_old
            b_new = b_old
            # calcularea noii pozitii
            x_new = (a_new * f(b_new) - b_new * f(a_new)) / (f(b_new) - f(a_new))
        if abs(x_new - x_old) / abs(x_old) < eps:  #Criteriul de oprire
            break

        x_old = x_new #actualizarea pozitiei precedente

    return x_new, counter

def secanta(f, a_old, b_old, x0, x1, eps):
    counter = 1  #Numar de pasi

    while(abs(x1 - x0) / abs(x0) >= eps):  #Criteriul de oprire
        counter += 1
        x2 = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0))  #calcularea noii pozitii
        #actualizare pentru cele 2 pozitii precedente
        x0 = x1
        x1 = x2

    return x2, counter

#FUNCTII
def f1_(x): #functia descrisa in exercitiul 1
    return x**2 - 13

def f2a_(x): #prima functie din egalitatea de la exercitiul 2
    return np.exp(x-2)

def f2b_(x): #a doua functie din egalitatea de la exercitiul 2
    return np.cos(np.exp((x-2))) + 1

def f2_(x): #functie folosita pentru a egala cu 0 diferenta celor 2 functii de mai sus
    return f2a_(x) - f2b_(x)

def df2_(x): #derivata functie precedente
    return np.exp(x-2) + np.exp(x-2) * np.sin(np.exp((x-2)))

def f3_(x): #functia de la exercitiul 3
    return x**3 + 3*(x**2) - 18*x - 40

def f4_(x): #functia de la exercitiul 4
    return x**3 + x**2 - 4*x - 4

#EXERCITII
def ex1():
    print("Exercitiul 1")
    #Am transformat ecuatia sqrt(13) = x in ecuatia x**2 - 13 = 0
    #Am ales intervalul [a,b] pentru ca solutia se afla in acest interval
    a = 0
    b = 5
    eps = 1e-7 #pentru o precizie de 7 zecimale

    x, N = metoda_bisectiei(f=f1_, a_old=a, b_old=b, eps=eps)

    print('Metoda Bisectiei')
    print('Intervalul: [{:.5f}, {:.5f}]'.format(a, b))
    print('Solutia numerica: x_num = {:.7f}'.format(x))
    print()

def ex2():
    print("Exercitiul 2")
    #Am ales intervalul [a,b] pentru ca solutia se afla in acest interval
    a = 0
    b = 3
    x0 = 2.5 #punctul de pornire pentru metoda cu proprietatea ca f(x0)*f''(x0) > 0
    eps = 1e-5

    x_ = np.linspace(a, b, 50)  # Discretizare a intervalului [a, b]
    y1_ = f2a_(x_) #valorile pentru prima functie in intervalul dat
    y2_ = f2b_(x_) #valorile pentru a doua functie in intervalul dat

    x, N = metoda_newton_raphson(f=f2_, df=df2_, x_old=x0, epsilon=eps) #rezultatul ecuatiei f1(x) - f2(x) = 0 este valoarea in care cele doua functii se intersecteaza
    y = f2a_(x) #valoarea primei functii pentru punctul de interscetie

    print('Metoda Newton_Raphson')
    print('Intervalul: [{:.5f}, {:.5f}]'.format(a, b))
    print('Solutia numerica: x_num = {:.5f}'.format(x))
    print('Valoarea functiilor in punctul de intersectie = {:.5f}'.format(y))
    print('Numarul de iteratii: N = {}'.format(N))
    print()

    plt.figure(0)  # Initializare figura
    plt.plot(x_, y1_, linestyle='-', linewidth=3)  #Plotarea pentru prima functie
    plt.plot(x_, y2_, linestyle='-', linewidth=3) #Plotarea pentru a doua functie
    plt.scatter(x, y, s=50, c='black', marker='o')
    plt.legend(['f(x)', 'x_num'])  # Adauga legenda
    plt.axvline(0, c='black')  # Adauga axa OY
    plt.axhline(0, c='black')  # Adauga axa OX
    plt.xlabel('x')  # Label pentru axa OX
    plt.ylabel('f(x)')  # Label pentru axa OY
    plt.title('Metoda Newton_Raphson')  # Titlul figurii
    plt.axis('scaled')
    plt.grid(b=True)
    plt.show()

def ex3():
    print("Exercitiul 3")
    A = -5
    B = 5
    eps = 1e-5

    x_ = np.linspace(A, B, 50)  # Discretizare a intervalului [A, B]
    y_ = f3_(x_) #valorile pentru functie in intervalul dat

    plt.figure(0)  # Initializare figura
    plt.plot(x_, y_, linestyle='-', linewidth=3)  # Plotarea functiei

    #punctele de separare in subintervale - aproximare pentru punctele in care derivata se anuleaza
    a1 = (-2 - np.sqrt(28)) / 2
    b1 = (-2 + np.sqrt(28)) / 2
    a = [-5., a1, b1]  # Capetele din stanga ale intervalelor
    b = [a1, b1, 5.]  # Capetele din dreapta ale intervalelor

    for i in range(len(a)):
        # Calculeaza solutia numerica si numarul de iteratii
        x_num, N = pozitie_falsa(f=f3_, a_old=a[i], b_old=b[i], eps=eps)

        # Printeaza la consola rezultatele
        print('Metoda Pozitiei false')
        print('Intervalul: [{:.5f}, {:.5f}]'.format(a[i], b[i]))
        print('Solutia numerica: x_num = {:.5f}'.format(x_num))
        print('Numarul de iteratii: N = {}'.format(N))
        print('-' * 72)

        plt.scatter(x_num, 0, s=50, c='black', marker='o')  # Adauga in grafic solutia numerica

    plt.legend(['f(x)', 'x_num'])  # Adauga legenda
    plt.axvline(0, c='black')  # Adauga axa OY
    plt.axhline(0, c='black')  # Adauga axa OX
    plt.xlabel('x')  # Label pentru axa OX
    plt.ylabel('f(x)')  # Label pentru axa OY
    plt.title('Metoda Pozitiei false')  # Titlul figurii
    plt.grid(b=True)
    plt.show()

def ex4():
    print("Exercitiul 4")
    A = -3
    B = 3
    eps = 1e-5

    x_ = np.linspace(A, B, 50)  # Discretizare a intervalului [A, B]
    y_ = f4_(x_) #valorile pentru functie in intervalul dat

    plt.figure(0)  # Initializare figura
    plt.plot(x_, y_, linestyle='-', linewidth=3)  # Plotarea functiei

    #punctele de separare in subintervale - aproximare pentru punctele in care derivata se anuleaza
    a1 = (-1 - np.sqrt(13)) / 3
    b1 = (-1 + np.sqrt(13)) / 3
    a = [-3., a1, b1]  # Capetele din stanga ale intervalelor
    b = [a1, b1, 3.]  # Capetele din dreapta ale intervalelor

    for i in range(len(a)):
        # Calculeaza solutia numerica si numarul de iteratii
        d = b[i] - a[i] #distanta dintre capetele subintervalului
        #am ales punctele de start pentru metoda secantei ca fiind la o cincime distanta de capetele subintervalului
        x_num, N = secanta(f=f4_, a_old=a[i], b_old=b[i], x0=a[i]+d/5, x1=b[i]-d/5, eps=eps)

        # Printeaza la consola rezultatele
        print('Metoda Secantei')
        print('Intervalul: [{:.5f}, {:.5f}]'.format(a[i], b[i]))
        print('Solutia numerica: x_num = {:.5f}'.format(x_num))
        print('Numarul de iteratii: N = {}'.format(N))
        print('-' * 72)

        plt.scatter(x_num, 0, s=50, c='black', marker='o')  # Adauga in grafic solutia numerica

    plt.legend(['f(x)', 'x_num'])  # Adauga legenda
    plt.axvline(0, c='black')  # Adauga axa OY
    plt.axhline(0, c='black')  # Adauga axa OX
    plt.xlabel('x')  # Label pentru axa OX
    plt.ylabel('f(x)')  # Label pentru axa OY
    plt.title('Metoda Secantei')  # Titlul figurii
    plt.grid(b=True)
    plt.show()

ex1()
ex2()
ex3()
ex4()


