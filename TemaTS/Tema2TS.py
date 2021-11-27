import numpy as np
import math as math
import matplotlib.pyplot as plt
#raportul dintre functia de repartitie a Gama(0,1,v) si functia de repartitie pentru Exponentiala(v)
def r_(x,v):
    return (v * (x**(v-1)) * np.exp(-x)) / (f_gama(v) * np.exp(-x/v))

def f_gama (v):
    return ((v - 1) ** (v - 1)) * (math.exp (1 - v)) * math.sqrt (2 * math.pi * (v - 1))

def media (x):
    return  np.sum(x) / len(x)


def dispersia (x):
    m = media(x)
    return np.sum(x**2) / len(x) - m ** 2


def var_gama(A,l,v):

    a = v**v * np.exp(1 - v) / f_gama(v) #punctul de maxim al functiei r_

    while(True):
        u = np.random.uniform(0,1)
        u_aux = np.random.uniform(0, 1)
        e = -np.log(u_aux) * v
        #e = np.random.exponential(v)

        r = r_(e,v) / a
        if u <= r:
            return e/l + A

def test_gama(A,l,v):
    nr_samples = 10000
    g_test = np.random.gamma(v, 1 / l, nr_samples)
    g = []
    for i in range(nr_samples):
        g.append(var_gama(A,l,v))

    E_teoretic = A + v/l
    Var_teoretic = v/(l*l)

    E_empiric = media(np.array(g))
    Var_empiric = dispersia(np.array(g))

    E_test = media(np.array(g_test))
    Var_test = dispersia(np.array(g_test))

    print(E_teoretic, Var_teoretic)
    print(E_empiric,Var_empiric)
    print(E_test,Var_test)

    count, bins, ignored = plt.hist(g, 30, density=True)
    plt.plot(bins, np.ones_like(bins), linewidth=2, color='r')
    plt.show()

    count, bins, ignored = plt.hist(g_test, 30, density=True)
    plt.plot(bins, np.ones_like(bins), linewidth=2, color='r')
    plt.show()

def var_hipergeometrica(N,p,n):
    A = math.ceil(N * p) #bile albe
    X = 0 #numar bile albe extrase
    j = 0 #numar incercari
    while(j < n):
        u = np.random.uniform(0, 1) #extragem o bila
        if u < p: #am extras o bila alba
            X += 1
            S = 1
        else: #am extras o bila neagra
            S = 0
        N -= 1
        A = A - S
        p = A / N
        j += 1
    return X

def test_hipergeometrica(N,p,n):
    A = math.ceil(N * p)
    B = N - A
    nr_samples = 10000
    h_test = np.random.hypergeometric(A,B,n,nr_samples)

    h = []
    for i in range(nr_samples):
        h.append(var_hipergeometrica(N,p,n))

    E_teoretic = n*p
    Var_teoretic = n*p*(1-p)*(N-n)/(N-1)

    E_empiric = media(np.array(h))
    Var_empiric = dispersia(np.array(h))

    E_test = media(np.array(h_test))
    Var_test = dispersia(np.array(h_test))

    print(E_teoretic, Var_teoretic)
    print(E_empiric, Var_empiric)
    print(E_test,Var_test)


    count, bins, ignored = plt.hist(h_test, 30, density=True)
    plt.plot(bins, np.ones_like(bins), linewidth=2, color='r')
    plt.show()
    count, bins, ignored = plt.hist(h, 30, density=True)
    plt.plot(bins, np.ones_like(bins), linewidth=2, color='r')
    plt.show()

def ex1():
    print("Gama:")
    A = 0
    l = 4
    v = 6
    #print(var_gama(A,l,v))
    test_gama(A,l,v)
def ex2():
    print("Hipergeometrica:")
    N = 40
    A = 20
    p = A/N
    n = 15
    #print(var_hipergeometrica(N,p,n))
    test_hipergeometrica(N,p,n)

ex1()
ex2()
