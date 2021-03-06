import numpy as np
import matplotlib.pyplot as plt

def f_(x):
    return x**3 + 6 * (x**2) - 19 * x -84

def df_(x):
    return 3 * (x**2) + 12 * x - 19

def metoda_newton_raphson(f, df, x_old, epsilon):

    x_new = x_old - f(x_old)/df(x_old)  # Prima aproximare
    counter = 1  # Numar de pasi
    while np.abs(f(x_new)) > epsilon:  # Primul criteriu de oprire
        x_old = x_new
        x_new = x_old - f(x_old) / df(x_old)
        counter += 1

    return x_new, counter

def subs_asc_fast(a, b):

    """ Initializarea vectorului solutiei numerice. """
    n = b.shape[0] - 1
    x_num = np.zeros(shape=n + 1)

    """ Determinarea solutiei numerice. """
    x_num[0] = b[0] / a[0, 0]  # Scrie ultima componenta a solutiei
    for k in range(1, n+1):
        s = np.dot(a[k, :k], x_num[:k])
        x_num[k] = (b[k] - s) / a[k, k]

    return x_num

def DescQR(A,b):
    N = A.shape[0]
    Q = np.zeros([N,N])
    R = np.zeros([N,N])
    for i in range (N):
        R[0,0] += A[i,0] ** 2
    R[0,0] = np.sqrt(R[0,0])
    for i in range (N):
        Q[i,0] = A[i,0] / R[0,0]
    for i in range(1,N):
        for j in range(N):
            R[0,i] += Q[j,0] * A[j,i]

    for k in range(1,N):
        s1 = 0
        for j in range (N):
            s1 += A[j,k] ** 2

        s2 = 0
        for j in range(i - 1):
            s2 += R[j,k] ** 2

        R[k,k] = np.sqrt(s1 - s2)

        for i in range (N):
            s3 = 0;
            for j in range(k - 1):
                s3 += R[j,k] * Q[i,j]
            Q[i,k] = (A[i,k] - s3) / R[k,k]

        for j in range(k + 1, N):
            for i in range(N):
                R[k,j] += Q[i,k] * A[i,j]

    x = subs_asc_fast(R, np.matmul(Q.T, b))
    return Q, R, x

def ex1():
    A = np.array([
        [1., 1., 0.],
        [1., 0., 1.],
        [0., 1., 1.]
    ])
    b = np.array([1., 2., 5.])
    (Q, R, x) = DescQR(A,b)
    #print(Q)
    #print(np.matmul(Q,Q.T))
    #print(R)
    print("Solutia sistemului este:")
    print(x)
    ok = False
    y = x == np.linalg.solve(A,b)
    for i in range(3):
        if y[i] == True:
            ok = True
    print(ok)

def ex2():
    # Am ales subintervalele [-8, -5], [-5, 0], [0, 5] pentru ca cele trei solutii ale ecuatiei sunt (-7, -3, 4)
    # si astfel fiecare se afla intr-unul dintre cele 3 subintervale
    A = [-8,-5,0]
    B = [-5,0,5]
    X0 = [-7.5, -2.8, 4.5]
    # Am ales aceste puncte de start pentru ca rspecta conditia f(x0)*f''(x0) > 0
    # care asigura ca sirurile aproximarilor raman in subintervalele respective
    # -25.85 * (-10.5) > 0 pentu x0 = -7.5
    # -5.712 * (-4.8) > 0 pentu x0 = -2.8
    # 43.125 * 15 > 0 pentu x0 = 4.5

    eps = 1e-3
    x_ = np.linspace(-8, 5, 50)
    y_ = f_(x_)
    plt.figure(0)  # Initializare figura
    plt.plot(x_, y_, linestyle='-', linewidth=3)  # Plotarea pentru  functie
    plt.legend(['f(x)', 'x_num'])  # Adauga legenda
    print("Solutiile ecuatiei sunt:")
    for i in range(3):
        a = A[i]
        b = B[i]
        x0 = X0[i]
        x, N = metoda_newton_raphson(f=f_, df=df_, x_old=x0, epsilon=eps)
        plt.scatter(x, 0, s=50, c='black', marker='o')
        print(x)
    plt.axvline(0, c='black')  # Adauga axa OY
    plt.axhline(0, c='black')  # Adauga axa OX
    plt.xlabel('x')  # Label pentru axa OX
    plt.ylabel('f(x)')  # Label pentru axa OY
    plt.title('Metoda Newton_Raphson')  # Titlul figurii
    plt.grid(b=True)
    plt.show()

#ex1()
ex2()