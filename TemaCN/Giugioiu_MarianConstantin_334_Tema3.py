import numpy as np
import matplotlib.pyplot as plt

eps = 10 **(-3) #pentru exercitiul 2 am gasit cel mai bun N pentru acest epsilon
eps1 = 10 ** (-10) #pentru primul exercitiu

def f1_(x, y):
    return 12.5 * (x**2) - 20.0 * x * y - 4 * x + 32.5 * (y**2) - 8 * y

def f2_(x):
    return 3 * np.sin(3 * x) + 3 * np.cos(5 * x) - 1.27 * x

def f3_(x):
    #am rescris functia, scotant afara din paranteze - pentru sin si cos
    return 9 * np.sin(2 * x) - 2 * np.cos(3 * x) + 9.04 * x
def df3_(x):
    #derivata functiei 3
    return 6 * np.sin(3 * x) + 18 * np.cos(2 * x) + 9.04

def pas_descendent(A, b):
    x = np.array([1/3, 1/3])
    #am ales acest punct de inceput pentru ca se afla in apropierea solutiei,
    #apoi am implementat algoritmul prezentat la laborator
    R = b - np.matmul(A, x)

    while(abs(np.prod(R)) > eps1):
        R = b - np.matmul(A, x)
        a = np.matmul(R, R) / np.matmul(np.matmul(R, A), R)
        x = x + a * R

    return x

def gradienti_conjugati(A, b):
    #aceeasi abordare ca la metoda pasului descendent
    x = np.array([1/3, 1/3])
    R_old = b - np.matmul(A, x)
    d = R_old

    while(abs(np.prod(R_old)) > eps1):
        a = np.matmul(R_old, R_old) / np.matmul(np.matmul(d, A), d)
        x = x + a * d
        R_new = R_old - a * np.matmul(A, d)
        B = np.matmul(R_new, R_new) / np.matmul(R_old, R_old)
        d = R_new + B * d
        R_old = R_new

    return x

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

def subs_desc_fast_index(a, b, index):
    """Este aceeasi ca subs_desc_fast, aplicand in plus permutarea data de index la solutie"""
    """ Initializarea vectorului solutiei numerice. """
    n = b.shape[0] - 1
    x_num = np.zeros(shape=n+1)

    """ Determinarea solutiei numerice. """
    x_num[n] = b[n] / a[n, n]  # Scrie ultima componenta a solutiei
    for k in range(n-1, -1, -1):
        s = np.dot(a[k, k+1:], x_num[k+1:])
        x_num[k] = (b[k] - s) / a[k, k]
    x_num_final = np.zeros(shape=n+1)

    for i in range(n+1):
        x_num_final[index[i]] = x_num[i]
    return x_num_final

def meg_pivotare_totala(a, b):
    N = b.ndim
    """In cazul in care dimensiunea lui b este 2(b este matrice) rezolv simultan toate sistemele date de matricea a si fiecare coloana a lui b"""
    if N == 1:
        a_ext = np.concatenate((a, b[:, None]), axis=1)
    else:
        a_ext = np.concatenate((a, b), axis=1)
    n = a.shape[0]
    """Index va retine interschimbarile dintre coloane pe parcursul algoritmului"""
    index = np.arange(0,n,1)

    for k in range(n-1):
        """Pozitia pivotului la pasul k va fi pozitia celui mai mare element din matricea formata din elementele de pe pozitii k-n pe linie si coloana"""
        pos = np.argmax(np.absolute(a_ext[k:n, k:n]))
        p = pos // (n - k)
        m = pos % (n - k)
        p += k
        m += k
        ''' Schimba linia 'k' cu 'p' si coloana 'k' cu 'm' daca pivotul nu se afla pe linia respectiv coloana k, iar pentru coloane modific index '''
        if k != p:
            a_ext[[p, k], :] = a_ext[[k, p], :]
        if k != m:
            a_ext[:, [m, k]] = a_ext[:, [k, m]]
            index[[m,k]] = index[[k,m]]

        """ Zero pe coloana sub pozitia pivotului. """
        for j in range(k+1, n):
            m = a_ext[j, k] / a_ext[k, k]
            a_ext[j, :] -= m * a_ext[k, :]

    """ Gaseste solutia folosind metoda substitutiei descendente si permutarea data de index. """
    if N == 1:
        """Pentru rezolvarea unui singur sistem"""
        x_num = subs_desc_fast_index(a_ext[:, :-1], a_ext[:, -1], index)

    else:
        """Pentru rezolvarea mai multor sisteme simultan"""
        x_num = np.empty([n, n])
        for j in range(n):
            x_num[:, j] = subs_desc_fast_index(a_ext[:, :n], a_ext[:, n+j], index)

    return x_num

def Lagrange_Newton(X, Y, N):
    #este format(asa cum este prezentat in curs) un sistem inferior triunghiular cu coeficientii cautati ca si necunoscute
    M = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,i+1):
            M[i,j] = 1
            for k in range(0,j):
                M[i,j] *= X[i] - X[k]
    x = subs_asc_fast(M,Y)
    return x

def interpolare_LN_elem(C, X, N, x):
    y = 0
    for i in range(0,N):
        s = C[i]
        for j in range(0,i):
            s *= x - X[j]
        y += s
    return y

def spline_cubice(X, Y, N):
    h = X[1] - X[0]
    a = np.zeros(N - 1)
    c = np.zeros(N - 1)
    d = np.zeros(N - 1)
    a[0] = Y[0]
    #a va fi valoarea prin functie a fiecarui element din diviziune
    #pentru a calcula b se formeaza sistemul descris la curs folosind matricile M si P
    M = np.zeros((N,N)) #coeficientii fiecarei ecuatii
    P = np.zeros(N) #termenul din dreapta al ecuatiei
    P[0] = df3_(X[0])
    P[N - 1] = df3_(X[N - 1])
    M[0,0] = 1
    M[N-1,N-1] = 1
    for i in range(1, N-1):
        a[i] = Y[i]
        P[i] = (3 / h) * (f3_(X[i+1]) - f3_(X[i-1]))
        M[i][i-1] = 1
        M[i][i] = 4
        M[i][i+1] = 1

    b = meg_pivotare_totala(M, P)
    #folosesc metoda de eliminare Gauss cu pivotare totala pentru a rezolva sistemul

    #se calculeaza c si d dupa b asa cum este prezentat in curs
    for i in range(0, N-1):
        c[i] = (3/(h**2)) * (f3_(X[i+1]) - f3_(X[i])) - (b[i+1] + 2*b[i])/h
        d[i] = ((-2)/(h**3)) * (f3_(X[i+1]) - f3_(X[i])) - (b[i+1] + 2*b[i])/(h**2)

    return (a,b,c,d)

def interpolare_sc_elem(X, a,b,c,d, N, x):
    #se cauta intervalul din care face parte valoare x si se calculeaza valoare
    #sa prin polinomul corespunzator
    for i in range(0,N-1):
        if((x >= X[i] and x < X[i+1]) or (i == N-2 and x == X[i+1])):
            y = a[i] + b[i] * (x - X[i]) + c[i] * ((x - X[i])**2) + d[i] * ((x - X[i])**3)
            return y


def ex1():
    #A si b sunt matricile asociate functiei de la exercitiul 1
    A = np.array([[12.5 * 2, -20.0], [-20.0, 32.5 * 2]])
    #Matricea A este simetrica si pozitiv definita, deci admite un punct de minim unic
    b = np.array([4,8])

    #in continuare este calculat punctul de minim cu cele 2 metode
    x1 = pas_descendent(A, b)
    print(x1)

    x2 = gradienti_conjugati(A, b)
    print(x2)

def ex2():
    #capetele intervalului
    a = -np.pi
    b = np.pi
    n = 100

    X_ = np.linspace(a, b, n)
    Y_ = f2_(X_) #imaginea intervalului prin functia data
    Y_c = np.zeros(n) #imaginea intervalului prin aproximare

    N = 3 #pasul de start
    x_ = np.linspace(a, b, N) #diviziune a intervalului dat cu N elemente
    y_ = f2_(x_) #valoarea prin functie a acestei diviziuni
    #Cautam N astef incat polinomul c1 + c2(x-x2) + c3(x-x1)(x-x2) + .. + cn(x-x1)..(x-xn-1)
    #sa fie o aproximare cat mai buna a functiei date

    C = Lagrange_Newton(x_,y_,N) #calculeaza coeficientii c1..cn folosind metodata Newton
    for i in range(0,n):
        Y_c[i] = (interpolare_LN_elem(C, x_, N, X_[i])) #calculeaza valoarea unui element prin polinomul descris mai sus

    #se incrementeaza N si se repeta procesul cat timp diferenta maxima dintre valoarea unui element
    #din interval prin functia data si valoarea acestuia prin aproximare este mai mare decat epsilon
    while(np.max(abs(Y_c - Y_)) > eps):
        N += 1
        x_ = np.linspace(a, b, N)
        y_ = f2_(x_)
        C = Lagrange_Newton(x_, y_, N)
        for i in range(0, n):
            Y_c[i] = (interpolare_LN_elem(C, x_, N, X_[i]))

    print(N)
    print(np.max(abs(Y_c - Y_))) #eroarea de trunchiere

    #reprezentarea grafica a functiei si cea a aproximarii
    plt.figure(0)  # Initializare figura
    plt.plot(X_, Y_, linestyle='-', linewidth=3)
    plt.plot(X_, Y_c, linestyle='-', linewidth=3)

    #nourile de interpolare
    for i in range(0,N):
        plt.scatter(x_[i], y_[i], s=50, c='black', marker='o')

    plt.axvline(0, c='black')  # Adauga axa OY
    plt.axhline(0, c='black')  # Adauga axa OX
    plt.xlabel('x')  # Label pentru axa OX
    plt.ylabel('f(x)')  # Label pentru axa OY
    plt.title('Interpolare Lagrange cu metoda Newton')  # Titlul figurii
    plt.grid(b=True)
    plt.show()

def ex3():
    #capetele intervalului
    A = -np.pi
    B = np.pi
    n = 100

    X_ = np.linspace(A, B, n)
    Y_ = f3_(X_) #imaginea intervalului prin functia data
    Y_c = np.zeros(n) #imaginea intervalului prin aproximare

    N = 200 #nu am reusit sa gasesc un N pentru a satisface eroarea de trunchiere asa ca am ales unul pentru care aproximarea este cat mai buna
    x_ = np.linspace(A, B, N) #diviziune a intervalului dat cu N elemente
    y_ = f3_(x_) #valoarea prin functie a acestei diviziuni

    #metodata calculeaza listete de coeficienti a,b,c,d astefel incat sa formeze
    #polinoamele ai + bi(x-xi) + ci(x-xi)^2 + di(x-xi)^3 cu i de la 0 la N
    (a,b,c,d) = spline_cubice(x_,y_,N)
    for i in range(0,n):
        Y_c[i] = (interpolare_sc_elem(x_, a,b,c,d, N, X_[i])) #calculeaza valoarea unui element prin polinomul descris mai sus

    print(np.max(abs(Y_c - Y_))) #eroarea de trunchiere

    plt.figure(0)  # Initializare figura
    plt.plot(X_, Y_, linestyle='-', linewidth=3)
    plt.plot(X_, Y_c, linestyle='-', linewidth=3)

    # nourile de interpolare
    for i in range(0, N):
        plt.scatter(x_[i], y_[i], s=50, c='black', marker='o')

    plt.axvline(0, c='black')  # Adauga axa OY
    plt.axhline(0, c='black')  # Adauga axa OX
    plt.xlabel('x')  # Label pentru axa OX
    plt.ylabel('f(x)')  # Label pentru axa OY
    plt.title('Interpolare cu functii spline cubice')  # Titlul figurii
    plt.grid(b=True)
    plt.show()

#ex1()
#ex2()
#ex3()