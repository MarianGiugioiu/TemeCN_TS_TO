import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(precision=2)


def subs_desc_fast(a, b):

    """ Initializarea vectorului solutiei numerice. """
    n = b.shape[0] - 1
    x_num = np.zeros(shape=n + 1)

    """ Determinarea solutiei numerice. """
    x_num[n] = b[n] / a[n, n]  # Scrie ultima componenta a solutiei
    for k in range(n - 1, -1, -1):
        s = np.dot(a[k, k + 1:], x_num[k + 1:])
        x_num[k] = (b[k] - s) / a[k, k]

    return x_num


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
        print(a_ext);
    print(a_ext);
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


def meg_inversa(A):
    """Pentru a obtine inversa matricei A (A este inversabila daca determinantul ei este nenul) se va folosi matricea identitate I cu acelasi numar de linii 'n' ca matricea A"""
    """Se vor rezolva 'n' sisteme in care matricea extinsa este formata din matricea A la care se adauga cate o coloana din I"""
    """Matricea cu coloanele formate din sirul de solutii calculate mai sus va fi inversa matricii A"""
    n = A.shape[0]
    I = np.identity(n)
    return meg_pivotare_totala(A,I)

def factorizare_LU_GPP(a,b):
    """Algoritmul urmeaza sa calculeze U(initial matricea data) si L(matricea identitate) astefel incat a = LU"""
    n = a.shape[0]
    U = np.array(a)
    L = np.identity(n)
    """La fiecare pas al algoritmului de eliminare Gauss se vor aplica asupra matricii identitate aceleasi schimbari care s-au aplicat si matricii date"""
    """Matricea rezultata descrisa mai sus se va salva la fiecare pas in G"""
    """In final U va fi matricea initiala peste care s-a aplicat algoritmul de eliminare Gauss iar L va fi produsul inverelor matricilor calculate in G"""
    G = np.empty([n - 1, n, n])
    """'w' va retine interschimbarile dintre linii pe parcursul algoritmului"""
    w = np.arange(0, n, 1)
    for k in range(n-1):
        """M este matricea adaugata in G la pasul curent"""
        M = np.identity(n)
        """Aflam pozitia pivotului de pe colona k"""
        p = np.argmax(np.absolute(U[k:, k]))  # Pozitia celui mai mare element de pe coloana k de la linia k in jos
        p += k
        """ Schimba linia 'k' cu 'p' daca pivotul nu se afla pe diagonala principala. """
        if k != p:
            U[[p, k], :] = U[[k, p], :]
            w[[p, k]] = w[[k, p]]
            M[[p, k], :] = M[[k, p], :]
        """ Zero pe coloana sub pozitia pivotului. """
        for j in range(k + 1, n):
            m = U[j, k] / U[k, k]
            U[j, :] -= m * U[k, :]
            M[j, :] -= m * M[k, :]
        G[k] = M
        L = np.matmul(L,meg_inversa(M))

    """Trebuie sa aplic permutarea calculata in w asupra lui L si b"""
    L=L[w, :]
    B = b[w]
    """In acest punct L este o matrice inferior triunghiulara iar U superior triunghiulara"""
    """Vor fi rezolvate sistemele L*y=B si U*x=y, x fiind solutia cautata"""
    y = subs_asc_fast(L,B)
    x = subs_desc_fast(U,y)
    return x

def metoda_Cholesky(a):
    """Va fi calculata matricea L astefel incat produsul dintre ea si transpusa ei sa dea matricea initiala folosind algoritmul de la curs"""
    n = a.shape[0]
    L = np.zeros([n,n])
    L[0,0] = np.sqrt(a[0,0])
    for k in range(1,n):
        L[k,0] = a[k,0] / L[0,0]
    for k in range(1, n):
        m = a[k,k] - np.dot(L[k, 0:k],L[k, 0:k])
        L[k,k] = np.sqrt(m)
        for i in range(k+1,n):
            L[i,k] = (1/L[k,k])*(a[i,k] - np.dot(L[i, 0:k],L[k, 0:k]))

    return L


def ex1():
    A = np.array([
        [3., 8., 5.],
        [3., 28., 23.],
        [3., 3., 1.]
    ])
    b = np.array([18., 76., 4.])
    assert A.shape[0] == A.shape[1]
    assert A.shape[0] == b.shape[0]
    """Sistemul este compotabil daca determinantul este nenul"""
    if np.linalg.det(A) != 0:
        x = meg_pivotare_totala(A,b)

        print(x)
        #print(np.matmul(A, x))
    else:
        raise AssertionError('Sistem incompatibil sau sistem compatibil nedeterminat')

def ex2():
    A = np.array([
            [0., 9., -4., 0.],
            [0., 4., -4., -5.],
            [7., 1., -6., -4.],
            [5., 2., -1., 1.]
        ])
    assert A.shape[0] == A.shape[1]
    if np.linalg.det(A) != 0:
        B = meg_inversa(A)
        print(B)
        #print(np.matmul(A, B))
    else:
        raise AssertionError('Matricea nu este inversabila')

def ex3():
    A = np.array([
        [0., -2., 0., -1.],
        [4., -1., -8., -5.],
        [2., -9., 0., -9.],
        [5., 2., 8., -7.]
    ])
    b = np.array([-14., -62., -84., 21.])
    assert A.shape[0] == A.shape[1]
    assert A.shape[0] == b.shape[0]
    if np.linalg.det(A) != 0:
        x = factorizare_LU_GPP(A,b)
        print(x)
        #print(np.matmul(A,x))
    else:
        raise AssertionError('Sistem incompatibil sau sistem compatibil nedeterminat')

def ex4():
    A = np.array([
        [4., 2., -10., 18.],
        [2., 101., 55., -51.],
        [-10., 55., 86., -126.],
        [18., -51., -126., 298.]
    ])
    assert A.shape[0] == A.shape[1]
    """Determinantul matricei trebuie sa fie pozitiv"""
    if np.linalg.det(A) > 0:
        L = metoda_Cholesky(A)
        print(L)
        #Lt = L.transpose()
        #print(Lt)
        #print(np.matmul(L,Lt))
    else:
        raise AssertionError('Matricea nu admite factorizare Cholesky ')

ex1()
#ex2()
#ex3()
#ex4()
