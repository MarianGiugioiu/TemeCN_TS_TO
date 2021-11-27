import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(precision=2)

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



def meg_inversa(A):
    """Pentru a obtine inversa matricei A (A este inversabila daca determinantul ei este nenul) se va folosi matricea identitate I cu acelasi numar de linii 'n' ca matricea A"""
    """Se vor rezolva 'n' sisteme in care matricea extinsa este formata din matricea A la care se adauga cate o coloana din I"""
    """Matricea cu coloanele formate din sirul de solutii calculate mai sus va fi inversa matricii A"""
    n = A.shape[0]
    I = np.identity(n)
    return meg_pivotare_totala(A,I)

def tema2():
    B = np.array([
        [1, 3, 2],
        [0, -2, 1],
        [2, 0, 1]
    ], dtype=float)

    k = 2
    C = np.array([-1, 2, 3], dtype=float)

    print("Matricea aleasa B este:")
    print(B)

    if np.linalg.det(B) == 0:
        raise AssertionError('Matricea nu este inversabila')

    B_inv = meg_inversa(B)
    print("Inversa matricii B este:")
    print(B_inv)

    print("Vectorul ales C este:")
    print(C)
    print("K =", k)

    B_tilde = B.copy()
    B_tilde[:, k] = C

    print("Matricea B_tilde este:")
    print(B_tilde)

    y = np.matmul(B_inv,C)
    print(y)
    if y[k] == 0:
        raise AssertionError('Matricea B_tilde nu este inversabila')

    eta = -y / y[k]
    eta[k] = 1 / y[k]

    print(eta)

    N = B.shape[0]
    E = np.eye(N)
    E[:, k] = eta

    B_tilde_inv = np.matmul(E,B_inv)

    print("Inversa lui B_tilde:")
    print(B_tilde_inv)
    print(np.linalg.inv(B_tilde))

tema2()

