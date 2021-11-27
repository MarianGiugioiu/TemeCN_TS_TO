import numpy as np
import matplotlib.pyplot as plt


def metoda_bisectiei(a, b, f, epsilon):
    """
    Implementare a metodei bisectiei.
    """
    assert a < b, 'a trebuie sa fie mai mic strict ca b!'
    x_num = (a + b) / 2  # Prima aproximare
    N = int(np.floor(np.log2((b - a) / epsilon)))  # Criteriul de oprire
        # int pentru ca `range` de mai jos nu accepta float

    # ==============================================================
    # Iteratiile algoritmului
    # ==============================================================
    for _ in range(1, N):
        if f(x_num) == 0:
            break
        elif np.sign(f(a)) * np.sign(f(x_num)) < 0:
            b = x_num
        else:
            a = x_num

        x_num = (a + b) / 2

    return x_num, N


# ==============================================================
# Aplicare a metodei bisectiei pentru un exemplu
# ==============================================================
f = lambda x: x**3 - 7*(x**2) + 14*x - 6  # Declararea functiei
A = 0  # Capatul din stanga al intervalului
B = 4  # Capatul din dreapta al intervalului
epsilon = 1e-5  # Eroarea maxima acceptata
x = np.linspace(A, B, 50)  # Discretizare a intervalului [A, B]
y = f(x)  # Valorile functiei pentru punctele din discretizare

a = [0., 1., 3.2]  # Capetele din stanga ale intervalelor
b = [1., 3.2, 4.]  # Capetele din dreapta ale intervalelor
assert len(a) == len(b), 'Capetele intervalelor trebuie sa aiba acelasi numar de elemente!'


# ==============================================================
# Afisare grafic
# ==============================================================
plt.figure(0)  # Initializare figura
plt.plot(x, y, linestyle='-', linewidth=3)  # Plotarea functiei

# Trece prin cele trei intervale
for i in range(len(a)):
    # Calculeaza solutia numerica si numarul de iteratii
    x_num, N = metoda_bisectiei(a=a[i], b=b[i], f=f, epsilon=epsilon)

    # Printeaza la consola rezultatele
    print('Metoda Bisectiei')
    print('Ecuatia: x^3 - 7*x^2 + 14*x - 6 = 0')
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
plt.title('Metoda Bisectiei')  # Titlul figurii
plt.axis('scaled')
plt.grid(b=True)
plt.show()  # Arata graficul
