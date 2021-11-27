import numpy as np
import matplotlib.pyplot as plt


def metoda_newton_raphson(f, df, x_old, epsilon):
    """ Varianta a implementarii metodei N-R. """

    x_new = x_old - f(x_old)/df(x_old)  # Prima aproximare
    counter = 1  # Numar de pasi
    while np.abs(f(x_new)) > epsilon:  # Primul criteriu de oprire
        x_old = x_new
        x_new = x_old - f(x_old) / df(x_old)
        counter += 1

    return x_new, counter


# ==========================================================================
# Aplicare pe exercitiul din laborator pentru eps=10^{-5}
# ==========================================================================
def f_(x):
    """ Functia data. """
    return np.cos(x) - x


def df_(x):
    """ Derivata functiei date. """
    return -np.sin(x) - 1


a = 0  # Capat stanga interval
b = np.pi/2  # Capat dreapta interval
x_ = np.linspace(a, b, 100)  # Discretizare interval [a, b]
eps = 1e-5  # Eroarea acceptata
y_ = f_(x_)  # Valorile functiei pe discretizarea intervalului

x0 = 1.5  # Punct de start

x_num, N = metoda_newton_raphson(f=f_, df=df_, x_old=x0, epsilon=eps)
print(x_num)
plt.figure(0)
plt.plot(x_, y_, linestyle='-', linewidth=3, c='orange')  # Plotarea functiei
plt.scatter(x_num, 0, s=50, c='black', marker='o')
plt.legend(['f(x)', 'x_num'])
plt.axvline(0, c='black')  # Adauga axa OY
plt.axhline(0, c='black')  # Adauga axa OX
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Metoda Newton-Raphson, N = {}'.format(N))
plt.grid(b=True)
plt.axis('scaled')
plt.show()

