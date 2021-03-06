import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-5, 5, 100)

def ex1():
    # max x - 2y
    plt.figure(dpi=100)

    plt.xlim(-5, 5)
    plt.ylim(-5, 5)

    # x + 10y >= 25
    y = -1/10 * x + 5/2
    plt.plot(x, y)
    plt.fill_between(x, y, y + 10, alpha=0.6)

    # 5x + y <= -10
    y = -5*x - 10
    plt.plot(x, y)
    plt.fill_between(x, y, y - 20, alpha=0.6)

    # -9x + 7y <= -17
    y = 9/7 * x - 17/7
    plt.plot(x, y)
    plt.fill_between(x, y, y - 10, alpha=0.6)

    #plt.savefig('plot1.png')
    plt.show()

def ex2():
    # max -2x + 3y
    plt.figure(dpi=100)

    plt.xlim(-5, 5)
    plt.ylim(-5, 5)

    # y >= -1
    y = 0*x - 1
    plt.plot(x, y)
    plt.fill_between(x, y, y + 10, alpha=0.6)

    # x + y <= 1
    y = (-x + 1)
    plt.plot(x, y)
    plt.fill_between(x, y, y - 20, alpha=0.6)

    # -x + y <= -3
    y =  x - 3
    plt.plot(x, y)
    plt.fill_between(x, y, y - 10, alpha=0.6)

    #plt.savefig('plot2.png')
    plt.show()

def ex3():
    # max 7x + 8y
    plt.figure(dpi=100)

    plt.xlim(-5, 5)
    plt.ylim(-5, 5)

    # 7x + 8Y <= -19
    y = (-7*x - 19)/8
    plt.plot(x, y)
    plt.fill_between(x, y, y - 20, alpha=0.6)

    # 7x + 8Y >= -19
    y = (-7*x - 19)/8
    plt.plot(x, y)
    plt.fill_between(x, y, y + 20, alpha=0.6)

    #plt.savefig('plot3.png')
    plt.show()

def ex4():
    # max x + y
    plt.figure(dpi=100)

    plt.xlim(-5, 5)
    plt.ylim(-5, 5)

    # 3x + y >= -2
    y = -3 * x - 2
    plt.plot(x, y)
    plt.fill_between(x, y, y + 30, alpha=0.6)

    #plt.savefig('plot4.png')
    plt.show()


#ex1()
#ex2()
#ex3()
#ex4()

def ex5() :
    x = np.linspace(0, 25, 100)
    plt.figure(dpi=100)

    plt.xlim(0, 25)
    #plt.ylim(-5, 5)

    # 3x + y >= -2
    y = np.pi / (1000 * 1.2 ** (-x))
    plt.plot(x, y)

    # plt.savefig('plot4.png')
    plt.show()

ex2()