import nonLinearClass as nl
import matplotlib.pyplot as plt
import numpy as np
import math as mt

def plot(a,b,x1,x2):

    N  = 100
    x  = np.linspace(x1,x2,N)
    gx = np.zeros(N,dtype=float)
    bx = np.zeros(N,dtype=float)
    fx = np.zeros(N,dtype=float)
    i = 0
    for xi in x:
      gx[i] = a.value(xi)*xi
      bx[i] = b.value(xi)
      fx[i] = bx[i] - gx[i]
      i+=1

    plt.xlabel('x')
    plt.xlim(0,x2)
    plt.ylim(0,max(max(gx),max(bx)))
    plt.plot(x ,gx, label='A(x)x')
    plt.plot(x ,bx, label='b(x)')
    plt.plot(x ,fx, label='b(x) - A(x)x')

    plt.legend()
    plt.show()


def main():

    defi1 = 'f(x) = a1/ (x ** a3 + a2)'
    def f1(a,x):

        a1 = a[0]
        a2 = a[1]
        a3 = a[2]

        return a1 / (x ** a3 + a2)

    def df1dx(a,x):

        a1 = a[0]
        a2 = a[1]
        a3 = a[2]

        return (-a1*a3*x**(a3-1.0))/((x**a3 + a2)**2)


    defi2 = 'f(x) = a1-a2*log(x)'
    def f2(a,x):

        a1 = a[0]
        a2 = a[1]
        a3 = a[2]

        return a1-a2*mt.log(a3*x)

    def df2dx(a,x):

        a1 = a[0]
        a2 = a[1]
        a3 = a[2]

        return -a2/x

    A=nl.Func((1.0,1.0,1.0), f1, df1dx, defi1)
    b=nl.Func((0.9,0.00,1.0), f2, df2dx, defi2)

    sist=nl.AxEqB(A,b)

    plot(A,b,0.1,20.0)

    print(sist)
    x0 = 0.01
    sist.picard(x0=x0,fPlot=False)
    sist.quasiNewton   (x0=x0,fPlot=False)
    sist.quasiNewtonMod(x0=x0,fPlot=False,maxIt=5000)
    sist.newton(x0=x0,fPlot=False,fExact=False)
    sist.newtonMod(x0=x0,fPlot=False,maxIt=5000)
    sist.newtonSec(x0=x0,fPlot=False)


if __name__ == '__main__':

    main()