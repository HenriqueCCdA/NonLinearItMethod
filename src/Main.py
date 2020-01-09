import nonLinearClass as nl


def main():

    defi1 = 'f(x) = a1/ (x ** a3 + a2)'
    def f1(a,x):

        a1 = a[0]
        a2 = a[1]
        a3 = a[2]

        return a1 / (x ** a3 + a2)


    defi2 = 'f(x) = a1'
    def f2(a,x):

        a1 = a[0]

        return a1

    A=nl.Func((1.0,1.0,2.0), f1, None, defi1)
    b=nl.Func((0.2,0.0,0.0), f2, None, defi2)

    sist=nl.AxEqB(A,b)

    print(sist)
    x0 = 1.0
    sist.picard        (x0=x0,fPlot=False)
    sist.quasiNewton   (x0=x0,fPlot=False)
    sist.quasiNewtonMod(x0=x0,fPlot=False)

if __name__ == '__main__':

    main()