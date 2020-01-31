class Func(object):
    '''
    Class que dedine uma funcao do tipo f(x)
    '''

# ---
    def __init__(self, cas,f,dfdx,defi):
        '''
        construtor da Class A
        '''
        self.cas  = cas  # coeficientes
        self.f    = f    # definicao da funcao
        self.dfdx = dfdx # definicao da derivada da funcao
        self.defi = defi # string que contem a forma matematica da funcao
# ----------------------------------------------------------------------------------------------------------------------

# ---
    def __str__(self):
        ss = ''
        for a in self.cas:
            ss = ss +' '+ str(a)

        s = '{0}\n a ={1}'.format(self.defi,ss)
        return s
# ----------------------------------------------------------------------------------------------------------------------

    def value(self, x):
        '''
        valor da funcao calculado em x
        '''
        return self.f(self.cas, x)
# ----------------------------------------------------------------------------------------------------------------------

# ---
    def deriv(self, x):
        '''
        valor da derivada da funcao calculado em x
        '''
        return self.dfdx(self.cas, x)
# ----------------------------------------------------------------------------------------------------------------------



class AxEqB(object):
    '''
    Class que dedine uma um sistema nao-linaer da forma A(x)x = b(x)
    '''

# ---
    def __init__(self, A, b):
        self.A = A
        self.b = b
# ----------------------------------------------------------------------------------------------------------------------

# ---
    def __str__(self):
        s1 = 'A(x) x = b(x):\n'
        s2 = 'A(x):\n' + self.A.__str__() + '\n'
        s3 = 'b(x):\n' + self.b.__str__()

        return s1 + s2 + s3
# ----------------------------------------------------------------------------------------------------------------------

    #   def jac(self, x):
    #       A = self.A.A(x)
    #       dadx = self.A.dAdx(self.A, x)
    #       dbdx = self.b.dBdx(x)
    #
    #       return dbdx - dadx * x - matrixA(x)

# ---
    def vetorF(self,x):
        return self.b.value(x) - self.A.value(x)*x
# ----------------------------------------------------------------------------------------------------------------------

# ---
    def tol(self,type=float):
        '''
        checa a precisao da maquina
        '''
        mq = type(1)
        while mq + type(1) != type(1):
            mq = mq/type(2)
        return mq*type(2)

# ----------------------------------------------------------------------------------------------------------------------

# ---
    def picard(self, x0=0.0, maxIt=1000, tol=1.e-14, fPlot = True):
        '''
        ***********************************************************
        data criacao:     09/01/2020
        data modificacao: 00/00/0000
        -----------------------------------------------------------
        metodo de picard: solucao do sistema nao-linear A(x)x=b(x)
        -----------------------------------------------------------
        Entrada:
        -----------------------------------------------------------
        x0 - chute inicial
        maxIt - numero maximo de iteracoes
        tol   - tolerancia desejada
        fPlot - plotagem do historico de convergencia
        -----------------------------------------------------------
        Saida:
        -----------------------------------------------------------
        x  - retorno o valor a raiz
        -----------------------------------------------------------
        OBS:
        -----------------------------------------------------------
        ***********************************************************
        '''

        if tol == 0.0:
            self.tol(type(x0))
        x = x0
        res0  = abs(self.vetorF(x))
        for i in range(0,maxIt+1):
            x   = self.b.value(x)/self.A.value(x)
            res = abs(self.vetorF(x))
            # ... polta o historico da convergencia do processo
            if fPlot :
                print('{0:4} {1:e} {2:e} {3:e}'.format(i+1,x,res,res/res0))
            # ... checa a convergencia
            if res/res0 < tol:
                break;


        res = abs(self.vetorF(x))

        print('{0:15s}: It = {1:3d}, x = {2:e}, A(x)x = {3:f}, b(x) = {4:f}, Res = {5:e}'.format('picard', i, x,
                                                                                                 self.A.value(x) * x,
                                                                                                 self.b.value(x), res))
# ----------------------------------------------------------------------------------------------------------------------

# ---
    def quasiNewton(self, x0=0.0, maxIt=1000, tol=1.e-14, fPlot = True):
        '''
        ***********************************************************
        data criacao:     09/01/2020
        data modificacao: 00/00/0000
        -----------------------------------------------------------
        metodo de quaisiNetwon: solucao do sistema nao-linear
        A(x)x=b(x)
        -----------------------------------------------------------
        Entrada:
        -----------------------------------------------------------
        x0 - chute inicial
        maxIt - numero maximo de iteracoes
        tol   - tolerancia desejada
        fPlot - plotagem do historico de convergencia
        -----------------------------------------------------------
        Saida:
        -----------------------------------------------------------
        x  - retorno o valor a raiz
        -----------------------------------------------------------
        OBS:
        -----------------------------------------------------------
        Matrix jacobina aproximada por -A(x)
        ***********************************************************
        '''

        if tol == 0.0:
            self.tol(type(x0))
        x = x0
        res  = self.vetorF(x)
        res0 = abs(res)
        for i in range(0,maxIt+1):
            # mariz jacobiana aproximada
            Jac = -self.A.value(x)
            dx  = -res/Jac
            x  += dx
            res = self.vetorF(x)
            # ... polta o historico da convergencia do processo
            if fPlot :
                print('{0:4} {1:e} {2:e} {3:e}'.format(i+1,x,res,res/res0))
            # ... checa a convergencia
            if abs(res)/res0 < tol:
                break;


        res = abs(self.vetorF(x))

        print('{0:15s}: It = {1:3d}, x = {2:e}, A(x)x = {3:f}, b(x) = {4:f}, Res = {5:e}'.format('quasiNewton', i, x,
                                                                                                 self.A.value(x) * x,
                                                                                                 self.b.value(x), res))
# ----------------------------------------------------------------------------------------------------------------------

# ---
    def quasiNewtonMod(self, x0=0.0, maxIt=1000, tol=1.e-14, fPlot = True):
        '''
        ***********************************************************
        data criacao:     09/01/2020
        data modificacao: 00/00/0000
        -----------------------------------------------------------
        metodo de quaisiNetwonMod: solucao do sistema nao-linear
        A(x)x=b(x)
        -----------------------------------------------------------
        Entrada:
        -----------------------------------------------------------
        x0 - chute inicial
        maxIt - numero maximo de iteracoes
        tol   - tolerancia desejada
        fPlot - plotagem do historico de convergencia
        -----------------------------------------------------------
        Saida:
        -----------------------------------------------------------
        x  - retorno o valor a raiz
        -----------------------------------------------------------
        OBS:
        -----------------------------------------------------------
        Matrix jacobina aproximada por -A(x0)
        ***********************************************************
        '''

        if tol == 0.0:
            self.tol(type(x0))
        x = x0
        res  = self.vetorF(x)
        res0 = abs(res)
        Jac = -self.A.value(x)
        for i in range(0,maxIt+1):
            # mariz jacobiana aproximada
            dx  = -res/Jac
            x  += dx
            res = self.vetorF(x)
            # ... polta o historico da convergencia do processo
            if fPlot :
                print('{0:4} {1:e} {2:e} {3:e}'.format(i+1,x,res,res/res0))
            # ... checa a convergencia
            if abs(res)/res0 < tol:
                break;

        res = abs(self.vetorF(x))

        print('{0:15s}: It = {1:3d}, x = {2:e}, A(x)x = {3:f}, b(x) = {4:f}, Res = {5:e}'.format('quasiNewton', i, x,
                                                                                                 self.A.value(x) * x,
                                                                                                 self.b.value(x), res))
# ----------------------------------------------------------------------------------------------------------------------
