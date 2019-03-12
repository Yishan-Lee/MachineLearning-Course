import random
import matplotlib.pyplot as plt
import numpy as np
from itertools import chain 

def LSE(A, XVec, data, target):
    err = 0
    XVec = [[XVec[i]] for i in range(len(XVec))]
    y = mult_matrix(A, XVec)
    for i in range(len(data)):
        err += (y[i][0] - target[i])**2
    print("Total error : ",err)
    return err

def Equation(w, N): 


    x=list(chain.from_iterable(X))
#    print(x)
#    x=[-5.0,-4.795918,-4.5918,-3.979594,-3.571428,-2.959183,-2.75510,-1.734693] 
    if N == 1:
        formula = "y = %d" % w[0]
        yLine = w[0]
        plt.axhline(y=w[0], color='r', linestyle='-')
        plt.show()     
    else:
        formula = "y = %.8fx^%d" % (w[0], N-1)
        yLine =  w[0]*(np.power(x,N-1))
        for i in range(1, N):
            if w[i] < 0:
                formula += " - %.8fx^%d" % (-w[i], N-1-i)
                #yLine += yLine + w[i]*(np.power(x,N-1-i))
            else:
                formula += " + %.8fx^%d" % (w[i], N-1-i)
            yLine +=  w[i]*(np.power(x,N-1-i))
        plt.ylim(-10,100)
        plt.plot(x, yLine, 'g-')
        plt.show() 
 
    print(formula)
    return yLine

     
def Invert(a, L, U):
    #LY=B
    Y = [[0.0]*len(a) for i in range(len(a))]
    B = [[0.0]*len(a) for i in range(len(a))]

    for i in range(len(a)):
        B[i][i] = 1

    for j in range(len(a)):
        for i in range(len(a)):
            yy=B[i][j]
            for k in range(i):
                yy=yy-L[i][k]*Y[k][j]
            Y[i][j]=yy

    #UX=Y
    X = [[0.0]*len(a) for i in range(len(a))]
    for j in range(len(a)):
        for i in range(len(a)-1,-1,-1):
            x=Y[i][j]
            for k in range(len(a)-i-1):
                x=x-U[i][len(a)-1-k]*X[len(a)-1-k][j]
            X[i][j]=x/U[i][i]

    return X


def mult_matrix(N, M):

    row_N, col_N = len(N), len(N[0])
    row_M, col_M = len(M), len(M[0])

    A = [[0 for i in range(col_M)] for j in range(row_N)]
    for i in range(row_N):
        for j in range(col_M):
            for k in range(col_N):
                A[i][j] = A[i][j] + N[i][k] * M[k][j]

    return A


def lu_decomposition(A): 

    n = len(A)
    # Create zero matrices for L and U
    L = [[0.0] * n for i in range(n)] #2*2
    U = [[0.0] * n for i in range(n)] #2*2

    # Perform the LU Decomposition
    for j in range(n):
#        print("j: ",j)
        L[j][j] = 1.0
#        print("L[j][j]: ",L)
        for i in range(j+1):
#            print("i:",i)
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
#            print("k: ",k)
#            print("U[k][j] ",U)
#            print("L[i][k] ",L)
#            print("s1",s1)
            U[i][j] = A[i][j] - s1
#            print("U[i][j]",U)
        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (A[i][j] - s2) / U[j][j]

    return (L, U)



def NewtonMethod(ATA, AT, b, N,epsilon=1e-10):

    x = [random.randint(1,100)*0.1 for i in range(N)]    
    
    while(True):

        x = [[x[i]] for i in range(len(x))]
#        print(x)
        if(type(b[0])!= list):
            b = [[b[i]] for i in range(len(b))]
        ATAX = mult_matrix(ATA, x)
        Atb = mult_matrix(AT, b)
        ATAX_2 = [[0]*len(ATAX[0]) for i in range(len(ATAX))]
        Atb_2 = [[0]*len(Atb[0]) for i in range(len(Atb))]
        for i in range(len(ATAX)):
            for j in range(len(mult_matrix(ATA, x)[0])):
                ATAX_2[i][j] = mult_matrix(ATA, x)[i][j] * 2
        for i in range(len(Atb)):
           for j in range(len(Atb[0])):
               Atb_2[i][j] = Atb[i][j] * 2       
        
        first_der = [[0]*len(ATAX_2[0]) for i in range(len(Atb_2))]
        for i in range(len(ATAX_2)):
            for j in range(len(ATAX_2[0])):
                first_der[i][j] = ATAX_2[i][j] - Atb_2[i][j]

        second_der = [[0]*len(ATA[0]) for i in range(len(ATA))]
        for i in range(len(ATA)):
            for j in range(len(ATA[0])):
                second_der[i][j] = ATA[i][j] * 2

        L, U = lu_decomposition(second_der)
        xx = mult_matrix(Invert(second_der, L, U), first_der)
        
        x_new = [[0]*len(x[0]) for i in range(len(x))]
        for i in range(len(x)):
            for j in range(len(x[0])):
                x_new[i][j] = x[i][j] - xx[i][j]

        x_new = [x_new[i][0] for i in range(len(x_new))]

        flag = 1
        for i in range(N):
            if abs(x_new[i] - x[i][0]) > epsilon:
                flag = 0
        if flag:
            return x_new
        else:
            x = x_new
  
    
    
if __name__ == '__main__':

    fo = open("testfile.txt", "r")  
    x,y=[],[]
    for line in fo.readlines(): 
        line = line.strip() 
        ele = line.split(",")
        x.append(float(ele[0]))
        y.append(float(ele[1])) 
    fo.close()

    N = int(input('輸入Polynomial basis：'))
    lambda1 = int(input('輸入lambda： ')) #newton-lambda=0
    
    #create A, At, AtA
    MatrixA = [[0 for i in range(N)] for j in range(len(x))]  #23*2
    MatrixAT = [[0 for i in range(len(x))] for j in range(N)] #2*23
    MatrixATA_LSE = [[0 for i in range(N)] for j in range(N)] #2*2
    MatrixATA_NEW = [[0 for i in range(N)] for j in range(N)] #2*2

    # A & AtA
    for i in range(len(x)): #列 23
        for j in range(N): #次 2
            MatrixA[i][j] = x[i] ** (N-j-1) #23*2
            MatrixAT[j][i] = MatrixA[i][j] #2*23
    #ATA + λI
    for i in range(N):
        for j in range(N):
            for k in range(len(x)):
                MatrixATA_LSE[i][j] = MatrixATA_LSE[i][j] + MatrixAT[i][k] * MatrixA[k][j]
                MatrixATA_NEW[i][j] = MatrixATA_NEW[i][j] + MatrixAT[i][k] * MatrixA[k][j]
            if i == j: #對角線
                MatrixATA_LSE[i][j] = MatrixATA_LSE[i][j] + lambda1 

    #OPERATION
    L,U = lu_decomposition(MatrixATA_LSE)
    MatrixIvert = Invert(MatrixATA_LSE, L, U)
    MatrixX1 = mult_matrix(MatrixIvert, MatrixAT)
    Y = [[y[i]] for i in range(len(y))]
    X = [[x[i]] for i in range(len(x))]
    W = mult_matrix(MatrixX1, Y)

    #A-1*AT*y
    plt.plot(x, y, 'ro', label='data point')
    print("\n"+'--LSE--')
    weights = [W[i][0] for i in range(len(W))]
    yLine = Equation(weights,int(N))
    error = LSE(MatrixA, weights, x, y)
    
        
    print('--Newton--')
    plt.plot(x, y, 'ro', label='data point')
    Newtonweight = NewtonMethod(MatrixATA_NEW, MatrixAT, y, int(N), epsilon=1e-10)
    formula_New = Equation(Newtonweight,int(N))
    Newton_error = LSE(MatrixA, Newtonweight, x, y)

