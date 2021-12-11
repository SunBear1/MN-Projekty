import math
import time
import matplotlib.pyplot as plt

#Łukasz Niedźwiadek 180102

def norm(vector):
    x = 0
    for i in range(N):
            x = x + math.pow(vector[i], 2)
    return math.sqrt(x)

def residuum(A,v,b,N):
    res = [0]*N
    for i in range(N):
        for j in range(N):
            res[i] += A[i][j] * v[j] #wykonanie mnożenia res = M*r
    for i in range(N):
        res[i] = res[i] - b[i] #wykonanie odejmowania od res = res-b
    return res

def Jacobi(A,b,N):
    liczba_iteracji = 0
    r = [0]*N
    r_stare = [1]*N
    tic = time.time()
    norm_res = 1
    prog = 10**-9
    while norm_res >= prog:
        for i in range(N):
            LUR = 0 #(L + U)*r(k)
            for j in range(N):
                if(j != i): #nie diagonal
                    LUR = LUR + (A[i][j] * r_stare[j]) #(L + U)*r(k)
            r[i] = (b[i] - LUR) / A[i][i] #r = (LUR - b) / D
        liczba_iteracji = liczba_iteracji + 1
        for k in range(N):
            r_stare[k] = r[k]
        res = residuum(A,r,b,N)
        norm_res = norm(res)
        if math.isinf(norm_res) == 1: #warunek zbieżności
            print("nie zbiega się")
            break
        elif math.isnan(norm_res) == 1: #zabezpieczający warunek zbieżności
            print("nie zbiega się")
            break
    toc = time.time()
    czas = toc - tic
    return czas,liczba_iteracji,norm_res

def GaussSeidl(A,b,N):
    liczba_iteracji = 0
    r = [0]*N
    r_stare = [1]*N
    res = [1]*N
    tic = time.time()
    norm_res = 1
    prog = 10**-9
    while norm_res >= prog:
        for i in range(N):
            UR = 0 #U*r(k)
            L = 0
            LUR = 0
            for j in range(i+1,N):
                UR = UR + A[i][j] * r_stare[j] #U*r(k)
            for j in range(i):
                L = L + A[i][j] * r[j] #L
            LUR = LUR + L + UR #L + Ur
            r[i] = (b[i] - LUR) / A[i][i] #r = (b - Ur) / DL
        liczba_iteracji = liczba_iteracji + 1
        for k in range(N):
            r_stare[k] = r[k]
        res = residuum(A,r,b,N)
        norm_res = norm(res)
        if math.isinf(norm_res) == 1: #warunek zbieżności
            print("nie zbiega się")
            break
        elif math.isnan(norm_res) == 1: #zabezpieczający warunek zbieżności
            print("nie zbiega się")
            break
    toc = time.time()
    czas = toc - tic
    return czas,liczba_iteracji,norm_res

def FaktoryzacjaLU(A,b,N):
    tic = time.time()
    L = [[0 for x in range(N)] for y in range(N)] 
    U = [[0 for x in range(N)] for y in range(N)] 
    for i in range(N):
        L[i][i] = 1 #początkowe stworzenie macierzy jednostkowej L
        for j in range(N):
            U[i][j] = A[i][j] #początkowe stworzenie macierzy U=A

    for i in range(N):
        for j in range(i+1,N):
            L[j][i] = U[j][i] / U[i][i] #wypełnienie macierzy L
            for k in range(i,N):
                U[j][k] = U[j][k] - L[j][i] * U[i][k] # wypełnienie macierzy U

    Y = [1]*N
    Y[0] = b[0]/L[0][0]
    for i in range(1,N): #rozwiązywanie Ly = b za pomocą podstawienia wprzód
        suma = 0
        for j in range(i):
            suma = suma + L[i][j]*Y[j]
        Y[i] = (b[i]-suma) / L[i][i]

    X = [1]*N
    X[N-1] = b[N-1]/U[N-1][N-1]
    for i in range(N-1,-1,-1): #rozwiązywanie Ux = y za pomocą podstawienia wstecz
        suma = 0
        for j in range(N-1,i,-1):
            suma = suma + (U[i][j]*X[j])
        X[i] = (Y[i]-suma) / U[i][i]

    toc = time.time()
    czas = toc - tic
    res = residuum(A,X,b,N)
    norm_res = norm(res)
    return czas,norm_res
        
idx = 180102
c = 0
d = 2
e = 1
f = 0
a1 = 5 + e
a2 = a3 = -1
N = 900 + c*10 + d
N_tab = [100,500,1000,2000,3000]
Osx = [1,2,3,4,5]
czasy_jacobi = []
czasy_gs = []
czasy_LU = []
b = []

""" ZADANIE C
a1 = 3
for i in range(N):
    b.append(math.sin(i*(f+1)))
A = [[0 for x in range(N)] for y in range(N)] 
for i in range(N):
    for j in range(N):
        if j == i:
            A[i][i] = a1
            if j > 0:
                A[i][j-1]=a2
            if j > 1:
                A[i][j-2]=a3
            if j < N-1:
                A[i][j+1]=a2
            if j < N-2:
                A[i][j+2]=a3
print("ZADANIE C(N = 902, a1 = 3, a2 = a3 = -1")
czas,iterations,norm_res = Jacobi(A,b,N)
print("Metoda Jacobiego")
print("Ilość iteracji: ", iterations)
print("Czas działania [s]: ", czas)
print("Norma z residuum: ", norm_res)
czas,iterations,norm_res = GaussSeidl(A,b,N)
print("Metoda GaussaSeidla")
print("Ilość iteracji: ", iterations)
print("Czas działania [s]: ", czas)
print("Norma z residuum: ", norm_res)
"""
#ZADANIE E
for k in range(len(N_tab)):
    N = N_tab[k]
    for i in range(N):
        b.append(math.sin(i*(f+1)))
    A = [[0 for x in range(N)] for y in range(N)] 
    for i in range(N):
        for j in range(N):
            if j == i:
                A[i][i] = a1
                if j > 0:
                    A[i][j-1]=a2
                if j > 1:
                    A[i][j-2]=a3
                if j < N-1:
                    A[i][j+1]=a2
                if j < N-2:
                    A[i][j+2]=a3
    czas,iterations,norm_res = Jacobi(A,b,N)
    czasy_jacobi.append(czas)
    czas,iterations,norm_res = GaussSeidl(A,b,N)
    czasy_gs.append(czas)
    czas,norm_res = FaktoryzacjaLU(A,b,N)
    czasy_LU.append(czas)

fig = plt.figure(num=None, figsize=(23, 10), dpi=80, facecolor='w', edgecolor='k')
plt.plot(N_tab, czasy_jacobi, '-', color='green')
plt.plot(N_tab, czasy_gs, '-', color='red')
plt.plot(N_tab, czasy_LU, '-', color='blue')
plt.legend(['Jacobi', 'Gauss-Seidl','FaktoryzacjaLU'])
plt.ylabel("Czas [s]")
plt.title("Czas działania poszczególnych algorytmów w zależności od N")
fig.savefig('ZADANIE_E.png', dpi=fig.dpi)
plt.show()
    
"""
print("ZADANIE D(N = 902, a1 = 3, a2 = a3 = -1")
czas,norm_res = FaktoryzacjaLU(A,b,N)
print("Metoda FaktoryzacjiLU")
print("Czas działania [s]: ", czas)
print("Norma z residuum: ", norm_res)


print("ZADANIE B(N = 902, a1 = 6, a2 = a3 = -1")
czas,iterations,norm_res = Jacobi(A,b,N)
print("Metoda Jacobiego")
print("Ilość iteracji: ", iterations)
print("Czas działania [s]: ", czas)
print("Norma z residuum: ", norm_res)

czas,iterations,norm_res = GaussSeidl(A,b,N)
print("Metoda GaussaSeidla")
print("Ilość iteracji: ", iterations)
print("Czas działania [s]: ", czas)
print("Norma z residuum: ", norm_res)



czas,iterations,norm_res = GaussSeidl(A,b,N)
print("Metoda GaussaSeidla")
print("Ilość iteracji: ", iterations)
print("Czas działania [s]: ", czas)
print("Norma z residuum: ", norm_res)
"""
