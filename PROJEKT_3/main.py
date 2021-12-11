import pandas as pd
import time
from matplotlib import pyplot as plt
import numpy as np

#Łukasz Niedźwiadek 180102
def rysuj_wykres(tryb, lokacja, wezly_dystans, wezly_wysokosc, ilosc_wezlow, dystans, wysokosc, dane_dystans, dane_wysokosc):
    if tryb == 1:
        metoda = 'Lagrange'
    else:
        metoda = 'Splines'

    fig = plt.figure()
    plt.title(str(metoda) + ' ' + str(ilosc_wezlow) +
              ' węzłów ' + str(lokacja))
    plt.xlabel("dystans [m]")
    plt.ylabel("wysokość [m]")
    plt.plot(dane_dystans, dane_wysokosc,
             color='red', label='interpolated values')
    plt.plot(dystans, wysokosc, color='blue', label='actual values')
    plt.plot(wezly_dystans, wezly_wysokosc, 'o', color='blue')
    plt.savefig(str(metoda) + ' ' + str(ilosc_wezlow) +
                'wezlow ' + str(lokacja) + '.png')
    plt.close(fig)


def pivoting(U, L, P, i):

    #szukanie pivota
    pivot = abs(U[i][i])
    pivot_index = i

    for j in range(i+1, len(U)):
        if abs(U[j][i]) > pivot:
            pivot = abs(U[j][i])
            pivot_index = j

    if pivot_index != i:
        for j in range(0, len(U)):
            if j >= i:
                tmp = U[i][j]
                U[i][j] = U[pivot_index][j]
                U[pivot_index][j] = tmp
            else:
                tmp = L[i][j]
                L[i][j] = L[pivot_index][j]
                L[pivot_index][j] = tmp

            tmp = P[i][j]
            P[i][j] = P[pivot_index][j]
            P[pivot_index][j] = tmp

    return U, L, P


def faktoryzacja_LU(A, b, x, N):

    L = [[0 for x in range(N)] for Y in range(N)]
    U = [[0 for x in range(N)] for Y in range(N)]
    P = [[0 for x in range(N)] for Y in range(N)]
    Y = [1 for x in range(N)]

    for i in range(N):
        for j in range(N):
            if i == j:
                P[i][j] = 1
                L[i][j] = 1
            U[i][j] = A[i][j]

    for i in range(N - 1): #tworzę macierze L oraz U (A = LU)
        U, L, P = pivoting(U, L, P, i)
        for j in range(i+1, N):
            L[j][i] = U[j][i]/U[i][i]
            for k in range(i, N):
                U[j][k] = U[j][k] - (L[j][i]*U[i][k])

    b = np.matmul(P, b)

    for i in range(N): #rozwiązywanie Ly = b za pomocą podstawienia wprzód
        suma = 0
        for j in range(i):
            suma = suma + L[i][j] * Y[j]
        Y[i] = (b[i]-suma) / L[i][i]

    for i in range(N-1, -1, -1): #rozwiązywanie Ux = y za pomocą podstawienia wstecz
        suma = 0
        for j in range(i+1, N):
            suma = suma + (U[i][j] * x[j])
        x[i] = (Y[i]-suma)/U[i][i]


def splines(dane, dystans, ilosc_wezlow):
    N = 4*(ilosc_wezlow-1)
    M = [[0 for i in range(N)] for j in range(N)]
    b = [0 for i in range(N)]
    x = [1 for i in range(N)]

    b = [0]*N
    x = [0]*N

    for i in range(N):
        b[i] = 0
        x[i] = 1

    #tworzę układ równań
    M[0][0] = 1
    b[0] = dane[0][1]

    h = dane[1][0] - dane[0][0]
    M[1][0] = 1
    M[1][1] = h
    M[1][2] = h**2
    M[1][3] = h**3
    b[1] = dane[1][1]

    M[2][2] = 1
    b[2] = 0

    h = dane[ilosc_wezlow - 1][0] - dane[ilosc_wezlow - 2][0]
    M[3][4*(ilosc_wezlow - 2) + 2] = 2
    M[3][4*(ilosc_wezlow - 2) + 3] = 6 * h
    b[3] = 0

    for i in range(1, ilosc_wezlow-1):
        h = dane[i][0] - dane[i-1][0]

        M[4*i][4*i] = 1
        b[4*i] = dane[i][1]

        M[4*i + 1][4*i] = 1
        M[4*i + 1][4*i + 1] = h
        M[4*i + 1][4*i + 2] = h ** 2
        M[4*i + 1][4*i + 3] = h ** 3
        b[4*i + 1] = dane[i + 1][1]

        M[4*i + 2][4*(i - 1) + 1] = 1
        M[4*i + 2][4*(i - 1) + 2] = 2 * h
        M[4*i + 2][4*(i - 1) + 3] = 3 * h**2
        M[4*i + 2][4*i + 1] = -1
        b[4*i + 2] = 0

        M[4*i + 3][4*(i - 1) + 2] = 2
        M[4*i + 3][4*(i - 1) + 3] = 6*h
        M[4*i + 3][4*i + 2] = -2
        b[4*i + 3] = 0

    faktoryzacja_LU(M, b, x, N)

    elevation = 0
    for i in range(ilosc_wezlow-1):
        elevation = 0
        if dystans >= dane[i][0] and dystans <= dane[i+1][0]:
            for j in range(4):
                h = dystans - dane[i][0]
                elevation += x[4*i + j] * h**j
            break
    return elevation


def lagrange(tuples_wezly, dystans, ilosc_wezlow):
    wynik = 0
    for i in range(ilosc_wezlow):
        a = 1
        for j in range(ilosc_wezlow):
            if i != j:
                a = a * (dystans - tuples_wezly[j][0]) / \
                    (tuples_wezly[i][0] - tuples_wezly[j][0])

        wynik = wynik + a * tuples_wezly[i][1]
    return wynik


def interpolacja(dane, ilosc_wezlow, tuples_wezly, nazwa, metoda):

    wyniki = []
    tic = time.time()
    start = tuples_wezly[0]
    stop = tuples_wezly[ilosc_wezlow-1]
    for i in range(int(start[0]), int(stop[0]), 8):
        flaga = True
        for j in range(ilosc_wezlow):
            if int(tuples_wezly[j][0]) == i: #przypadek w którym wartości nie muszą być interpolowane
                flaga = False
                wynik = tuples_wezly[j][1]
                break
        if flaga == True:
            if metoda == 1:
                wynik = lagrange(tuples_wezly, i, ilosc_wezlow)
            else:
                wynik = splines(tuples_wezly, i, ilosc_wezlow)

        wyniki.append(wynik)

    toc = time.time()
    czas = toc - tic

    dystans = []
    for i in range(len(wyniki)):
        dystans.append(i*8)
    wezly_dystans, wezly_wysokosc = zip(*tuples_wezly)
    dane_dystans, dane_wysokosc = zip(*dane)
    rysuj_wykres(metoda, nazwa, wezly_dystans, wezly_wysokosc,
                 ilosc_wezlow, dystans, wyniki, dane_dystans, dane_wysokosc)
    return czas


def wczytaj_dane(nazwa):
    probki = pd.read_csv('2018_paths\\' + nazwa + '.csv',
                         sep=',', usecols=['Dystans (m)', 'Wysokość (m)'])
    dane = []
    for i in range(len(probki)):
        tmp_x = float(probki['Dystans (m)'][i])
        tmp_y = float(probki['Wysokość (m)'][i])
        wspl_tuple = (tmp_x, tmp_y)
        dane.append(wspl_tuple)
    return dane


odstepy = [100, 50, 33, 17]
# wezly = [5,10,15,29] dla takich wartości tablicy odstępy uzyskam takie węzły
ilosc_zestawow_wezlow = 4

tuples_dane = []

wspolrzedne = (0, 0)
nazwy = ['SpacerniakGdansk', 'MountEverest', 'WielkiKanionKolorado']
ilosc_danych = 3
l_danych = 500

for j in range(0, ilosc_danych):  # iteracja po lokacjach
    tuples_dane = wczytaj_dane(nazwy[j])
    print(nazwy[j])
    for i in range(0, ilosc_zestawow_wezlow):  # iteracja po ilości węzłów
        krok = odstepy[i]
        ilosc_wezlow = int(l_danych/krok)
        tuples_wezly = []
        k = 0
        for p in range(ilosc_wezlow):  # iteracja dla węzła
            wspolrzedne = tuples_dane[k]
            tuples_wezly.append(wspolrzedne)
            k = k + krok
        czas = interpolacja(tuples_dane, ilosc_wezlow,
                            tuples_wezly, nazwy[j], 1)
        print('Ilość węzłów: ' + str(ilosc_wezlow) +
              ' Czas działania [s] ' + str(czas))
