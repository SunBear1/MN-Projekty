import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np


# Łukasz Niedźwiadek 180102

def EMA(N, start, tab):
    licznik = 0.0
    mianownik = 0.0
    for i in range(start - N, start):
        alpha = math.pow(1 - 2 / (N + 1), i)
        p = tab[i]
        licznik = licznik + p * alpha
        mianownik = mianownik + alpha
    return licznik / mianownik


def Symulacja(MACD, SIGNAL, dane):  # Algorytm inwestujący
    stan_konta = 1000
    ilosc_akcji = 0
    start = False
    koniec = True
    for i in range(0, 964):
        if MACD[i] > SIGNAL[i] and MACD[i + 1] <= SIGNAL[i + 1]:  # MACD przecina od dołu - kupno
            while stan_konta - dane[i] >= 0:
                stan_konta = stan_konta - dane[i]
                ilosc_akcji = ilosc_akcji + 1
            start = True  # trzeba najpierw akcję kupić, aby je sprzedać
            koniec = False
        elif MACD[i] < SIGNAL[i] and MACD[i + 1] >= SIGNAL[i + 1] and start:  # MACD przecina od góry - sprzedaż
            for j in range(0, ilosc_akcji):
                stan_konta = stan_konta + dane[i]
            koniec = True
            ilosc_akcji = 0
    if not koniec:
        for j in range(0, ilosc_akcji):
            stan_konta = stan_konta + dane[i]  # sprzedanie na koniec wszystkich aktualnie posiadanych akcji
    return stan_konta


def Analiza(dane):
    MACD = []
    for i in range(26, 1000):
        MACD.append(EMA(12, i, dane) - EMA(26, i, dane))

    SIGNAL = []
    for i in range(9, 974):
        SIGNAL.append(EMA(9, i, MACD))

    Roznica = []
    for i in range(0, 965, 5):
        Roznica.append(abs(MACD[i] - SIGNAL[i]))

    for i in range(1, 10):  # wyrównanie ilośc próbek MACD i SIGNAL w celu wyświetlenia ich na wykresach
        MACD.pop()

    Osx1 = np.arange(0, 965)
    OsX2 = np.arange(0, 965, 5)
    OsX3 = np.arange(0, 1000)

    print(Symulacja(MACD, SIGNAL, dane))

    fig1 = plt.figure(num=None, figsize=(23, 10), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(Osx1, MACD, '-', color='blue')
    plt.plot(Osx1, SIGNAL, '-', color='red')
    idx = np.argwhere(
        np.diff(np.sign(np.array(MACD) - np.array(SIGNAL)))).flatten()  # wyznaczenie miejsc przecięcia krzywych
    plt.plot(np.array(Osx1)[idx], np.array(MACD)[idx], '.', color='green')  # narysowanie ich na wykresie
    plt.ylabel("wartość MACD")
    plt.title("MACD / SIGNAL")
    plt.legend(['MACD', 'SIGNAL'])
    fig1.savefig('MACD.png', dpi=fig1.dpi)

    fig2 = plt.figure(num=None, figsize=(23, 10), dpi=80, facecolor='w', edgecolor='k')
    plt.bar(OsX2, Roznica, color='maroon', width=4)
    plt.ylabel("MACD-SIGNAL")
    plt.title("Różnica między MACD, a SIGNAL")
    plt.legend(['MACD-SIGNAL'])
    fig2.savefig('Roznica.png', dpi=fig2.dpi)

    fig3 = plt.figure(num=None, figsize=(23, 10), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(OsX3, dane, '-', color='green')
    plt.ylabel("Wartość")
    plt.title("Dane wejściowe")
    fig3.savefig('Dane wejsciowe.png', dpi=fig1.dpi)
    plt.show()


print("1 - USD/PLN")
print("2 - CDPROJEKT")
print("3 - MBANK")
opt = int(input())
dane = []
if opt == 1:
    probki = pd.read_csv('USD_PLN.csv', sep=',', usecols=['otwarcie'])
    for i in range(1000, 0, -1):
        tmp = probki['otwarcie'][i]
        dane.append(float(tmp.replace(",", ".")))  # obróbka formatu danych wejściowych
elif opt == 2:
    probki = pd.read_csv('CDPROJEKT.csv', sep=';', usecols=['otwarcie'])
    for i in range(1000, 0, -1):
        dane.append(float(probki['otwarcie'][i]))
elif opt == 3:
    probki = pd.read_csv('MBANK.csv', sep=',', usecols=['otwarcie'])
    for i in range(1000, 0, -1):
        tmp = probki['otwarcie'][i]
        dane.append(float(tmp.replace(",", ".")))  # obróbka formatu danych wejściowych

Analiza(dane)
