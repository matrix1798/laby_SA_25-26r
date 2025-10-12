import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def odpowiedz_skokowa(K=1.0, T=0.1, t_end=1.0, dt=0.0005, szum=0.0002):
    """
    Generuje dane odpowiedzi skokowej układu inercyjnego I rzędu z szumem.

    Parametry:
    -----------
    K : float
        Wzmocnienie układu
    T : float
        Stała czasowa [s]
    t_end : float
        Czas symulacji [s]
    dt : float
        Krok próbkowania [s]
    szum : float
        Poziom szumu (0.02 = 2% amplitudy)

    Zwraca:
    --------
    pandas.DataFrame z kolumnami ['czas', 'in', 'out']
    """

    # Wektor czasu
    t = np.arange(0, t_end, dt)

    # Sygnał wejściowy – skok jednostkowy
    u = np.ones_like(t)

    # Idealna odpowiedź skokowa układu 1. rzędu
    y_czyste = K * (1 - np.exp(-t / T))

    # Dodaj losowy szum (realistyczny)
    np.random.seed(42)  # żeby wynik był powtarzalny
    noise = np.random.normal(0, szum * K, size=len(y_czyste))
    y_szum = y_czyste + noise

    # Tworzymy DataFrame
    df = pd.DataFrame({
        "czas": t,
        "in": u * K,   # sygnał wejściowy (skok)
        "out": y_szum  # odpowiedź z szumem
    })

    return df


def main():
    # Parametry układu
    K = 2.0     # wzmocnienie
    T = 0.15    # stała czasowa [s]
    t_end = 1.0 # czas symulacji
    dt = 0.0005 # krok próbkowania
    szum = 0.003 # 3% szumu

    df = odpowiedz_skokowa(K, T, t_end, dt, szum)

    # Zapis do pliku CSV (w tym samym folderze)
    filename = os.path.join(os.path.dirname(__file__), "symulacja.csv")
    df.to_csv(filename, index=False, header=False)
    print(f"✅ Dane zapisane do: {filename}")

    # Wizualizacja
    plt.figure(figsize=(10,6))
    plt.plot(df["czas"], df["in"], label="Wejście (IN) – skok")
    plt.plot(df["czas"], df["out"], label="Wyjście (OUT) – z szumem")
    plt.title(f"Odpowiedź skokowa układu 1. rzędu (K={K}, T={T}s, szum={szum*100:.0f}%)")
    plt.xlabel("Czas [s]")
    plt.ylabel("Amplituda [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
