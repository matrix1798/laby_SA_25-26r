import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


# --- Funkcja pomocnicza do wczytania i przygotowania danych ---
def wczytaj_csv(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()

    # --- Dla plików z nagłówkiem typu CH1,CH2,Start,Increment ---
    if first_line.startswith("CH1"):
        parts = second_line.split(",")
        start = float(parts[2])
        increment = float(parts[3])
        df = pd.read_csv(filepath, skiprows=2, header=None)
        df = df.iloc[:, :2]
        df.columns = ["in", "out"]
        df["czas"] = [start + i * increment for i in range(len(df))]
        df = df[["czas", "in", "out"]]
    else:
        # Dla zwykłego CSV z kolumnami czas, in, out
        df = pd.read_csv(filepath, skiprows=2, header=None)
        df.columns = ["czas", "in", "out"]

    # --- Oczyszczenie i konwersja ---
    df = df.dropna().astype(float)

    # --- Normalizacja wejścia i wyjścia względem początku pomiaru ---
    df["in"] = df["in"] - df["in"].iloc[0]
    df["out"] = df["out"] - df["out"].iloc[0]

    # --- 🔧 Kluczowa poprawka ---
    # Wyrównaj wyjście tak, by w chwili t≈0 zaczynało od 0 V
    idx_0 = (df["czas"] - 0).abs().idxmin()
    wartosc_0 = df.loc[idx_0, "out"]
    df["out"] = df["out"] - wartosc_0

    # Wyzeruj wszystkie próbki przed czasem 0 s
    df.loc[df["czas"] <= 0, "out"] = 0

    return df


# --- Funkcja do obliczenia Ti (T0) ---
def oblicz_Ti(df):
    delta_Uwe = df["in"].max() - df["in"].min()
    delta_Uwy = df["out"].iloc[-1] - df["out"].iloc[0]
    delta_t = df["czas"].iloc[-1] - df["czas"].iloc[0]

    if delta_Uwy == 0:
        return np.nan

    Ti = (delta_Uwe * delta_t) / delta_Uwy
    return Ti


# --- Główna część programu ---
def main():
    folder = os.path.dirname(__file__) if "__file__" in globals() else os.getcwd()
    pliki = [
        "fol_gabriel/NewFile4.csv",
        "fol_gabriel/NewFile5.csv",
        "fol_gabriel/NewFile6.csv",
    ]

    # Kolory
    kolory_we = ["red", "orange", "gold"]
    kolory_wy = ["blue", "green", "purple"]

    # Opisy napięć wejściowych
    napiecia_we = [2, 4, 8]

    dane = []
    Ti_wartosci = []

    plt.figure(figsize=(12, 7))

    for i, plik in enumerate(pliki, start=1):
        sciezka = os.path.join(folder, plik)
        if not os.path.exists(sciezka):
            print(f"❌ Brak pliku: {plik}")
            continue

        df = wczytaj_csv(sciezka)
        dane.append(df)

        # Rysowanie wejść i wyjść
        plt.plot(
            df["czas"],
            df["in"],
            color=kolory_we[i - 1],
            linewidth=2,
            label=f"Sygnał wejściowy {i} [{napiecia_we[i - 1]}V]",
        )
        plt.plot(
            df["czas"],
            df["out"],
            color=kolory_wy[i - 1],
            linewidth=2.5,
            label=f"Sygnał wyjściowy {i}",
        )

        # Oblicz Ti
        Ti = oblicz_Ti(df)
        Ti_wartosci.append(Ti)
        print(f"📘 {plik}: Ti = {Ti:.6f} s")

    # Oblicz średnią z Ti
    Ti_srednie = np.nanmean(Ti_wartosci)
    print("\n====================================")
    print(f"📊 Średnia wartość Ti = {Ti_srednie:.6f} s")
    print("====================================")

    # --- Rysowanie ---
    plt.title("Odpowiedzi układu całkującego dla trzech napięć wejściowych")
    plt.xlabel("Czas [s]")
    plt.ylabel("Napięcie [V]")
    plt.grid(True)
    plt.axhline(0, color="black", linestyle="--", linewidth=1)  # poziom 0 V
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
