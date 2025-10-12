import pandas as pd
import matplotlib.pyplot as plt
import os
from gen_dane import odpowiedz_skokowa as odp_skok
import numpy as np

def main():
    # --- USTAWIENIA ---
    # Ustaw na True, aby podmienić sygnał 'out' z pliku na symulację.
    # Ustaw na False, aby analizować oryginalne dane z pliku CSV.
    PODMIEN_WYJSCIE_NA_SYMULACJE = True

    # Parametry do symulacji (używane tylko, gdy flaga powyżej to True)
    K_symulacji = 0.74  # Docelowe wzmocnienie obiektu
    T_symulacji = 0.00079  # Docelowa stała czasowa obiektu
    szum_symulacji = 0.004 # Poziom szumu do dodania do symulacji

    # --- KROK 1: WCZYTYWANIE DANYCH Z PLIKU CSV ---
    filename = "fol_gabriel/NewFile1.csv"
    filepath = os.path.join(os.path.dirname(__file__), filename)

    if not os.path.exists(filepath):
        print(f"❌ Plik '{filename}' nie został znaleziony w folderze:\n{os.path.dirname(__file__)}")
        return

    with open(filepath, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()

    if first_line.startswith("CH1"):
        print("📗 Wykryto format typu 2 (CH1, CH2, Start, Increment)")
        df = pd.read_csv(filepath, skiprows=2, header=None)
        parts = second_line.split(",")
        start = float(parts[2])
        increment = float(parts[3])
        df = df.iloc[:, :2]
        df.columns = ["in", "out"]
        df["czas"] = [start + i * increment for i in range(len(df))]
        df = df[["czas", "in", "out"]]
        df = df.dropna()
        df = df.astype(float)
    else:
        print("📘 Wykryto format typu 1 (czas, sygnał_1, sygnał_2)")
        df = pd.read_csv(filepath, skiprows=2, header=None)
        df.columns = ["czas", "in", "out"]
        df = df.dropna()
        df = df.astype(float)

    # --- KROK 2: PRZESUNIĘCIE SYGNAŁÓW DO ZERA ---
    in_offset = -df["in"].iloc[0]
    out_offset = -df["out"].iloc[0]
    df["in"] = df["in"] + in_offset
    df["out"] = df["out"] + out_offset
    
    # --- KROK 3: ZEROWANIE SYGNAŁU WYJŚCIOWEGO PRZED CZASEM t=0 (NOWA SEKCJA) ---
    print("🧹 Zerowanie sygnału 'out' dla czasu <= 0...")
    df.loc[df['czas'] <= 0, 'out'] = 0
    
    # --- KROK 4: PODMIANA SYGNAŁU WYJŚCIOWEGO NA SYMULACJĘ ---
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        print("\n⚙️  Wykonywanie podmiany sygnału wyjściowego na symulację...")

        A = df["in"].iloc[-1]
        t = df["czas"].to_numpy()

        # Generujemy odpowiedź tylko dla t > 0
        y_symulowane = np.zeros_like(t)
        maska_dodatnia = t > 0
        t_dodatnie = t[maska_dodatnia]
        
        y_czyste = K_symulacji * A * (1 - np.exp(-t_dodatnie / T_symulacji))
        np.random.seed(42)
        noise = np.random.normal(0, szum_symulacji * K_symulacji * A, size=len(y_czyste))
        
        y_symulowane[maska_dodatnia] = y_czyste + noise

        df["out"] = y_symulowane
        print(f"✅  Kolumna 'out' została zastąpiona symulacją z K={K_symulacji} i T={T_symulacji}s.\n")


    # --- KROK 5: OBLICZANIE PARAMETRÓW ODPOWIEDZI SKOKOWEJ ---
    kp = df["out"].iloc[-1] / df["in"].iloc[-1]
    print(f"Wyznaczone wzmocnienie K_p = {kp:.4f}")
    
    wartosc_koncowa = kp * df['in'].iloc[-1]
    wartosc_w_punkcie_T = wartosc_koncowa * 0.6321
    
    df['roznica'] = abs(df['out'] - wartosc_w_punkcie_T)
    idx = df['roznica'].idxmin()
    Tp = df.loc[idx, 'czas']
    print(f'Wyznaczona stała czasowa T_p = {Tp:.4f}s')

    # --- KROK 6: RYSOWANIE ---
    plt.figure(figsize=(10, 6))
    plt.plot(df["czas"], df["in"], label="Wejście (IN) z pliku CSV")
    
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        plt.plot(df["czas"], df["out"], label=f"Wyjście (Symulowane, K={K_symulacji}, T={T_symulacji}s)")
        plt.title("Sygnały IN (z pliku) / OUT (symulowane)")
    else:
        plt.plot(df["czas"], df["out"], label="Wyjście (OUT) z pliku CSV")
        plt.title("Sygnały IN/OUT z pliku CSV")

    plt.xlabel("Czas [s]")
    plt.ylabel("Amplituda [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()