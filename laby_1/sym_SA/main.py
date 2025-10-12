import pandas as pd
import matplotlib.pyplot as plt
import os
from gen_dane import odpowiedz_skokowa as odp_skok
import numpy as np

def main():
    # --- USTAWIENIA GŁÓWNE ---
    # Ustaw na True, aby podmienić sygnał 'out' z pliku na podstawową symulację.
    PODMIEN_WYJSCIE_NA_SYMULACJE = True

    # Parametry symulacji podstawowej (używane, gdy flaga powyżej to True)
    K_symulacji = 0.74
    T_symulacji = 0.00079
    szum_symulacji = 0.004

    # --- USTAWIENIA DODATKOWYCH SYGNAŁÓW DO PORÓWNANIA ---
    # Włącz/wyłącz rysowanie sygnału z modelu uzyskanego z charakterystyki częstotliwościowej
    RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY = True
    K_czest = 0.72  # Wzmocnienie z analizy częstotliwościowej
    T_czest = 0.0008  # Stała czasowa z analizy częstotliwościowej

    # Włącz/wyłącz rysowanie sygnału z modelu uzyskanego metodą optymalizacji
    RYSUJ_SYGNAL_OPTYMALIZACJI = True
    K_opt = 0.75    # Wzmocnienie z optymalizacji
    T_opt = 0.0007    # Stała czasowa z optymalizacji

    # --- KROK 1: WCZYTYWANIE DANYCH Z PLIKU CSV ---
    filename = "fol_gabriel/NewFile1.csv"
    filepath = os.path.join(os.path.dirname(__file__), filename)

    if not os.path.exists(filepath):
        print(f"❌ Plik '{filename}' nie został znaleziony.")
        return

    # ... (reszta kodu wczytującego plik pozostaje bez zmian) ...
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
    else:
        print("📘 Wykryto format typu 1 (czas, sygnał_1, sygnał_2)")
        df = pd.read_csv(filepath, skiprows=2, header=None)
        df.columns = ["czas", "in", "out"]

    df = df.dropna().astype(float)

    # --- KROK 2: PRZESUNIĘCIE I CZYSZCZENIE SYGNAŁÓW ---
    in_offset = -df["in"].iloc[0]
    out_offset = -df["out"].iloc[0]
    df["in"] = df["in"] + in_offset
    df["out"] = df["out"] + out_offset
    
    print("🧹 Zerowanie sygnału 'out' dla czasu <= 0...")
    df.loc[df['czas'] <= 0, 'out'] = 0
    
    # --- KROK 3: PODMIANA GŁÓWNEGO SYGNAŁU WYJŚCIOWEGO NA SYMULACJĘ ---
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        # ... (ta sekcja pozostaje bez zmian) ...
        print("\n⚙️  Wykonywanie podmiany sygnału wyjściowego na symulację...")
        A = df["in"].iloc[-1]
        t = df["czas"].to_numpy()
        y_symulowane = np.zeros_like(t)
        maska_dodatnia = t > 0
        t_dodatnie = t[maska_dodatnia]
        y_czyste = K_symulacji * A * (1 - np.exp(-t_dodatnie / T_symulacji))
        np.random.seed(42)
        noise = np.random.normal(0, szum_symulacji * K_symulacji * A, size=len(y_czyste))
        y_symulowane[maska_dodatnia] = y_czyste + noise
        df["out"] = y_symulowane
        print(f"✅  Główna kolumna 'out' została zastąpiona symulacją z K={K_symulacji} i T={T_symulacji}s.\n")

    # --- KROK 4: GENEROWANIE DODATKOWYCH SYGNAŁÓW DO PORÓWNANIA (NOWA SEKCJA) ---
    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY or RYSUJ_SYGNAL_OPTYMALIZACJI:
        print("📊 Generowanie dodatkowych sygnałów symulacyjnych do porównania...")
        A = df["in"].iloc[-1]
        t = df["czas"].to_numpy()
        maska = t > 0
        t_dodatnie = t[maska]

        if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
            y_sym_czest = np.zeros_like(t)
            y_czyste_czest = K_czest * A * (1 - np.exp(-t_dodatnie / T_czest))
            y_sym_czest[maska] = y_czyste_czest
            df['out_czest'] = y_sym_czest
            print(f"  -> Wygenerowano sygnał dla K={K_czest}, T={T_czest}s (częstotliwościowy)")

        if RYSUJ_SYGNAL_OPTYMALIZACJI:
            y_sym_opt = np.zeros_like(t)
            y_czyste_opt = K_opt * A * (1 - np.exp(-t_dodatnie / T_opt))
            y_sym_opt[maska] = y_czyste_opt
            df['out_opt'] = y_sym_opt
            print(f"  -> Wygenerowano sygnał dla K={K_opt}, T={T_opt}s (optymalizacja)")

    # --- KROK 5: OBLICZANIE PARAMETRÓW DLA GŁÓWNEGO SYGNAŁU ---
    # Analiza jest zawsze przeprowadzana na sygnale z kolumny 'out'
    kp = df["out"].iloc[-1] / df["in"].iloc[-1] if df["in"].iloc[-1] != 0 else 0
    print(f"\nWyznaczone wzmocnienie K_p dla głównego sygnału 'out' = {kp:.4f}")
    
    wartosc_koncowa = kp * df['in'].iloc[-1]
    wartosc_w_punkcie_T = wartosc_koncowa * 0.6321
    
    df['roznica'] = abs(df['out'] - wartosc_w_punkcie_T)
    idx = df['roznica'].idxmin()
    Tp = df.loc[idx, 'czas']
    print(f"Wyznaczona stała czasowa T_p dla głównego sygnału 'out' = {Tp:.4f}s")

    # --- KROK 6: RYSOWANIE ---
    plt.figure(figsize=(12, 7))
    
    # Rysowanie sygnałów głównych
    plt.plot(df["czas"], df["in"], label="Wejście (IN) z pliku", color='black', linewidth=2)
    label_out = "Wyjście (z pliku CSV)" if not PODMIEN_WYJSCIE_NA_SYMULACJE else f"Wyjście (Symulacja bazowa K={K_symulacji}, T={T_symulacji})"
    plt.plot(df["czas"], df["out"], label=label_out, color='blue', linewidth=2.5)

    # Rysowanie sygnałów porównawczych
    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
        plt.plot(df["czas"], df["out_czest"], label=f"Model (częstotliwościowy, K={K_czest}, T={T_czest}s)", linestyle='--', color='green')

    if RYSUJ_SYGNAL_OPTYMALIZACJI:
        plt.plot(df["czas"], df["out_opt"], label=f"Model (optymalizacja, K={K_opt}, T={T_opt}s)", linestyle=':', color='red')

    plt.title("Porównanie odpowiedzi skokowej obiektu z modelami")
    plt.xlabel("Czas [s]")
    plt.ylabel("Amplituda [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()