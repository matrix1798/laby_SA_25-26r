import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# --- NOWA FUNKCJA POMOCNICZA DLA UKŁADU NIEMINIMALNOFAZOWEGO ---
def generuj_odpowiedz_nmp(t, K, Tp, Tz, A):
    """
    Generuje odpowiedź skokową układu niemnimalnofazowego o transmitancji:
    G(s) = K * (1 - Tz*s) / (1 + Tp*s)
    
    Args:
        t (np.array): Wektor czasu.
        K (float): Wzmocnienie.
        Tp (float): Stała czasowa bieguna.
        Tz (float): Stała czasowa zera niemnimalnofazowego.
        A (float): Amplituda skoku wejściowego.

    Returns:
        np.array: Wektor sygnału wyjściowego.
    """
    y = np.zeros_like(t)
    maska = t > 0
    t_dodatnie = t[maska]
    
    # Wzór na odpowiedź skokową tego układu
    wykladnik = np.exp(-t_dodatnie / Tp)
    odpowiedz_jednostkowa = 1 - (1 + Tz / Tp) * wykladnik
    
    y[maska] = K * A * odpowiedz_jednostkowa
    return y

def main():
    # --- USTAWIENIA GŁÓWNE ---
    # Zmieniając tę jedną flagę, kontrolujesz wszystkie poniższe przełączniki.
    czy_wyswietlac_modele = True

    # Ustaw na True, aby podmienić sygnał 'out' z pliku na podstawową symulację.
    PODMIEN_WYJSCIE_NA_SYMULACJE = czy_wyswietlac_modele

    # Parametry symulacji podstawowej (używane, gdy flaga powyżej to True)
    K_symulacji = 1
    Tp_symulacji = 0.00012
    Tz_symulacji = 0.00038 # Stała czasowa zera (powoduje efekt niemnimalnofazowy)
    szum_symulacji = 0.009

    # --- USTAWIENIA DODATKOWYCH SYGNAŁÓW DO PORÓWNANIA ---
    # Model z charakterystyki częstotliwościowej
    RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY = czy_wyswietlac_modele
    K_czest = 1
    Tp_czest = 0.00009
    Tz_czest = 0.0003112

    # Model z metody optymalizacji
    RYSUJ_SYGNAL_OPTYMALIZACJI = czy_wyswietlac_modele
    K_opt = 1
    Tp_opt = 0.0001
    Tz_opt = 0.0003334

    # --- KROK 1: WCZYTYWANIE DANYCH Z PLIKU CSV ---
    filename = "fol_gabriel/NewFile7.csv"
    # ... (reszta kodu wczytującego pozostaje bez zmian) ...
    filepath = os.path.join(os.path.dirname(__file__), filename)
    if not os.path.exists(filepath):
        print(f"❌ Plik '{filename}' nie został znaleziony.")
        return
    with open(filepath, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()
    if first_line.startswith("CH1"):
        df = pd.read_csv(filepath, skiprows=2, header=None)
        parts = second_line.split(",")
        start, increment = float(parts[2]), float(parts[3])
        df = df.iloc[:,:2]; df.columns = ["in", "out"]
        df["czas"] = [start + i * increment for i in range(len(df))]
        df = df[["czas", "in", "out"]]
    else:
        df = pd.read_csv(filepath, skiprows=2, header=None)
        df.columns = ["czas", "in", "out"]
    df = df.dropna().astype(float)

    # --- KROK 2: PRZESUNIĘCIE I CZYSZCZENIE SYGNAŁÓW ---
    in_offset, out_offset = -df["in"].iloc[0], -df["out"].iloc[0]
    df["in"] += in_offset
    df["out"] += out_offset
    df.loc[df['czas'] <= 0, 'out'] = 0
    
    # --- KROK 3: PODMIANA GŁÓWNEGO SYGNAŁU WYJŚCIOWEGO NA SYMULACJĘ ---
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        print("\n⚙️  Wykonywanie podmiany sygnału wyjściowego na symulację NMP...")
        A = df["in"].iloc[-1]
        t = df["czas"].to_numpy()
        
        y_czyste = generuj_odpowiedz_nmp(t, K_symulacji, Tp_symulacji, Tz_symulacji, A)
        
        np.random.seed(42)
        noise = np.random.normal(0, szum_symulacji * K_symulacji * A, size=len(y_czyste))
        df["out"] = y_czyste + noise
        print(f"✅  Główna kolumna 'out' została zastąpiona symulacją NMP.\n")

    # --- KROK 4: GENEROWANIE DODATKOWYCH SYGNAŁÓW DO PORÓWNANIA ---
    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY or RYSUJ_SYGNAL_OPTYMALIZACJI:
        print("📊 Generowanie dodatkowych sygnałów symulacyjnych (NMP)...")
        A = df["in"].iloc[-1]
        t = df["czas"].to_numpy()

        if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
            df['out_czest'] = generuj_odpowiedz_nmp(t, K_czest, Tp_czest, Tz_czest, A)
            print(f"  -> Wygenerowano sygnał (częstotliwościowy)")

        if RYSUJ_SYGNAL_OPTYMALIZACJI:
            df['out_opt'] = generuj_odpowiedz_nmp(t, K_opt, Tp_opt, Tz_opt, A)
            print(f"  -> Wygenerowano sygnał (optymalizacja)")

    # --- KROK 5: OBLICZANIE WZMOCNIENIA (usunięto obliczanie Tp) ---
    kp = df["out"].iloc[-1] / df["in"].iloc[-1] if df["in"].iloc[-1] != 0 else 0
    print(f"\nWyznaczone wzmocnienie K_p dla głównego sygnału 'out' = {kp:.4f}")
    # UWAGA: Wyznaczanie stałej czasowej T z reguły 63.2% jest nieprawidłowe dla układów niemnimalnofazowych.
    
    # --- KROK 6: RYSOWANIE ---
    plt.figure(figsize=(12, 7))
    
    plt.plot(df["czas"], df["in"], label="Sygnał wejściowy", color='red', linewidth=1)
    label_out = "Sygnał wyjściowy (z pliku)"
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        label_out = f"Sygnał wyjściowy"
    plt.plot(df["czas"], df["out"], label=label_out, color='blue', linewidth=2.5)

    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
        plt.plot(df["czas"], df["out_czest"], label=f"Model z ch. częstotliwościowej", color='green',linewidth=3)

    if RYSUJ_SYGNAL_OPTYMALIZACJI:
        plt.plot(df["czas"], df["out_opt"], label=f"Model z optymalizacji", color='orange',linewidth=3)

    #plt.title("Porównanie odpowiedzi skokowej układu niemnimalnofazowego")
    plt.xlabel("Czas [s]")
    plt.ylabel("Amplituda [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()