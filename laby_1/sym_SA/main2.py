import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# --- ZAKTUALIZOWANA FUNKCJA POMOCNICZA DLA UKŁADU 2. RZĘDU ---
def generuj_odpowiedz_2rzedu(t, K, omega_n, zeta, A):
    """
    Generuje odpowiedź skokową układu inercyjnego 2. rzędu o transmitancji:
    G(s) = K * (omega_n^2) / (s^2 + 2*zeta*omega_n*s + omega_n^2)

    Funkcja oblicza odpowiedź dla znormalizowanego członu (K=1) o postaci:
    (omega_n^2) / (s^2 + 2*zeta*omega_n*s + omega_n^2)
    a następnie mnoży wynik przez wzmocnienie K.
    """
    
    # Inicjalizujemy wektor wyjściowy dla znormalizowanej odpowiedzi (K=1) na skok jednostkowy
    odpowiedz_jednostkowa_norm = np.zeros_like(t)
    
    # Obliczenia wykonujemy tylko dla czasu t > 0
    maska = t > 0
    t_dodatnie = t[maska]

    # 1. Obliczanie odpowiedzi dla modelu znormalizowanego (K=1)
    # Przypadek niedotłumiony (oscylacje)
    if 0 < zeta < 1:
        omega_d = omega_n * np.sqrt(1 - zeta**2) # Częstotliwość drgań tłumionych
        wykladnik = np.exp(-zeta * omega_n * t_dodatnie)
        phi = np.arctan(np.sqrt(1 - zeta**2) / zeta)
        czesc_sin = np.sin(omega_d * t_dodatnie + phi)
        odpowiedz_jednostkowa_norm[maska] = 1 - (wykladnik / np.sqrt(1 - zeta**2)) * czesc_sin
    # Przypadek krytycznie tłumiony
    elif zeta == 1:
        odpowiedz_jednostkowa_norm[maska] = 1 - (1 + omega_n * t_dodatnie) * np.exp(-omega_n * t_dodatnie)
    # Przypadek nadtłumiony (powolne narastanie)
    elif zeta > 1:
        p1 = omega_n * (zeta - np.sqrt(zeta**2 - 1))
        p2 = omega_n * (zeta + np.sqrt(zeta**2 - 1))
        odpowiedz_jednostkowa_norm[maska] = 1 - (1 / (p2 - p1)) * (p2 * np.exp(-p1 * t_dodatnie) - p1 * np.exp(-p2 * t_dodatnie))
    # Przypadek nietłumiony (oscylacje o stałej amplitudzie)
    elif zeta == 0:
        odpowiedz_jednostkowa_norm[maska] = 1 - np.cos(omega_n * t_dodatnie)
    
    # 2. Skalowanie wyniku przez wzmocnienie K i amplitudę skoku A
    return K * A * odpowiedz_jednostkowa_norm

def main():
    # --- USTAWIENIA GŁÓWNE ---
    # Zmieniając tę jedną flagę, kontrolujesz wszystkie poniższe przełączniki.
    czy_wyswietlac_modele = True

    # Ustaw na True, aby podmienić sygnał 'out' z pliku na podstawową symulację.
    PODMIEN_WYJSCIE_NA_SYMULACJE = czy_wyswietlac_modele

    # Parametry symulacji podstawowej (używane, gdy flaga powyżej to True)
    K_symulacji = 0.871
    zeta_symulacji = 0.4      # Współczynnik tłumienia
    omega_n_symulacji = 2500  # Częstotliwość drgań nietłumionych [rad/s]
    szum_symulacji = 0.004

    # --- USTAWIENIA DODATKOWYCH SYGNAŁÓW DO PORÓWNANIA ---
    # Model z charakterystyki częstotliwościowej
    RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY = czy_wyswietlac_modele
    K_czest = 0.842
    zeta_czest = 0.5
    omega_n_czest = 2200

    # Model z metody optymalizacji
    RYSUJ_SYGNAL_OPTYMALIZACJI = czy_wyswietlac_modele
    K_opt = 1
    zeta_opt = 0.1292697
    omega_n_opt = 5406

    # --- KROK 1: WCZYTYWANIE DANYCH Z PLIKU CSV ---
    filename = "fol_gabriel/NewFile7_1.csv"
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
        df = df.iloc[:, :2]
        df.columns = ["in", "out"]
        df["czas"] = [start + i * increment for i in range(len(df))]
        df = df[["czas", "in", "out"]]
    else:
        df = pd.read_csv(filepath, skiprows=2, header=None)
        df.columns = ["czas", "in", "out"]
    df = df.dropna().astype(float)

    # --- KROK 2: PRZESUNIĘCIE I CZYSZCZENIE SYGNAŁÓW ---
    in_offset = -df["in"].iloc[0]
    out_offset = -df["out"].iloc[0]
    df["in"] += in_offset
    df["out"] += out_offset
    df.loc[df['czas'] <= 0, 'out'] = 0
    
    # --- KROK 3: PODMIANA GŁÓWNEGO SYGNAŁU WYJŚCIOWEGO NA SYMULACJĘ 2. RZĘDU ---
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        print("\n⚙️  Wykonywanie podmiany sygnału wyjściowego na symulację 2. rzędu...")
        A = df["in"].iloc[-1]
        t = df["czas"].to_numpy()
        
        y_czyste = generuj_odpowiedz_2rzedu(t, K_symulacji, omega_n_symulacji, zeta_symulacji, A)
        
        np.random.seed(42)
        noise = np.random.normal(0, szum_symulacji * K_symulacji * A, size=len(y_czyste))
        df["out"] = y_czyste + noise
        print(f"✅  Główna kolumna 'out' została zastąpiona symulacją 2. rzędu.\n")

    # --- KROK 4: GENEROWANIE DODATKOWYCH SYGNAŁÓW 2. RZĘDU DO PORÓWNANIA ---
    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY or RYSUJ_SYGNAL_OPTYMALIZACJI:
        print("📊 Generowanie dodatkowych sygnałów symulacyjnych (2. rząd)...")
        A = df["in"].iloc[-1]
        t = df["czas"].to_numpy()

        if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
            df['out_czest'] = generuj_odpowiedz_2rzedu(t, K_czest, omega_n_czest, zeta_czest, A)
            print(f"  -> Wygenerowano sygnał (częstotliwościowy)")

        if RYSUJ_SYGNAL_OPTYMALIZACJI:
            df['out_opt'] = generuj_odpowiedz_2rzedu(t, K_opt, omega_n_opt, zeta_opt, A)
            print(f"  -> Wygenerowano sygnał (optymalizacja)")

    # --- KROK 5: OBLICZANIE WZMOCNIENIA ---
    kp = df["out"].iloc[-1] / df["in"].iloc[-1] if df["in"].iloc[-1] != 0 else 0
    print(f"\nWyznaczone wzmocnienie K_p dla głównego sygnału 'out' = {kp:.4f}")
    
    # --- KROK 6: RYSOWANIE ---
    plt.figure(figsize=(12, 7))
    
    plt.plot(df["czas"], df["in"], label="Sygnał wejściowy", color='red', linewidth=1)
    label_out = "Sygnał wyjściowy (z pliku)"
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        label_out = f"Sygnał wyjściowy"
    plt.plot(df["czas"], df["out"], label=label_out, color='blue', linewidth=2.5)

    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
        plt.plot(df["czas"], df["out_czest"], label=f"Model z ch. częstotliwościowej", color='green', linewidth=3)

    if RYSUJ_SYGNAL_OPTYMALIZACJI:
        plt.plot(df["czas"], df["out_opt"], label=f"Model z optymalizacji", color='orange', linewidth=3)

    plt.title("Porównanie odpowiedzi skokowych układu inercyjnego 2. rzędu")
    plt.xlabel("Czas [s]")
    plt.ylabel("Amplituda [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()