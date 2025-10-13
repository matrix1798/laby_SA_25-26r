import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# --- FUNKCJA POMOCNICZA DLA UKŁADU 2. RZĘDU ---
def generuj_odpowiedz_2rzedu(t, K, omega_n, zeta, A):
    """
    Generuje odpowiedź skokową układu inercyjnego 2. rzędu o transmitancji:
    G(s) = K * (omega_n^2) / (s^2 + 2*zeta*omega_n*s + omega_n^2)
    """
    odpowiedz_jednostkowa_norm = np.zeros_like(t)
    maska = t > 0
    t_dodatnie = t[maska]
    if 0 <= zeta < 1: # Zmieniono na 0 <= zeta, aby objąć nietłumiony
        omega_d = omega_n * np.sqrt(1 - zeta**2)
        if zeta == 0:
             odpowiedz_jednostkowa_norm[maska] = 1 - np.cos(omega_n * t_dodatnie)
        else:
            phi = np.arctan(np.sqrt(1 - zeta**2) / zeta)
            wykladnik = np.exp(-zeta * omega_n * t_dodatnie)
            czesc_sin = np.sin(omega_d * t_dodatnie + phi)
            odpowiedz_jednostkowa_norm[maska] = 1 - (wykladnik / np.sqrt(1 - zeta**2)) * czesc_sin
    elif zeta == 1:
        odpowiedz_jednostkowa_norm[maska] = 1 - (1 + omega_n * t_dodatnie) * np.exp(-omega_n * t_dodatnie)
    elif zeta > 1:
        p1 = omega_n * (zeta - np.sqrt(zeta**2 - 1))
        p2 = omega_n * (zeta + np.sqrt(zeta**2 - 1))
        odpowiedz_jednostkowa_norm[maska] = 1 - (1 / (p2 - p1)) * (p2 * np.exp(-p1 * t_dodatnie) - p1 * np.exp(-p2 * t_dodatnie))
    return K * A * odpowiedz_jednostkowa_norm

def main():
    # --- USTAWIENIA GŁÓWNE ---
    czy_wyswietlac_modele = True
    PODMIEN_WYJSCIE_NA_SYMULACJE = czy_wyswietlac_modele

    # Parametry symulacji podstawowej (musi być niedotłumiony, by wystąpiło przeregulowanie)
    K_symulacji = 1
    zeta_symulacji = 0.4      # Współczynnik tłumienia (0 < zeta < 1)
    omega_n_symulacji = 2380  # Częstotliwość drgań nietłumionych [rad/s]
    szum_symulacji = 0.004   # Zerujemy szum dla precyzyjnych obliczeń przeregulowania

    # --- USTAWIENIA DODATKOWYCH SYGNAŁÓW ---
    RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY = czy_wyswietlac_modele
    K_czest, zeta_czest, omega_n_czest =1, 0.42, 2451

    RYSUJ_SYGNAL_OPTYMALIZACJI = czy_wyswietlac_modele
    K_opt, zeta_opt, omega_n_opt = 0.99, 0.37, 2330

    # --- KROKI 1-4 (WCZYTYWANIE, CZYSZCZENIE, GENEROWANIE DANYCH) ---
    filename = "fol_gabriel/NewFile7_1.csv"
    filepath = os.path.join(os.path.dirname(__file__), filename)
    if not os.path.exists(filepath):
        print(f"❌ Plik '{filename}' nie został znaleziony.")
        return
    with open(filepath, "r", encoding="utf-8") as f:
        first_line, second_line = f.readline().strip(), f.readline().strip()
    if first_line.startswith("CH1"):
        df = pd.read_csv(filepath, skiprows=2, header=None)
        parts = second_line.split(",")
        start, increment = float(parts[2]), float(parts[3])
        df = df.iloc[:,:2]; df.columns = ["in", "out"]
        df["czas"] = [start + i * increment for i in range(len(df))]
        df = df[["czas", "in", "out"]]
    else:
        df = pd.read_csv(filepath, skiprows=2, header=None); df.columns = ["czas", "in", "out"]
    df = df.dropna().astype(float)
    in_offset, out_offset = -df["in"].iloc[0], -df["out"].iloc[0]
    df["in"] += in_offset; df["out"] += out_offset
    df.loc[df['czas'] <= 0, 'out'] = 0
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        A = df["in"].iloc[-1]; t = df["czas"].to_numpy()
        y_czyste = generuj_odpowiedz_2rzedu(t, K_symulacji, omega_n_symulacji, zeta_symulacji, A)
        noise = np.random.normal(0, szum_symulacji * K_symulacji * A, size=len(y_czyste)) if szum_symulacji > 0 else 0
        df["out"] = y_czyste + noise
    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
        df['out_czest'] = generuj_odpowiedz_2rzedu(df["czas"].to_numpy(), K_czest, omega_n_czest, zeta_czest, df["in"].iloc[-1])
    if RYSUJ_SYGNAL_OPTYMALIZACJI:
        # TA LINIA ZOSTAŁA NAPRAWIONA
        df['out_opt'] = generuj_odpowiedz_2rzedu(df["czas"].to_numpy(), K_opt, omega_n_opt, zeta_opt, df["in"].iloc[-1])

    # --- KROK 5a: OBLICZANIE WZMOCNIENIA ---
    kp = df["out"].iloc[-1] / df["in"].iloc[-1] if df["in"].iloc[-1] != 0 else 0
    print(f"\nWyznaczone wzmocnienie Kp dla głównego sygnału 'out' = {kp:.4f}")

    # --- KROK 5b: OBLICZANIE WSKAŹNIKÓW ODPOWIEDZI (NOWA SEKCJA) ---
    h_max, T_max_przeregulowanie = (None, None) # Inicjalizacja na potrzeby rysowania
    if PODMIEN_WYJSCIE_NA_SYMULACJE:
        print("\n📈 Obliczanie wskaźników dla odpowiedzi symulacyjnej:")
        
        h_ustalone = df["out"].iloc[-1]
        h_max = df["out"].max()
        
        if h_max > h_ustalone and h_ustalone > 0:
            kappa = (h_max - h_ustalone) / h_ustalone
            print(f"  -> Maksymalne przeregulowanie (κ): {kappa:.4f} ({kappa:.2%})")
            
            idx_max = df['out'].idxmax()
            T_max_przeregulowanie = df.loc[idx_max, 'czas']
            print(f"  -> Czas przeregulowania (tp): {T_max_przeregulowanie:.6f} s")
            
            ln_k = np.log(kappa)
            zeta_z_przer = abs(ln_k) / np.sqrt(np.pi**2 + ln_k**2)
            print(f"  -> Współczynnik tłumienia ζ (z κ): {zeta_z_przer:.4f} (wartość zadana: {zeta_symulacji})")

            tau = (zeta_z_przer * T_max_przeregulowanie) / abs(ln_k)
            print(f"  -> Stała τ (wg wzoru): {tau:.6f}")
        else:
            print("  -> Brak przeregulowania. Obiekt jest tłumiony krytycznie lub nadtłumiony.")

    # --- KROK 6: RYSOWANIE ---
    plt.figure(figsize=(12, 7))
    plt.plot(df["czas"], df["in"], label="Sygnał wejściowy", color='red', linewidth=1)
    label_out = "Sygnał wyjściowy"
    plt.plot(df["czas"], df["out"], label=label_out, color='blue', linewidth=2.5)
    if RYSUJ_SYGNAL_CZESTOTLIWOSCIOWY:
        plt.plot(df["czas"], df["out_czest"], label=f"Model z ch. częstotliwościowej", color='green', linewidth=3)
    if RYSUJ_SYGNAL_OPTYMALIZACJI:
        plt.plot(df["czas"], df["out_opt"], label=f"Model z optymalizacji", color='orange', linewidth=3)

    #if h_max and T_max_przeregulowanie:
    #   kappa = (h_max - h_ustalone) / h_ustalone
    #   plt.axhline(y=h_max, color='purple', linestyle='--', label=f'Maksymalne przeregulowanie κ = {kappa:.2%}')
    #   plt.axvline(x=T_max_przeregulowanie, color='purple', linestyle='--', label=f'Czas przeregulowania tp = {T_max_przeregulowanie:.4f}s')

    plt.title("Porównanie odpowiedzi skokowych układu inercyjnego 2. rzędu")
    plt.xlabel("Czas [s]"); plt.ylabel("Amplituda [V]")
    plt.legend(); plt.grid(True); plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()