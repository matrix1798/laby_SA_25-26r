import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Generuje i wyświetla charakterystyki Bodego dla układu niemnimalnofazowego
    o transmitancji: G(s) = (1 - Tz*s) / (1 + Tp*s)
    """
    
    # --- USTAWIENIA MODELI TEORETYCZNYCH ---
    # Wprowadź parametry dla każdego z trzech modeli.
    # Tz - stała czasowa zera w prawej półpłaszczyźnie
    # Tp - stała czasowa bieguna w lewej półpłaszczyźnie
    
    params_zestaw1 = {
        "Tz": 0.1, "Tp": 1.0,
        "kolor": "blue", "styl": "-", "nazwa": "Model 1"
    }

    params_zestaw2 = {
        "Tz": 0.5, "Tp": 1.0,
        "kolor": "green", "styl": "--", "nazwa": "Model 2 (większy wpływ zera)"
    }

    params_zestaw3 = {
        "Tz": 0.1, "Tp": 0.2,
        "kolor": "red", "styl": ":", "nazwa": "Model 3 (szybszy biegun)"
    }
    
    zestawy_parametrow = [params_zestaw1, params_zestaw2, params_zestaw3]

    # --- GENEROWANIE DANYCH ---
    # Definiujemy zakres pulsacji (rad/s) na skali logarytmicznej
    omega = np.logspace(-2, 2, 1000) 

    fig, (ax_amplituda, ax_faza) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    print("📊 Obliczanie charakterystyk Bodego dla modeli niemnimalnofazowych...")

    for params in zestawy_parametrow:
        Tz, Tp = params["Tz"], params["Tp"]

        # 1. Charakterystyka amplitudowa
        # Amplituda jest taka sama jak dla układu minimalnofazowego (1+Tz*s)/(1+Tp*s)
        licznik_amp = np.sqrt(1 + (omega * Tz)**2)
        mianownik_amp = np.sqrt(1 + (omega * Tp)**2)
        amplituda_db = 20 * np.log10(licznik_amp / mianownik_amp)
        
        # 2. Charakterystyka fazowa
        # Faza = faza z zera (1-jωTz) - faza z bieguna (1+jωTp)
        # Faza = -arctan(ωTz) - arctan(ωTp)
        faza_zera_rad = -np.arctan(omega * Tz)
        faza_bieguna_rad = np.arctan(omega * Tp)
        faza_calkowita_stopnie = np.degrees(faza_zera_rad - faza_bieguna_rad)

        # Rysowanie
        label_modelu = f'{params["nazwa"]} (Tz={Tz}, Tp={Tp})'
        ax_amplituda.plot(omega, amplituda_db, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        ax_faza.plot(omega, faza_calkowita_stopnie, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        
        print(f"  -> Narysowano model: Tz={Tz}s, Tp={Tp}s")

    # --- KONFIGURACJA WYKRESÓW ---
    ax_amplituda.set_ylabel("Amplituda [dB]")
    ax_amplituda.set_title("Charakterystyki Bodego dla układu niemnimalnofazowego")
    ax_amplituda.grid(which='both', linestyle='--')
    ax_amplituda.legend()
    
    ax_faza.set_xlabel("Pulsacja ω [rad/s]")
    ax_faza.set_ylabel("Faza [°]")
    ax_faza.grid(which='both', linestyle='--')
    
    # Ustawienie "ticków" na osi Y wykresu fazowego
    ax_faza.set_yticks(np.arange(0, -181, -45))
    
    plt.xscale('log')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()