import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Generuje i wyświetla charakterystyki Bodego dla układu inercyjnego 2. rzędu
    o transmitancji: G(s) = K * (omega_n^2) / (s^2 + 2*zeta*omega_n*s + omega_n^2)
    """
    
    # --- USTAWIENIA MODELI TEORETYCZNYCH ---
    # Wprowadź parametry dla każdego z trzech modeli.
    # K - wzmocnienie
    # omega_n - częstotliwość drgań nietłumionych [rad/s]
    # zeta - współczynnik tłumienia
    
    params_zestaw1 = {
        "K": 1, "omega_n": 5340, "zeta": 0.14,
        "kolor": "blue", "styl": "-", "nazwa": "Model 1 (np. symulacyjny)"
    }

    params_zestaw2 = {
        "K": 1, "omega_n": 5512, "zeta": 0.149,
        "kolor": "green", "styl": "--", "nazwa": "Model 2 (np. z ch. częstotliwościowej)"
    }

    params_zestaw3 = {
        "K": 1, "omega_n": 5406, "zeta": 0.129,
        "kolor": "red", "styl": ":", "nazwa": "Model 3 (np. z optymalizacji)"
    }
    
    zestawy_parametrow = [params_zestaw1, params_zestaw2, params_zestaw3]

    # --- USTAWIENIA POMIARÓW ---
    # Włącz/wyłącz rysowanie punktów pomiarowych
    RYSUJ_PUNKTY_POMIAROWE = True
    
    # Wklej tutaj swoje dane pomiarowe.
    punkty_pomiarowe = {
        "czestotliwosc_hz": [10,35,126,450,1000,1600],
        "amplituda_db": [0,-0.0852,0.0859,1.327,10.603,-12.517],# dla 10.601 Uout = 6.78V
        "faza_stopnie": [-3.6,-5.6,-5.6,-4.42,-85,-153],
        "kolor": "black", "marker": "x", "nazwa": "Dane pomiarowe"
    }

    # --- GENEROWANIE DANYCH ---
    # Definiujemy zakres pulsacji (rad/s) na skali logarytmicznej
    omega = np.logspace(2, 5, 1000) 

    fig, (ax_amplituda, ax_faza) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    print("📊 Obliczanie charakterystyk Bodego dla modeli 2. rzędu...")

    for params in zestawy_parametrow:
        K, omega_n, zeta = params["K"], params["omega_n"], params["zeta"]

        # 1. Charakterystyka amplitudowa
        # M(omega) = 20 * log10( K * omega_n^2 / sqrt((omega_n^2 - omega^2)^2 + (2*zeta*omega_n*omega)^2) )
        licznik_amp = K * omega_n**2
        mianownik_amp = np.sqrt((omega_n**2 - omega**2)**2 + (2 * zeta * omega_n * omega)**2)
        amplituda_db = 20 * np.log10(licznik_amp / mianownik_amp)
        
        # 2. Charakterystyka fazowa
        # phi(omega) = -arctan(2*zeta*omega_n*omega / (omega_n^2 - omega^2))
        faza_stopnie = -np.arctan2(2 * zeta * omega_n * omega, omega_n**2 - omega**2) * 180 / np.pi

        # Rysowanie
        label_modelu = f'{params["nazwa"]} (K={K}, ωn={omega_n}, ζ={zeta})'
        ax_amplituda.plot(omega, amplituda_db, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        ax_faza.plot(omega, faza_stopnie, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        
        print(f"  -> Narysowano model: K={K}, ωn={omega_n} rad/s, ζ={zeta}")

    # --- NANOSZENIE PUNKTÓW POMIAROWYCH ---
    if RYSUJ_PUNKTY_POMIAROWE:
        print("📍 Nanoszenie punktów pomiarowych na wykres...")
        omega_pomiarowe = 2 * np.pi * np.array(punkty_pomiarowe["czestotliwosc_hz"])
        ax_amplituda.scatter(omega_pomiarowe, punkty_pomiarowe["amplituda_db"],
                             color=punkty_pomiarowe["kolor"], marker=punkty_pomiarowe["marker"],
                             label=punkty_pomiarowe["nazwa"], zorder=5)
        ax_faza.scatter(omega_pomiarowe, punkty_pomiarowe["faza_stopnie"],
                        color=punkty_pomiarowe["kolor"], marker=punkty_pomiarowe["marker"],
                        label=punkty_pomiarowe["nazwa"], zorder=5)
        print("  -> Punkty pomiarowe zostały dodane.")

    # --- KONFIGURACJA WYKRESÓW ---
    ax_amplituda.set_ylabel("Amplituda [dB]")
    ax_amplituda.set_title("Charakterystyki Bodego dla układu inercyjnego 2. rzędu")
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