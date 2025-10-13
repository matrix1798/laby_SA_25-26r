import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Generuje i wyświetla charakterystyki Bodego dla układu niemnimalnofazowego
    o transmitancji: G(s) = (1 - Tz*s) / (1 + Tp*s)
    """
    
    # --- USTAWIENIA MODELI TEORETYCZNYCH ---
    # Wprowadź parametry dla każdego z trzech modeli.
    params_zestaw1 = {
        "Tz":  0.00038, "Tp":0.00012,
        "kolor": "blue", "styl": "-", "nazwa": "Model z odp. skokowej"
    }
    params_zestaw2 = {
        "Tz": 0.0003112, "Tp":  0.00009,
        "kolor": "green", "styl": "--", "nazwa": "Model z char. częstotliwościowej"
    }
    params_zestaw3 = {
        "Tz": 0.0003334, "Tp": 0.0001,
        "kolor": "orange", "styl": ":", "nazwa": "Model z optymalizacji"
    }
    zestawy_parametrow = [params_zestaw1, params_zestaw2, params_zestaw3]

    # --- USTAWIENIA POMIARÓW (NOWA SEKCJA) ---
    # Włącz/wyłącz rysowanie punktów pomiarowych
    RYSUJ_PUNKTY_POMIAROWE = True
    
    # Wklej tutaj swoje dane pomiarowe.
    punkty_pomiarowe = {
        "czestotliwosc_hz": [10, 100, 1000],
        "amplituda_db": [0.214,-0.0868,6.59],
        "faza_stopnie": [-6, -17.5, -87],
        "kolor": "black",
        "marker": "x",
        "nazwa": "Dane pomiarowe"
    }

    # --- GENEROWANIE DANYCH ---
    omega = np.logspace(1, 6, 1000) 
    fig, (ax_amplituda, ax_faza) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    print("📊 Obliczanie charakterystyk Bodego dla modeli niemnimalnofazowych...")

    for params in zestawy_parametrow:
        Tz, Tp = params["Tz"], params["Tp"]

        licznik_amp = np.sqrt(1 + (omega * Tz)**2)
        mianownik_amp = np.sqrt(1 + (omega * Tp)**2)
        amplituda_db = 20 * np.log10(licznik_amp / mianownik_amp)
        
        faza_zera_rad = -np.arctan(omega * Tz)
        faza_bieguna_rad = np.arctan(omega * Tp)
        faza_calkowita_stopnie = np.degrees(faza_zera_rad - faza_bieguna_rad)

        label_modelu = f'{params["nazwa"]}'
        ax_amplituda.plot(omega, amplituda_db, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        ax_faza.plot(omega, faza_calkowita_stopnie, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        print(f"  -> Narysowano model: Tz={Tz}s, Tp={Tp}s")

    # --- NANOSZENIE PUNKTÓW POMIAROWYCH NA WYKRES (NOWA SEKCJA) ---
    if RYSUJ_PUNKTY_POMIAROWE:
        print("📍 Nanoszenie punktów pomiarowych na wykres...")
        # Konwersja częstotliwości z Hz na pulsację w rad/s (ω = 2 * pi * f)
        omega_pomiarowe = 2 * np.pi * np.array(punkty_pomiarowe["czestotliwosc_hz"])
        
        # Rysowanie punktów na wykresie amplitudowym
        ax_amplituda.scatter(omega_pomiarowe, punkty_pomiarowe["amplituda_db"],
                             color=punkty_pomiarowe["kolor"],
                             marker=punkty_pomiarowe["marker"],
                             label=punkty_pomiarowe["nazwa"],
                             zorder=5) # zorder=5 sprawia, że punkty są na wierzchu

        # Rysowanie punktów na wykresie fazowym
        ax_faza.scatter(omega_pomiarowe, punkty_pomiarowe["faza_stopnie"],
                        color=punkty_pomiarowe["kolor"],
                        marker=punkty_pomiarowe["marker"],
                        label=punkty_pomiarowe["nazwa"],
                        zorder=5)
        print("  -> Punkty pomiarowe zostały dodane.")


    # --- KONFIGURACJA WYKRESÓW ---
    ax_amplituda.set_ylabel("Amplituda [dB]")
    ax_amplituda.set_title("Charakterystyki Bodego dla układu niemnimalnofazowego")
    ax_amplituda.grid(which='both', linestyle='--')
    ax_amplituda.legend()
    
    ax_faza.set_xlabel("Pulsacja ω [rad/s]")
    ax_faza.set_ylabel("Faza [°]")
    ax_faza.grid(which='both', linestyle='--')
    
    ax_faza.set_yticks(np.arange(0, -181, -45))
    
    plt.xscale('log')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()