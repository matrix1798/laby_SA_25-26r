import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Główna funkcja programu, która generuje i wyświetla
    charakterystyki Bodego dla obiektów inercyjnych I rzędu
    oraz nanosi na nie punkty pomiarowe.
    """
    
    # --- USTAWIENIA MODELI TEORETYCZNYCH ---
    # Wprowadź tutaj parametry K i T dla modeli teoretycznych.
    
    params_pomiarowe = {
        "K": 0.871,
        "T": 0.00078,
        "kolor": "blue",
        "styl": "-",  # Linia ciągła
        "nazwa": "Model z odp. skokowej"
    }

    # Zestaw 2: Parametry z charakterystyki częstotliwościowej
    params_czestotliwosciowe = {
        "K": 0.842,
        "T": 0.000758,
        "kolor": "green",
        "styl": "--", # Linia przerywana
        "nazwa": "Model z ch. częstotliwościowej"
    }

    # Zestaw 3: Parametry z optymalizacji
    params_optymalizacja = {
        "K": 0.8366,
        "T": 0.00074,
        "kolor": "orange",
        "styl": ":",  # Linia kropkowana
        "nazwa": "Model z optymalizacji"
    }
    
    zestawy_parametrow = [params_pomiarowe, params_czestotliwosciowe, params_optymalizacja]

    # --- USTAWIENIA POMIARÓW (NOWA SEKCJA) ---
    # Włącz/wyłącz rysowanie punktów pomiarowych
    RYSUJ_PUNKTY_POMIAROWE = True
    
    # Wklej tutaj swoje dane pomiarowe.
    # Upewnij się, że każda lista ma tyle samo elementów!
    punkty_pomiarowe = {
        "czestotliwosc_hz": [10,25,63,159,400],
        "amplituda_db": [-1.49,-1.52,-1.98,-2.67,-12.06],
        "faza_stopnie": [-3.5, -7.98, -15.9, -29.9, -71.54],
        "kolor": "black",
        "marker": "x", # Możesz użyć 'o', 's', '^' itp.
        "nazwa": "Dane pomiarowe"
    }

    # --- GENEROWANIE DANYCH DLA MODELI ---
    omega = np.logspace(0, 5, 500) # Zakres pulsacji [rad/s]
    fig, (ax_amplituda, ax_faza) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    print("📊 Obliczanie charakterystyk Bodego dla poszczególnych modeli...")

    for params in zestawy_parametrow:
        K, T = params["K"], params["T"]
        amplituda_db = 20 * np.log10(K / np.sqrt(1 + (omega * T)**2))
        faza_stopnie = -np.arctan(omega * T) * 180 / np.pi
        
        ax_amplituda.plot(omega, amplituda_db, label=f'{params["nazwa"]}', color=params["kolor"], linestyle=params["styl"])
        ax_faza.plot(omega, faza_stopnie, label=f'{params["nazwa"]} ', color=params["kolor"], linestyle=params["styl"])
        print(f"  -> Narysowano model: K={K}, T={T}s")

    # --- NANOSZENIE PUNKTÓW POMIAROWYCH NA WYKRES ---
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
    ax_amplituda.set_title("Porównanie modeli z danymi pomiarowymi na wykresach Bodego")
    ax_amplituda.grid(which='both', linestyle='--')
    ax_amplituda.legend()
    
    ax_faza.set_xlabel("Pulsacja ω [rad/s]")
    ax_faza.set_ylabel("Faza [°]")
    ax_faza.grid(which='both', linestyle='--')
    # Legenda dla wykresu fazowego jest opcjonalna, bo etykiety są te same
    # ax_faza.legend() 

    ticks = list(ax_faza.get_yticks())
    if -90 not in ticks:
        ticks.append(-90)
    ticks.sort()
    ax_faza.set_yticks(ticks)
        
    plt.xscale('log')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()