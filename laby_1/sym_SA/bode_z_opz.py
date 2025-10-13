import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Generuje i wyświetla charakterystyki Bodego dla obiektu 
    inercyjnego I rzędu Z OPÓŹNIENIEM TRANSPORTOWYM.
    """
    
    # --- USTAWIENIA MODELI TEORETYCZNYCH ---
    # Wprowadź parametry: K (wzmocnienie), T (stała czasowa) 
    # oraz Td (opóźnienie transportowe [s]) dla każdego modelu.
    
    params_zestaw1 = {
        "K": 0.843, "T": 0.00078, "Td": 0.00026,
        "kolor": "blue", "styl": "-", "nazwa": "Model z odp. skokowej"
    }

    params_zestaw2 = {
        "K": 0.843, "T": 0.00078, "Td": 0.00035,
        "kolor": "green", "styl": "--", "nazwa": "Model z ch. częstotliwościowej"
    }

    params_zestaw3 = {
        "K": 0.871, "T": 0.00077, "Td": 0.000412,
        "kolor": "orange"
        "", "styl": ":", "nazwa": "Model z optymalizacji"
    }
    
    zestawy_parametrow = [params_zestaw1, params_zestaw2, params_zestaw3]

    # --- USTAWIENIA POMIARÓW ---
    # Włącz/wyłącz rysowanie punktów pomiarowych.
    RYSUJ_PUNKTY_POMIAROWE = True
    
    # Wklej tutaj swoje dane pomiarowe.
    punkty_pomiarowe = {
        "czestotliwosc_hz": [10, 28 ,77, 159 ,600],
        "amplituda_db": [-1.18, -1.248, -1.735, -3.188,-10.91],
        "faza_stopnie": [-4.25, -11.87 ,-31.85 ,-61.15 ,-159.98], # Przykładowe dane z opóźnieniem
        "kolor": "black", "marker": "x", "nazwa": "Dane pomiarowe"
    }

    # --- GENEROWANIE DANYCH ---
    # Definiujemy zakres pulsacji (rad/s) na skali logarytmicznej
    omega = np.logspace(-1, 4, 1000) 

    fig, (ax_amplituda, ax_faza) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    print("📊 Obliczanie charakterystyk Bodego dla modeli z opóźnieniem...")

    for params in zestawy_parametrow:
        K, T, Td = params["K"], params["T"], params["Td"]

        # 1. Charakterystyka amplitudowa (pozostaje BEZ ZMIAN)
        # Opóźnienie transportowe nie wpływa na amplitudę.
        amplituda_db = 20 * np.log10(K / np.sqrt(1 + (omega * T)**2))
        
        # 2. Charakterystyka fazowa (NOWA FORMUŁA)
        # Faza = faza z członu inercyjnego + faza z opóźnienia
        # phi(omega) = -arctan(omega*T) - omega*Td
        faza_inercyjna_rad = -np.arctan(omega * T)
        faza_opoznienia_rad = -omega * Td
        faza_calkowita_stopnie = np.degrees(faza_inercyjna_rad + faza_opoznienia_rad)

        # Rysowanie
        label_modelu = f'{params["nazwa"]} '
        ax_amplituda.plot(omega, amplituda_db, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        ax_faza.plot(omega, faza_calkowita_stopnie, label=label_modelu, color=params["kolor"], linestyle=params["styl"])
        
        print(f"  -> Narysowano model: K={K}, T={T}s, Td={Td}s")

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
    ax_amplituda.set_title("Charakterystyki Bodego dla obiektu I rzędu z opóźnieniem")
    ax_amplituda.grid(which='both', linestyle='--')
    ax_amplituda.legend()
    
    ax_faza.set_xlabel("Pulsacja ω [rad/s]")
    ax_faza.set_ylabel("Faza [°]")
    ax_faza.grid(which='both', linestyle='--')
    
    # Dodanie ważnych wartości (-90, -180) do osi Y wykresu fazowego
    ticks = list(ax_faza.get_yticks())
    for val in [-90, -180]:
        if val not in ticks:
            ticks.append(val)
    ticks.sort()
    ax_faza.set_yticks(ticks)
    
    plt.xscale('log')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()