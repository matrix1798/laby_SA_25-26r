import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Generuje wykresy Bodego (amplituda i faza) na podstawie punktów pomiarowych
    dla trzech układów: k1, k2, k3.
    """

    # --- DANE POMIAROWE ---

    # (k3)
    f_k3 = np.array([10, 100, 200, 500, 795, 1000])
    Uin_k3 = np.array([2, 2, 2, 2, 2, 2])
    Uout_k3 = np.array([2.06, 2.06, 2.12, 3.58, 1.42, 0.58])
    phi_k3 = np.array([1.8, 1.8, 30.24, 83.8, -157.4, -133.2])

    # (k2)
    f_k2 = np.array([10, 100, 200, 500, 680, 1000])
    Uin_k2 = np.array([2, 2, 2, 2, 2, 2])
    Uout_k2 = np.array([2.06, 2.04, 1.96, 1.19, 1.40, 0.38])
    phi_k2 = np.array([2.16, 22.32, 39.6, 112.5, -177.5, -137.5])

    # (k1)
    f_k1 = np.array([10, 100, 200, 250, 500, 800])
    Uin_k1 = np.array([2, 2, 2, 2, 2, 2])
    Uout_k1 = np.array([2.06, 1.86, 1.58, 1.46, 1.1, 0.46])
    phi_k1 = np.array([3.9, 32.8, 59.7, 72.9, 132.3, -162.4])

    # --- OBLICZENIA ---
    # Konwersja częstotliwości na pulsację ω = 2πf
    omega_k1 = 2 * np.pi * f_k1
    omega_k2 = 2 * np.pi * f_k2
    omega_k3 = 2 * np.pi * f_k3

    # Wzmocnienie w dB
    A_k1 = 20 * np.log10(Uout_k1 / Uin_k1)
    A_k2 = 20 * np.log10(Uout_k2 / Uin_k2)
    A_k3 = 20 * np.log10(Uout_k3 / Uin_k3)

    # --- WYKRESY BODEGO ---
    fig, (ax_amp, ax_phase) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # --- Amplituda ---
    ax_amp.semilogx(omega_k1, A_k1, 'o-', color='green', label='k1 = 0.52')
    ax_amp.semilogx(omega_k2, A_k2, 'o-', color='orange', label='k2 = 1.12')
    ax_amp.semilogx(omega_k3, A_k3, 'o-', color='red', label='k3 = 1.67')
    ax_amp.set_ylabel('Amplituda [dB]')
    ax_amp.set_title('Charakterystyki Bodego (punkty pomiarowe)')
    ax_amp.grid(True, which='both', linestyle='--')
    ax_amp.legend()

    # --- Faza ---
    ax_phase.semilogx(omega_k1, phi_k1, 'o-', color='green', label='k1 = 0.52')
    ax_phase.semilogx(omega_k2, phi_k2, 'o-', color='orange', label='k2 = 1.12')
    ax_phase.semilogx(omega_k3, phi_k3, 'o-', color='red', label='k3 = 1.67')
    ax_phase.set_xlabel('Pulsacja ω [rad/s]')
    ax_phase.set_ylabel('Faza [°]')
    ax_phase.grid(True, which='both', linestyle='--')
    ax_phase.legend()

    plt.tight_layout()
    plt.savefig('bode_punkty.png', dpi=300)
    plt.show()

    print("✅ Zapisano wykres Bodego jako 'bode_punkty.png'")

if __name__ == "__main__":
    main()
