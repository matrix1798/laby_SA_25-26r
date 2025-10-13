import pandas as pd
import numpy as np
import os

def znajdz_rampe(df, kolumna="CH2", prog_min=0.15, prog_max=0.75):
    """
    Znajduje początek i koniec rampy sygnału wyjściowego.
    Wersja odporna – zawsze coś zwraca.
    """
    y = df[kolumna].values
    if len(y) < 10:
        return 0, len(y) - 1

    y_min, y_max = np.min(y), np.max(y)
    zakres = y_max - y_min

    # Jeżeli sygnał jest płaski – zwracamy cały zakres
    if zakres < 1e-6:
        return 0, len(y) - 1

    start_val = y_min + prog_min * zakres
    end_val = y_min + prog_max * zakres

    # Indeksy, gdzie przekracza progi
    above_start = np.where(y > start_val)[0]
    above_end = np.where(y > end_val)[0]

    if len(above_start) == 0 or len(above_end) == 0:
        # W razie problemu — bierz środek sygnału
        start_index = int(len(y) * 0.25)
        end_index = int(len(y) * 0.75)
    else:
        start_index = above_start[0]
        end_index = above_end[0]

    # Upewnij się, że Δx > 0
    if end_index <= start_index:
        end_index = start_index + max(10, len(y) // 20)

    # Ogranicz maksymalny czas rampy
    end_index = min(end_index, len(y) - 1)

    return start_index, end_index

def process_file(filename, uin):
    """
    Wczytuje plik CSV, automatycznie wykrywa rampę i oblicza Ti = (Uin * Δx) / Δy.
    """
    try:
        with open(filename, 'r') as f:
            f.readline()
            meta_line = f.readline().strip().split(',')
            increment = float(meta_line[3])
    except (IOError, ValueError) as e:
        print(f"Błąd podczas odczytu pliku {filename}: {e}")
        return None

    df = pd.read_csv(filename, skiprows=2, header=None, usecols=[0, 1], names=['CH1', 'CH2'], dtype=np.float64)

    # Usuwamy ewentualne NaN
    df = df.dropna()

    # Wykrycie rampy
    idx_start, idx_end = znajdz_rampe(df, "CH2", prog_min=0.10, prog_max=0.90)

    # Dane rampy
    delta_y = df["CH2"].iloc[idx_end] - df["CH2"].iloc[idx_start]
    delta_x = (idx_end - idx_start) * increment

    # Obliczenie Ti ze wzoru Ti = (Uin * Δx) / Δy
    Ti_s = (uin * delta_x) / delta_y
    Ti_ms = Ti_s * 1000

    print(f"--- Wyniki dla pliku {filename} ---")
    print(f"Zadane napięcie wejściowe Uin: {uin:.4f} V")
    print(f"Δx = {delta_x:.6f} s, Δy = {delta_y:.4f} V")
    print(f"Wyznaczony czas całkowania Ti = {Ti_ms:.4f} ms")
    print("-----------------------------------")

    return Ti_s, delta_x, delta_y


def main():
    """
    Główna funkcja analizy dla wszystkich plików.
    """
    files_and_uins = {
        'fol_gabriel/NewFile4.csv': 2.0,
        'fol_gabriel/NewFile5.csv': 4.0,
        'fol_gabriel/NewFile6.csv': 8.0
    }

    wyniki = []

    for file, uin_value in files_and_uins.items():
        if os.path.exists(file):
            Ti_s, dx, dy = process_file(file, uin_value)
            wyniki.append([file, uin_value, dx, dy, Ti_s * 1000])
        else:
            print(f"❌ Nie znaleziono pliku: {file}")

    if wyniki:
        df = pd.DataFrame(wyniki, columns=["Plik", "Uin [V]", "Δx [s]", "Δy [V]", "Ti [ms]"])
        print("\n===== Tabela wyników =====")
        print(df.to_string(index=False))

        Ti_mean = np.mean(df["Ti [ms]"])
        Ti_std = np.std(df["Ti [ms]"])
        print("\nŚrednia Ti = {:.4f} ms, odchylenie standardowe = {:.4f} ms".format(Ti_mean, Ti_std))
    else:
        print("\nBrak danych do analizy.")


if __name__ == "__main__":
    main()
