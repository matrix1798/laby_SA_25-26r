import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    filename = "fol_gabriel/NewFile1.csv"
    filepath = os.path.join(os.path.dirname(__file__), filename)

    if not os.path.exists(filepath):
        print(f"❌ Plik '{filename}' nie został znaleziony w folderze:\n{os.path.dirname(__file__)}")
        return

    # Sprawdzenie formatu pliku
    with open(filepath, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()

    # --- FORMAT 2: z nagłówkami CH1, CH2, Start, Increment
    if first_line.startswith("CH1"):
        print("📗 Wykryto format typu 2 (CH1, CH2, Start, Increment)")
        df = pd.read_csv(filepath, skiprows=2, header=None)

        # Odczytaj Start i Increment z drugiej linii
        parts = second_line.split(",")
        start = float(parts[2])
        increment = float(parts[3])

        # Dane: CH1, CH2 (pierwsze 2 kolumny)
        df = df.iloc[:, :2]
        df.columns = ["in", "out"]

        # Dodaj kolumnę czasu
        df["czas"] = [start + i * increment for i in range(len(df))]
        df = df[["czas", "in", "out"]]
        # Usuń wiersze, w których brakuje danych
        df = df.dropna()

        # Dla pewności — upewnij się, że kolumny są numeryczne
        df = df.astype(float)

    # --- FORMAT 1: klasyczny z trzema kolumnami liczbowymi
    else:
        print("📘 Wykryto format typu 1 (czas, sygnał_1, sygnał_2)")
        df = pd.read_csv(filepath, skiprows=2, header=None)
        df.columns = ["czas", "in", "out"]
        # Usuń wiersze, w których brakuje danych
        df = df.dropna()

        # Dla pewności — upewnij się, że kolumny są numeryczne
        df = df.astype(float)

    # --- PRZESUNIĘCIE DO ZERA ---
    in_offset = -df["in"].iloc[0]   # różnica pierwszej próbki od zera
    out_offset = -df["out"].iloc[0] # analogicznie dla OUT

    df["in"] = df["in"] + in_offset
    df["out"] = df["out"] + out_offset

    #print(f"✅ Sygnały przesunięte: IN +({in_offset:+.3e}), OUT +({out_offset:+.3e})")

    #Parametry odpowiedzi skokowej
    kp = df["out"][len(df["out"])-1]/df["in"][len(df["in"])-1]
    print(f"K_p = {kp:4f}")
    
    df['roznica'] = abs(df['out'] - kp*0.6321)
    idx = df['roznica'].idxmin()
    Tp = df.loc[idx,'czas']
    print(f'T_p = {Tp}')

    # --- RYSOWANIE ---
    plt.figure(figsize=(10, 6))
    plt.plot(df["czas"], df["in"], label="Wejście (IN)")
    plt.plot(df["czas"], df["out"], label="Wyjście (OUT)")

    plt.title("Wizualizacja sygnałów IN/OUT (przesunięte tak, by startowały od 0 V)")
    plt.xlabel("Czas [s]")
    plt.ylabel("Amplituda [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
