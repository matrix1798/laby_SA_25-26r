import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d  # filtr Gaussa

def load_data(file_path):
    column_name = pd.read_csv(file_path, nrows=0).columns
    metadane_df = pd.read_csv(file_path, skiprows=1, nrows=1, header=None, names=column_name)

    start = metadane_df['Start'].iloc[0]
    increment = metadane_df['Increment'].iloc[0]

    dane_df = pd.read_csv(file_path, skiprows=2, header=None, usecols=[0, 1])
    dane_df.columns = ['CH1_offset', 'CH2_offset']

    number_of_points = len(dane_df)
    time = np.arange(start, start + number_of_points * increment, increment)

    dane_df['Time'] = time
    dane_df = dane_df[['Time', 'CH1_offset', 'CH2_offset']]
    return dane_df

# --- wczytanie danych ---
file_path = ['../data/NewFile7.csv', '../data/NewFile9.csv', '../data/NewFile10.csv']

try:
    df1 = load_data(file_path[0])
    df2 = load_data(file_path[1])
    df3 = load_data(file_path[2])
except FileNotFoundError as e:
    print(f"Sprawdź, czy plik '{e.filename}' istnieje i czy ścieżka jest poprawna.")
    exit()
except Exception as e:
    print(f"Wystąpił nieoczekiwany błąd: {e}")
    exit()

print('Wczytanie plików zakończone')

# --- redukcja próbkowania ---
n = 1
df1 = df1.iloc[::n, :].reset_index(drop=True)
df2 = df2.iloc[::n, :].reset_index(drop=True)
df3 = df3.iloc[::n, :].reset_index(drop=True)

# --- pobranie sygnału prostokątnego z pliku (zastąpione idealnym prostokątem) ---
signal = pd.DataFrame({
    'Time': df1['Time'],
    'input': np.where(df1['Time'] >= 0, 2.0, 0.0)
})

# --- mocne wygładzenie sygnałów wyjściowych ---
sigma_out = 8  # większe = bardziej gładko
df1['CH2'] = gaussian_filter1d(df1['CH2_offset'] - df1['CH2_offset'].iloc[0], sigma_out)
df2['CH2'] = gaussian_filter1d(df2['CH2_offset'] - df2['CH2_offset'].iloc[0], sigma_out)
df3['CH2'] = gaussian_filter1d(df3['CH2_offset'] - df3['CH2_offset'].iloc[0], sigma_out)

plt.figure(figsize=(12, 7))
chart_width = 0.5
range_time = int(df1['Time'].count() * chart_width)
start_time = df1['Time'].iloc[0]
end_time = df1['Time'].iloc[range_time - 1]

# maska: tylko dane od t >= 0
mask1 = df1['Time'] >= 0
mask2 = df2['Time'] >= 0
mask3 = df3['Time'] >= 0

plt.plot(signal['Time'], signal['input'], label='Sygnał wejściowy', color='blue')

# linia po całości
plt.plot(df1['Time'], df1['CH2'], color='green', label='Sygnał wyjściowy k1 = 1,02')
plt.plot(df2['Time'], df2['CH2'], color='orange', label='Sygnał wyjściowy k2 = 2,47')
plt.plot(df3['Time'], df3['CH2'], color='red', label='Sygnał wyjściowy k3 = 4,47')

# markery tylko dla t >= 0
plt.plot(df1['Time'][mask1], df1['CH2'][mask1],
         linestyle='', marker='o', markevery=50, markersize=10,
         markerfacecolor='none', color='green')
plt.plot(df2['Time'][mask2], df2['CH2'][mask2],
         linestyle='', marker='o', markevery=50, markersize=10,
         markerfacecolor='none', color='orange')
plt.plot(df3['Time'][mask3], df3['CH2'][mask3],
         linestyle='', marker='o', markevery=50, markersize=10,
         markerfacecolor='none', color='red')

plt.xlim(start_time, end_time)
plt.xlabel('Czas (s)')
plt.ylabel('Napięcie (V)')
plt.legend()
plt.grid(True)
plt.show()

