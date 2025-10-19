import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def load_data(file_path):

    column_name = pd.read_csv(file_path,nrows=0).columns
    metadane_df = pd.read_csv(file_path,skiprows=1,nrows=1,header=None,names=column_name)

    start = metadane_df['Start'].iloc[0]
    increment = metadane_df['Increment'].iloc[0]

    dane_df = pd.read_csv(file_path, skiprows=2, header =None, usecols = [0,1])
    dane_df.columns = ['CH1_offset','CH2_offset']

    number_of_points = len(dane_df)
    time = np.arange(start,start+ number_of_points * increment , increment)

    dane_df['Time'] = time
    dane_df = dane_df[['Time','CH1_offset','CH2_offset']]

    return dane_df

# load data
file_path = ['../data/NewFile2.csv','../data/NewFile3.csv','../data/NewFile4.csv']

try:
    df1 = load_data(file_path[0])
    df2 = load_data(file_path[1])
    df3 = load_data(file_path[2])

except FileNotFoundError as e:
    # 'e' to obiekt błędu, który zawiera informacje, np. nazwę brakującego pliku
    print(f"Sprawdź, czy plik '{e.filename}' na pewno istnieje i czy ścieżka jest poprawna.")
    # Zakończ program, ponieważ bez danych dalsze operacje nie mają sensu
    exit()
except Exception as e:
    # Dobra praktyka to łapanie też innych, nieprzewidzianych błędów
    print(f"Wystąpił nieoczekiwany błąd: {e}")
    exit()

print('Wczytanie plikow')

# Shift the offset

offset = df1['CH1_offset'].iloc[0]

signal = df1[['Time','CH1_offset']].copy()
signal['CH1_offset'] = signal['CH1_offset'] - offset
signal.rename(columns={'CH1_offset':'input'},inplace=True)

df1['CH2'] = df1['CH2_offset'] - offset
df2['CH2'] = df2['CH2_offset'] - offset
df3['CH2'] = df3['CH2_offset'] - offset

# visualization data

plt.figure(figsize=(12,7))

# range
chart_width = 0.5
range_time = int(df1['Time'].count() * chart_width)
start_time = df1['Time'].iloc[0]
end_time = df1['Time'].iloc[range_time-1]

plt.plot(signal['Time'],signal['input'],label='Sygnał wejściowy', color='blue',marker='o',markevery = 50,markersize=10,markerfacecolor='none')
plt.plot(df1['Time'],df1['CH2'],label='Sygnał wyjściowy k1 = ...',color = 'green',marker='o',markevery = 50,markersize=10,markerfacecolor='none')
plt.plot(df2['Time'],df2['CH2'],label='Sygnał wyjściowy k2 = ...',color = 'orange',marker='o',markevery = 50,markersize=10,markerfacecolor='none')
plt.plot(df3['Time'],df3['CH2'],label='Sygnał wyjściowy k3 = ...',color='red',marker='o',markevery = 50,markersize=10,markerfacecolor='none')

plt.xlim(start_time,end_time)

plt.xlabel('Czas (s)')
plt.ylabel('Napięcie (V)')
plt.legend()
plt.grid(True)

plt.show()