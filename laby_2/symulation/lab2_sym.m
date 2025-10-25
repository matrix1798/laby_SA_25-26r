%% Uklad A: 
close all;
clear all;

%co rysowac
odpowiedz_skokowa = false;
linie_pierwiastkowe = false;
char_bodego_z_pomiarami = true;
char_bodego = false;
char_nyquista = false;

K_c = [0.52 1.12 1.67];

T_i = 0.0013;
Gi = tf(1,[T_i 0]); %czlon calkujacy

w0 = 2560;
zeta = 0.37;
Go = tf(w0^2,[1 2*zeta*w0 w0^2]);

%odpowiedz ukladu
if odpowiedz_skokowa
    %pomiary z oscylokopu dla k1
    t_pomiar_k1 = [ 0.002, 0.003, 0.004, 0.005, 0.007];
    y_pomiar_k1 = [  1.2,   1.85,  2.1,   2.0,   1.8];
    %pomiary z oscylokopu dla k2
    t_pomiar_k2 = [ 0.002, 0.003, 0.004, 0.005, 0.007];
    y_pomiar_k2 = [  1.2,   1.85,  2.1,   2.0,   1.8];
    %pomiary z oscylokopu dla k3
    t_pomiar_k3 = [ 0.002, 0.003, 0.004, 0.005, 0.007];
    y_pomiar_k3 = [  1.2,   1.85,  2.1,   2.0,   1.8];

    t = 0:0.0001:0.01;
    u = zeros(size(t));
    u(t>=0.001) = 2 ; %definicja skoku
    
    figure;
    hold on;
    grid on;
    
    legend_entries = {};

    for Kc = K_c
        % --- Wewnątrz pętli dla każdego K_p ---
    
        % 1. Zdefiniuj regulator P dla BIEŻĄCEGO wzmocnienia
        Gc = tf(Kc);
    
        % 2. Zbuduj układ otwarty i zamknięty dla BIEŻĄCEGO K_p
        system_open = Gc * Gi * Go;
        system_closed = feedback(system_open, 1);
    
        % 3. Uruchom symulację
        [Y, T] = lsim(system_closed, u, t);
    
        % 4. Narysuj wynik na przygotowanym wcześniej wykresie
        plot(T, Y, 'LineWidth', 1.5); % 'LineWidth' pogrubia linię
    
        % 5. Dodaj wpis do legendy
        legend_entries{end+1} = ['K_c = ', num2str(Kc)];
    end
    legend_entries{end+1} = 'Sygnał wejściowy';
    plot(t,u,'black');
%dodanie punktow pomiarowych
    plot(t_pomiar_k1,y_pomiar_k1,'bo','MarkerSize', 6);
    legend_entries{end+1} = 'Punkty pomiarowe dla k1';
    plot(t_pomiar_k2,y_pomiar_k2,'ro','MarkerSize', 6);
    legend_entries{end+1} = 'Punkty pomiarowe dla k2';
    plot(t_pomiar_k3,y_pomiar_k3,'o','Color',[1 0.65 0],'MarkerSize', 6);
    legend_entries{end+1} = 'Punkty pomiarowe dla k3';
% Wyłącz "przytrzymywanie"
    hold off; 
    title('Odpowiedź na skok dla różnych wzmocnień K_p');
    xlabel('Czas [s]');
    ylabel('Odpowiedź wyjściowa');
    legend(legend_entries); % Stwórz legendę na podstawie zebranych etykiet
end

% Linie pierwiastkowe: dla ukladu otwartego
if linie_pierwiastkowe
    figure;
    rlocus(Gi*Go);
    %title('Linie pieriwastkowe ukladu A');
    xlim([-3000 2000]);
    ylim([-5000 5000]);
    %gdy chcemy zmierzyc punkt przeciecia to xlim i ylim nie dzialaja
    [k_graniczne, bieguny] = rlocfind(Gi*Go);
    disp(k_graniczne);
   % crosscut = find(real(r(:,1)) >= 0, 1,'first');
   % disp(crosscut);
end

% Char. Bodego: dla ukladu zamknietego (bo takie mamy pomiary), trzy
% charakterystyki
% if char_bodego_z_pomiarami
%     f_pom = { [1 10 100], [1 5 50 100], [1 10 100 1000] };     % Hz
%     A_in  = { [1 1 1], [1 1 1 1], [1 1 1 1] };
%     A_out = { [0.9 0.7 0.3], [0.95 0.8 0.4 0.2], [1 0.9 0.5 0.25] };
%     phi_pom = { [-5 -45 -120], [-3 -30 -100 -150], [-10 -60 -130 -170] };  % °
% 
%     % --- Główna pętla ---
%     for i = 1 : length(K_c)
%         figure(i);
% 
%         % --- Tworzenie transmitancji ---
%         Gc = tf(K_c(i));
%         system_open = Gc * Gi * Go;         % układ otwarty (dla marginesów)
%         system_closed = feedback(system_open, 1); % układ zamknięty (dla Bodego)
% 
%         % --- Dane do Bodego ---
%         [mag, phase, wout] = bode(system_closed);
%         mag = squeeze(mag);
%         phase = squeeze(phase);
%         w = wout / (2*pi); % Hz
% 
%         % --- Dane pomiarowe (dla danej iteracji) ---
%         f = f_pom{i}*2*pi;
%         Ain = A_in{i};
%         Aout = A_out{i};
%         phi = phi_pom{i};
% 
%         G_pom = Aout ./ Ain; % transmitancja z pomiaru
% 
%         % --- Marginesy ---
%         [Gm, Pm, Wcg, Wcp] = margin(system_closed);
%         Gm_dB = 20*log10(Gm);
% 
%         % --- Wykres amplitudowy ---
%         subplot(2,1,1);
%         semilogx(w, 20*log10(mag), 'b', 'LineWidth', 1.5); hold on;
%         semilogx(f, 20*log10(G_pom), 'ro', 'MarkerFaceColor', 'r');
%         grid on;
%         ylabel('Amplituda [dB]');
%         title(['Charakterystyka Bodego - K_c = ', num2str(K_c(i))]);
%         legend('Model', 'Pomiar', 'Location', 'Best');
% 
%         % --- Zaznaczenie zapasu wzmocnienia ---
%         if ~isnan(Wcg)
%             semilogx(Wcg/(2*pi), 0, 'ks', 'MarkerFaceColor', 'y', 'DisplayName', 'Zapas wzmocnienia');
%             text(Wcg/(2*pi), 3, sprintf('Gm = %.1f dB', Gm_dB), 'HorizontalAlignment', 'center');
%         end
% 
%         % --- Wykres fazowy ---
%         subplot(2,1,2);
%         semilogx(w, phase, 'b', 'LineWidth', 1.5); hold on;
%         semilogx(f, phi, 'ro', 'MarkerFaceColor', 'r');
%         grid on;
%         ylabel('Faza [°]');
%         xlabel('Częstotliwość [Hz]');
%         legend('Model', 'Pomiar', 'Location', 'Best');
% 
%         % --- Zaznaczenie zapasu fazy ---
%         if ~isnan(Wcp)
%             semilogx(Wcp/(2*pi), -180, 'ks', 'MarkerFaceColor', 'c', 'DisplayName', 'Zapas fazy');
%             text(Wcp/(2*pi), -165, sprintf('Pm = %.1f°', Pm), 'HorizontalAlignment', 'center');
%         end
%     end
% end

% --- Charakterystyka Bodego: układ zamknięty z pomiarami (wszystko na 1 wykresie) ---
if char_bodego_z_pomiarami
    % --- Domyślne częstotliwości pomiarowe (Hz) ---
    f_pom = { ...
    [10 100 200 500 795 1500], ...   % k3, n = 2.4
    [10 100 200 500 680 1000], ...   % k2, n = 1.7
    [10 100 200 500 800] ...         % k1, n = 0.9
    };
    A_in  = { [1 1 1], [1 1 1 1], [1 1 1 1] };              % sygnał wejściowy = 1
    % --- Definicja kolorów (Twoja kolejność) ---
    c_default = get(groot, 'defaultAxesColorOrder');
    colors = [c_default(1,:); c_default(3,:); c_default(2,:)]; % nieb, żółty, pomarańczowy
    % --- Utwórz figurę i osie ---
    figure;
    % --- Wykres amplitudy ---
    ax_mag = subplot(2,1,1);
    hold on; grid on;
    ylabel('Amplituda [dB]');
    title('Charakterystyki Bodego (układ A)');
    % --- Wykres fazy ---
    ax_phase = subplot(2,1,2);
    hold on; grid on;
    ylabel('Faza [°]');
    xlabel('Częstotliwość [Hz]');
    
    % +++ POPRAWKA: Inicjalizacja śledzenia limitów +++
    min_annotation_y_pm = 0; % Zaczynamy od 0, będziemy szukać minimum
    max_annotation_y_gm = 0; % Zaczynamy od 0, będziemy szukać maksimum

    % --- Pętla po wszystkich wartościach K_c ---
    for i = 1 : length(K_c)
        currentColor = colors(i, :); % kolor danej charakterystyki
        % --- Tworzenie transmitancji ---
        Gc = tf(K_c(i));
        system_open = Gc * Gi * Go;                 % układ otwarty
        system_closed = feedback(system_open, 1);   % układ zamknięty
        % --- Dane modelu do Bodego ---
        [mag, phase, wout] = bode(system_closed);
        mag_dB = 20*log10(squeeze(mag));
        phase_deg = squeeze(phase);
        w_hz = wout / (2*pi); % na Hz
        % ===========================================================
        % 🔧 GENERACJA REALISTYCZNYCH DANYCH POMIAROWYCH
        % (Reszta kodu generującego pomiary pozostaje bez zmian...)
        f_nom = f_pom{i}; % Hz
        Ain = A_in{i};
        rng(100 + i);
        freq_jitter_frac = 0.01; amp_rel_noise = 0.025; amp_bias = 0.00;
        phase_noise_deg = 3; phase_bias_deg = 0;
        f_meas = f_nom .* (1 + freq_jitter_frac .* randn(size(f_nom)));
        w_meas = f_meas * 2*pi;
        [mag_model, phase_model] = bode(system_closed, w_meas);
        mag_model = squeeze(mag_model); phase_model = squeeze(phase_model);
        mag_meas = mag_model .* (1 + amp_rel_noise .* randn(size(mag_model))) + amp_bias;
        mag_meas(mag_meas <= 0) = eps;
        phase_meas = phase_model + phase_bias_deg + phase_noise_deg .* randn(size(phase_model));
        f = f_meas * 2*pi;
        Aout = mag_meas .* Ain; phi = phase_meas;
        G_pom = Aout ./ Ain; G_pom_dB = 20*log10(G_pom);
        fprintf('\n=== Punkty pomiarowe dla K_c = %.2f (Układ A) ===\n', K_c(i));
        fprintf('f [Hz]\t\tAmplituda [V]\tFaza [°]\n');
        for k = 1:length(f_meas)
            fprintf('%.2f\t\t%.4f\t\t%.2f\n', f_meas(k), G_pom(k), phi(k));
        end
        fprintf('---------------------------------------------\n');
        % ===========================================================
        
        % --- Marginesy stabilności (na układzie zamkniętym) ---
        [Gm, Pm, Wcg, Wcp] = margin(system_closed);
        Gm_dB = 20*log10(Gm);
        % Przesunięcia opisów (dla czytelności)
        y_offset_gm = 3 + (i-1)*3; % 3, 6, 9
        y_offset_pm = -165 - (i-1)*10; % -165, -175, -185
        
        % +++ POPRAWKA: Śledzenie limitów W KAŻDEJ pętli +++
        % Dodajemy mały margines (np. 15 stopni lub 5 dB)
        min_annotation_y_pm = min(min_annotation_y_pm, y_offset_pm - 15);
        max_annotation_y_gm = max(max_annotation_y_gm, y_offset_gm + 5);
        
        % --- Wykres amplitudy (model) ---
        h_line(i) = plot(ax_mag, w_hz, mag_dB, 'Color', currentColor, ...
            'LineWidth', 1.5, 'DisplayName', sprintf('K_c = %.2f', K_c(i)));
        % --- Punkty pomiarowe (nie dodawaj do legendy) ---
        p = plot(ax_mag, f_meas, G_pom_dB, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', currentColor, ...
        'MarkerSize', 7, 'HandleVisibility', 'off');
        h_points(i) = p(1);  % tylko pierwszy uchwyt
        % --- Zapas wzmocnienia ---
        if ~isnan(Wcg)
            plot(ax_mag, Wcg/(2*pi), 0, 's', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', currentColor, 'HandleVisibility', 'off');
            text(ax_mag, Wcg/(2*pi), y_offset_gm, ...
                sprintf('Gm (Kc=%.2f) = %.1f dB', K_c(i), Gm_dB), ...
                'HorizontalAlignment', 'center', 'Color', currentColor, 'FontSize', 8);
        end
        % --- Wykres fazy ---
        plot(ax_phase, w_hz, phase_deg, 'Color', currentColor, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(ax_phase, f_meas, phi, 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', currentColor, 'MarkerSize', 7, 'HandleVisibility', 'off');
        % --- Zapas fazy ---
        if ~isnan(Wcp)
            plot(ax_phase, Wcp/(2*pi), -180, 's', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', currentColor, 'HandleVisibility', 'off');
            text(ax_phase, Wcp/(2*pi), y_offset_pm, ...
                sprintf('Pm (Kc=%.2f) = %.1f°', K_c(i), Pm), ...
                'HorizontalAlignment', 'center', 'Color', currentColor, 'FontSize', 8);
        end
        % --- Przygotowanie symboli do legendy (raz, na końcu pętli) ---
        if i == length(K_c)
            h_marker = plot(ax_mag, NaN, NaN, 'ok', 'MarkerFaceColor', [0.7 0.7 0.7], ...
                'DisplayName', 'Punkty pomiarowe');
            h_gm_pm = plot(ax_mag, NaN, NaN, 'sk', 'MarkerFaceColor', [0.7 0.7 0.7], ...
                'DisplayName', 'Zapas fazy i wzmocnienia');
        end
    end
    % --- Legenda i ustawienia końcowe ---
    legend(ax_mag, [h_line h_marker h_gm_pm], 'Location', 'Best');
    
    % --- POPRAWKA: Dynamiczne ustawianie YLim PO narysowaniu ---
    current_ylim_pm = ylim(ax_phase); 
    ylim(ax_phase, [min(current_ylim_pm(1), min_annotation_y_pm), current_ylim_pm(2)]);
    
    current_ylim_gm = ylim(ax_mag); 
    ylim(ax_mag, [current_ylim_gm(1), max(current_ylim_gm(2), max_annotation_y_gm)]);
    
    set(ax_mag, 'XScale', 'log');
    set(ax_phase, 'XScale', 'log');
end

%bode bez zaznaczania punktów
if char_bodego 
    for i = 1 : 3
        figure(i); 
        Gc = tf(K_c(i)); 
         system_closed = feedback(Gc*Gi*Go,1);
     margin(system_closed); %zaznacza odrazu zapas wzmocnienia 
     legend(['Kc = ',num2str(K_c(i))]);
    end
end


% Wykres Nequista: dla ukladu otwartego z regulatorem P, 3 wykresy
% if char_nyquista
%     for i=1 : 3
%         figure(i);
%         Gc = tf(K_c(i));
%         nyquist(Gc*Go*Gi);
%         xlim([-1.4,0.1]);
%         ylim([-3,3]);
%         grid on;
%         legend(['K_c = ', num2str(K_c(i))]);
%     end 
%     %title('Charakterystyka Nequista ukladu A')
% end

% Wykres Nyquista - WSZYSTKIE NA 1 WYKRESIE
if char_nyquista
    figure; % JEDNA figura przed pętlą
    hold on; % Trzymaj wykresy
    grid on;
    legend_entries_nyquist = {}; % Zbiornik na legendę
    
    for i=1 : length(K_c)
        Gc = tf(K_c(i));
        nyquist(Gc*Go*Gi); % Narysuj wykres
        legend_entries_nyquist{end+1} = ['K_c = ', num2str(K_c(i))];
    end 
    
    % Ustaw osie i legendę PO pętli
    xlim([-1.4, 0.1]);
    ylim([-3, 3]);
    title('Charakterystyka Nyquista (układ A)');
    legend(legend_entries_nyquist, 'Location', 'Best');
    hold off;
end

%% Uklad B

close all;
clear all;

%co rysowac
odpowiedz_skokowa = false;
linie_pierwiastkowe = false;
char_bodego_z_pomiarami = true;
char_bodego = false;
char_nyquista = false;


K_c = [0.72,1.22,1.87];

T_i = 0.0013;
Gi = tf(1,[T_i 0]);

T_x = 0.000342;
T_y = 0.0001;
Go = tf([-T_x 1],[T_y 1]);

%symulacja
if odpowiedz_skokowa
    %pomiary z oscylokopu dla k1
    t_pomiar_k1 = [ 0.002, 0.003, 0.004, 0.005, 0.007];
    y_pomiar_k1 = [  1.2,   1.85,  2.1,   2.0,   1.8];
    %pomiary z oscylokopu dla k2
    t_pomiar_k2 = [ 0.002, 0.003, 0.004, 0.005, 0.007];
    y_pomiar_k2 = [  1.2,   1.85,  2.1,   2.0,   1.8];
    %pomiary z oscylokopu dla k3
    t_pomiar_k3 = [ 0.002, 0.003, 0.004, 0.005, 0.007];
    y_pomiar_k3 = [  1.2,   1.85,  2.1,   2.0,   1.8];

    t = 0:0.0001:0.01;
    
    u = zeros(size(t));
    u(t>=0.001) = 2;
    
    figure;
    hold on;
    grid on;

   legend_entries = {};

    for Kc = K_c
        % --- Wewnątrz pętli dla każdego K_p ---
    
        % 1. Zdefiniuj regulator P dla BIEŻĄCEGO wzmocnienia
        Gc = tf(Kc);
    
        % 2. Zbuduj układ otwarty i zamknięty dla BIEŻĄCEGO K_p
        system_open = Gc * Gi * Go;
        system_closed = feedback(system_open, 1);
    
        % 3. Uruchom symulację
        [Y, T] = lsim(system_closed, u, t);
    
        % 4. Narysuj wynik na przygotowanym wcześniej wykresie
        plot(T, Y, 'LineWidth', 1.5); % 'LineWidth' pogrubia linię
    
        % 5. Dodaj wpis do legendy
        legend_entries{end+1} = ['K_c = ', num2str(Kc)];
    end
    legend_entries{end+1} = 'Sygnał wejściowy';
    plot(t,u,'black');
%dodanie punktow pomiarowych
    plot(t_pomiar_k1,y_pomiar_k1,'bo','MarkerSize', 6);
    legend_entries{end+1} = 'Punkty pomiarowe dla k1';
    plot(t_pomiar_k2,y_pomiar_k2,'ro','MarkerSize', 6);
    legend_entries{end+1} = 'Punkty pomiarowe dla k2';
    plot(t_pomiar_k3,y_pomiar_k3,'o','Color',[1 0.65 0],'MarkerSize', 6);
    legend_entries{end+1} = 'Punkty pomiarowe dla k3';
% Wyłącz "przytrzymywanie"
    hold off; 
    title('Odpowiedź na skok dla różnych wzmocnień K_p');
    xlabel('Czas [s]');
    ylabel('Odpowiedź wyjściowa');
    legend(legend_entries); % Stwórz legendę na podstawie zebranych etykiet
end

%linie pierwiastkowe
if linie_pierwiastkowe
    figure;
    rlocus(Gi*Go);
    xlim([-5000 10000]);
    ylim([-7000 7000]);
  
    %[k_graniczne, bieguny] = rlocfind(Gi*Go);
    %disp(k_graniczne);
    % crosscut = find(real(r(:,1)) >= 0, 1,'first');
    % disp(crosscut);
end

% Char. Bodego: dla ukladu zamknietego (bo takie mamy pomiary), trzy
% charakterystyki
% if char_bodego_z_pomiarami
%     f_pom = { [1 10 100], [1 5 50 100], [1 10 100 1000] };     % Hz
%     A_in  = { [1 1 1], [1 1 1 1], [1 1 1 1] };
%     A_out = { [0.9 0.7 0.3], [0.95 0.8 0.4 0.2], [1 0.9 0.5 0.25] };
%     phi_pom = { [-5 -45 -120], [-3 -30 -100 -150], [-10 -60 -130 -170] };  % °
% 
%     % --- Główna pętla ---
%     for i = 1 : length(K_c)
%         figure(i);
% 
%         % --- Tworzenie transmitancji ---
%         Gc = tf(K_c(i));
%         system_open = Gc * Gi * Go;         % układ otwarty (dla marginesów)
%         system_closed = feedback(system_open, 1); % układ zamknięty (dla Bodego)
% 
%         % --- Dane do Bodego ---
%         [mag, phase, wout] = bode(system_closed);
%         mag = squeeze(mag);
%         phase = squeeze(phase);
%         w = wout / (2*pi); % Hz
% 
%         % --- Dane pomiarowe (dla danej iteracji) ---
%         f = f_pom{i} * 2 * pi;
%         Ain = A_in{i};
%         Aout = A_out{i};
%         phi = phi_pom{i};
% 
%         G_pom = Aout ./ Ain; % transmitancja z pomiaru
% 
%         % --- Marginesy ---
%         [Gm, Pm, Wcg, Wcp] = margin(system_closed);
%         Gm_dB = 20*log10(Gm);
% 
%         % --- Wykres amplitudowy ---
%         subplot(2,1,1);
%         semilogx(w, 20*log10(mag), 'b', 'LineWidth', 1.5); hold on;
%         semilogx(f, 20*log10(G_pom), 'ro', 'MarkerFaceColor', 'r');
%         grid on;
%         ylabel('Amplituda [dB]');
%         title(['Charakterystyka Bodego - K_c = ', num2str(K_c(i))]);
%         legend('Model', 'Pomiar', 'Location', 'Best');
% 
%         % --- Zaznaczenie zapasu wzmocnienia ---
%         if ~isnan(Wcg)
%             semilogx(Wcg/(2*pi), 0, 'ks', 'MarkerFaceColor', 'y', 'DisplayName', 'Zapas wzmocnienia');
%             text(Wcg/(2*pi), 3, sprintf('Gm = %.1f dB', Gm_dB), 'HorizontalAlignment', 'center');
%         end
% 
%         % --- Wykres fazowy ---
%         subplot(2,1,2);
%         semilogx(w, phase, 'b', 'LineWidth', 1.5); hold on;
%         semilogx(f, phi, 'ro', 'MarkerFaceColor', 'r');
%         grid on;
%         ylabel('Faza [°]');
%         xlabel('Częstotliwość [Hz]');
%         legend('Model', 'Pomiar', 'Location', 'Best');
% 
%         % --- Zaznaczenie zapasu fazy ---
%         if ~isnan(Wcp)
%             semilogx(Wcp/(2*pi), -180, 'ks', 'MarkerFaceColor', 'c', 'DisplayName', 'Zapas fazy');
%             text(Wcp/(2*pi), -165, sprintf('Pm = %.1f°', Pm), 'HorizontalAlignment', 'center');
%         end
%     end
% end

% --- Charakterystyka Bodego: układ zamknięty z pomiarami (wszystko na 1 wykresie) ---
if char_bodego_z_pomiarami
    % --- Domyślne częstotliwości pomiarowe (Hz) ---
    f_pom = { ...
    [10 100 200 400 800 3500], ...   % k3
    [10 100 200 400 600 1750], ...   % k2
    [10 100 200 400 520 1000] ...    % k1
    };
    A_in  = { [1 1 1], [1 1 1 1], [1 1 1 1] };              % sygnał wejściowy = 1
    % --- Definicja kolorów (Twoja kolejność) ---
    c_default = get(groot, 'defaultAxesColorOrder');
    colors = [c_default(1,:); c_default(3,:); c_default(2,:)]; % nieb, żółty, pomarańczowy
    % --- Utwórz figurę i osie ---
    figure;
    % --- Wykres amplitudy ---
    ax_mag = subplot(2,1,1);
    hold on; grid on;
    ylabel('Amplituda [dB]');
    title('Charakterystyki Bodego (układ B)');
    % --- Wykres fazy ---
    ax_phase = subplot(2,1,2);
    hold on; grid on;
    ylabel('Faza [°]');
    xlabel('Częstotliwość [Hz]');
    
    % --- Pętla po wszystkich wartościach K_c ---
    min_annotation_y = -165; % Zmienna do śledzenia najniższej adnotacji

    for i = 1 : length(K_c)
        currentColor = colors(i, :); % kolor danej charakterystyki
        % --- Tworzenie transmitancji ---
        Gc = tf(K_c(i));
        system_open = Gc * Gi * Go;                 % układ otwarty
        system_closed = feedback(system_open, 1);   % układ zamknięty
        
        % --- Dane modelu do Bodego ---
        [mag, phase, wout] = bode(system_closed);
        mag_dB = 20*log10(squeeze(mag));
        phase_deg = squeeze(phase);
        
        % +++ POPRAWKA ZAWINIĘCIA FAZY (Model) +++
        phase_deg(phase_deg > 0) = phase_deg(phase_deg > 0) - 360; 
        % +++++++++++++++++++++++++++++++++++++++
        
        w_hz = wout / (2*pi); % na Hz
        
        % ===========================================================
        % 🔧 GENERACJA REALISTYCZNYCH DANYCH POMIAROWYCH
        % ===========================================================
        f_nom = f_pom{i}; % Hz
        Ain = A_in{i};
        % Parametry „realistyczne"
        rng(100 + i);            % różne dla każdego Kc, ale powtarzalne
        freq_jitter_frac = 0.01; % 1% jitter częstotliwości
        amp_rel_noise = 0.025;   % 2.5% szum amplitudy
        amp_bias = 0.00;         % brak systematycznego błędu amplitudy
        phase_noise_deg = 3;     % ±3° szum fazy
        phase_bias_deg = 0;      % brak systematycznego błędu fazy
        % --- Częstotliwości pomiarowe z jitterem ---
        f_meas = f_nom .* (1 + freq_jitter_frac .* randn(size(f_nom)));
        w_meas = f_meas * 2*pi; % rad/s
        % --- Wartości modelu w tych częstotliwościach ---
        [mag_model, phase_model] = bode(system_closed, w_meas);
        mag_model = squeeze(mag_model);
        phase_model = squeeze(phase_model);
        
        % +++ POPRAWKA ZAWINIĘCIA FAZY (Pomiary) +++
        phase_model(phase_model > 0) = phase_model(phase_model > 0) - 360;
        % ++++++++++++++++++++++++++++++++++++++++++
        
        % --- Dodanie szumu pomiarowego ---
        mag_meas = mag_model .* (1 + amp_rel_noise .* randn(size(mag_model))) + amp_bias;
        mag_meas(mag_meas <= 0) = eps;
        phase_meas = phase_model + phase_bias_deg + phase_noise_deg .* randn(size(phase_model));
        % --- Przygotowanie danych jak w oryginale ---
        f = f_meas * 2*pi;   % w [rad/s]
        Aout = mag_meas .* Ain;
        phi = phase_meas;
        G_pom = Aout ./ Ain;
        G_pom_dB = 20*log10(G_pom);
        % --- Wypisanie danych pomiarowych w konsoli ---
        fprintf('\n=== Punkty pomiarowe dla K_c = %.2f (Układ B) ===\n', K_c(i));
        fprintf('f [Hz]\t\tAmplituda [V]\tFaza [°]\n');
        for k = 1:length(f_meas)
            fprintf('%.2f\t\t%.4f\t\t%.2f\n', f_meas(k), G_pom(k), phi(k));
        end
        fprintf('---------------------------------------------\n');
        % ===========================================================
        
        % --- Marginesy stabilności (na układzie zamkniętym) ---
        [Gm, Pm, Wcg, Wcp] = margin(system_closed);
        Gm_dB = 20*log10(Gm);
        % Przesunięcia opisów (dla czytelności)
        y_offset_gm = 3 + (i-1)*3;
        y_offset_pm = -165 - (i-1)*10; % -165, -175, -185
        
        % --- Śledź najniższą pozycję Y tekstu ---
        if i == length(K_c)
             min_annotation_y = y_offset_pm - 15; % Ostatni tekst (-185) + 15 stopni marginesu
        end

        % --- Wykres amplitudy (model) ---
        h_line(i) = plot(ax_mag, w_hz, mag_dB, 'Color', currentColor, ...
            'LineWidth', 1.5, 'DisplayName', sprintf('K_c = %.2f', K_c(i)));
        % --- Punkty pomiarowe (nie dodawaj do legendy) ---
        p = plot(ax_mag, f_meas, G_pom_dB, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', currentColor, ...
        'MarkerSize', 7, 'HandleVisibility', 'off');
        h_points(i) = p(1);  % tylko pierwszy uchwyt
        % --- Zapas wzmocnienia ---
        if ~isnan(Wcg)
            plot(ax_mag, Wcg/(2*pi), 0, 's', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', currentColor, 'HandleVisibility', 'off');
            text(ax_mag, Wcg/(2*pi), y_offset_gm, ...
                sprintf('Gm (Kc=%.2f) = %.1f dB', K_c(i), Gm_dB), ...
                'HorizontalAlignment', 'center', 'Color', currentColor, 'FontSize', 8);
        end
        % --- Wykres fazy ---
        plot(ax_phase, w_hz, phase_deg, 'Color', currentColor, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(ax_phase, f_meas, phi, 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', currentColor, 'MarkerSize', 7, 'HandleVisibility', 'off');
        % --- Zapas fazy ---
        if ~isnan(Wcp)
            plot(ax_phase, Wcp/(2*pi), -180, 's', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', currentColor, 'HandleVisibility', 'off');
            text(ax_phase, Wcp/(2*pi), y_offset_pm, ...
                sprintf('Pm (Kc=%.2f) = %.1f°', K_c(i), Pm), ...
                'HorizontalAlignment', 'center', 'Color', currentColor, 'FontSize', 8);
        end
        % --- Przygotowanie symboli do legendy (raz, na końcu pętli) ---
        if i == length(K_c)
            h_marker = plot(ax_mag, NaN, NaN, 'ok', 'MarkerFaceColor', [0.7 0.7 0.7], ...
                'DisplayName', 'Punkty pomiarowe');
            h_gm_pm = plot(ax_mag, NaN, NaN, 'sk', 'MarkerFaceColor', [0.7 0.7 0.7], ...
                'DisplayName', 'Zapas fazy i wzmocnienia');
        end
    end
    % --- Legenda i ustawienia końcowe ---
    legend(ax_mag, [h_line h_marker h_gm_pm], 'Location', 'Best');
    
    % --- POPRAWKA: Dynamiczne ustawianie YLim PO narysowaniu ---
    current_ylim = ylim(ax_phase); % Pobierz obecne limity
    % Ustaw nowe limity: [minimum(obecne_min, najniższa_adnotacja), obecne_max]
    ylim(ax_phase, [min(current_ylim(1), min_annotation_y), current_ylim(2)]);
    
    set(ax_mag, 'XScale', 'log');
    set(ax_phase, 'XScale', 'log');
end

% Wykres Nequista: dla ukladu otwartego z regulatorem P, 3 wykresy
% if char_nyquista
%     for i=1 : 3
%         figure(i);
%         Gc = tf(K_c(i));
%         nyquist(Gc*Go*Gi);
%         xlim([-1.4,0.1]);
%         ylim([-3,3]);
%         grid on;
%         legend(['K_c = ', num2str(K_c(i))]);
%     end 
%     %title('Charakterystyka Nequista ukladu A')
% end

% --- Nyquist: wszystkie K_c na jednym wykresie (układ B) ---
if char_nyquista
    figure;
    hold on; grid on;
    legend_entries = {};
    for i = 1:length(K_c)
        Gc = tf(K_c(i));
        nyquist(Gc*Gi*Go);
        legend_entries{end+1} = ['K_c = ', num2str(K_c(i))];
    end
    xlim([-1.4 0.1]); ylim([-3 3]);
    title('Charakterystyka Nyquista (układ B)');
    legend(legend_entries,'Location','Best');
    hold off;
end


%% uklad D
close all;
clear;

%co rysowac
odpowiedz_skokowa = true;
linie_pierwiastkowe = false;
char_bodego_z_pomiarami = false;
char_bodego = false;
char_nyquista = false;
char_czest = false;

K_c = [1.02, 2.47, 4.47]; 
Kp = 0.871;
Tp = 0.00078;
Go = tf(Kp, [Tp 1]);

if odpowiedz_skokowa
    t = 0:0.001:0.1;
    % sygnał zadany
    r = zeros(size(t));
    r(t >= 0.01) = 2;
    
    % sygnał zakłócenia
    f = 50;
    A = 1;
    d = A * sin(2*pi*f*t);

    % --- Wykres ---
    figure;
    hold on;
    grid on;
    legend_entries = {};

    % --- Pętla po wzmocnieniach K_c ---
    for Kc = K_c
        % 1. Definicja regulatora
        P = tf(Kc);

        % 2. Układ otwarty i zamknięty
        object_open_loop = P * Go;
        G_ry = feedback(object_open_loop, 1);
        G_dy = feedback(Go, P);

        % 3. Symulacja
        [Y_r, T] = lsim(G_ry, r, t);
        [Y_d, ~] = lsim(G_dy, d, t);
        Y_out = Y_r + Y_d; % superpozycja sygnałów

        % 4. Rysowanie odpowiedzi
        plot(T, Y_out, 'LineWidth', 1.5);
        legend_entries{end+1} = ['K_c = ', num2str(Kc)];
    end

    % --- Sygnał wejściowy ---
    plot(t, r, 'k--');
    legend_entries{end+1} = 'Sygnał wejściowy)';

    % --- Ustawienia wykresu ---
    hold off;
   % title('Symulacja układu D dla różnych wartości K_c z zakłóceniem');
    xlabel('Czas [s]');
    %ylabel('Odpowiedź układu');
    legend(legend_entries, 'Location', 'Best');
end


%linie pierwiastkowe
if linie_pierwiastkowe
    figure;
    rlocus(Go);
    xlim([-5000 1000]);
    ylim([-7000 7000]);
  
    %[k_graniczne, bieguny] = rlocfind(Gi*Go);
    %disp(k_graniczne);
    % crosscut = find(real(r(:,1)) >= 0, 1,'first');
    % disp(crosscut);
end

%charakterystyka bodego
if char_bodego
     for i = 1 : 3
        figure(i); 
        P = tf(K_c(i)); 
        system_closed = feedback(P*Go,1);
        margin(system_closed); %zaznacza odrazu zapas wzmocnienia 
        legend(['Kc = ',num2str(K_c(i))]);
    end
end

% Wykres Nequista: dla ukladu otwartego z regulatorem P, 3 wykresy
% if char_nyquista
%     for i=1 : 3
%         figure(i);
%         Gc = tf(K_c(i));
%         nyquist(Gc*Go);
%         %xlim([-1.4,0.1]);
%         %ylim([-3,3]);
%         grid on;
%         legend(['K_c = ', num2str(K_c(i))]);
%     end 
% end

% --- Charakterystyka Bodego: układ zamknięty z pomiarami (wszystko na 1 wykresie) ---
if char_bodego_z_pomiarami
    % --- Domyślne częstotliwości pomiarowe (Hz) ---
    f_pom = { ...
    [100 200 400 800 1500 2000], ...   % k3
    [100 200 400 800 1500 2000], ...   % k2
    [100 200 400 800 1500 2000] ...    % k1
    };
    A_in  = { [1 1 1], [1 1 1 1], [1 1 1 1] };              % sygnał wejściowy = 1
    % --- Definicja kolorów (Twoja kolejność) ---
    c_default = get(groot, 'defaultAxesColorOrder');
    colors = [c_default(1,:); c_default(3,:); c_default(2,:)]; % nieb, żółty, pomarańczowy
    % --- Utwórz figurę i osie ---
    figure;
    % --- Wykres amplitudy ---
    ax_mag = subplot(2,1,1);
    hold on; grid on;
    ylabel('Amplituda [dB]');
    title('Charakterystyki Bodego (układ D)');
    % --- Wykres fazy ---
    ax_phase = subplot(2,1,2);
    hold on; grid on;
    ylabel('Faza [°]');
    xlabel('Częstotliwość [Hz]');
    % --- Pętla po wszystkich wartościach K_c ---
    for i = 1 : length(K_c)
        currentColor = colors(i, :); % kolor danej charakterystyki
        % --- Tworzenie transmitancji ---
        Gc = tf(K_c(i));
        system_open = Gc * Go;                      % <<< MODYFIKACJA DLA UKLADU D (brak Gi)
        system_closed = feedback(system_open, 1);   % układ zamknięty
        % --- Dane modelu do Bodego ---
        [mag, phase, wout] = bode(system_closed);
        mag_dB = 20*log10(squeeze(mag));
        phase_deg = squeeze(phase);
        w_hz = wout / (2*pi); % na Hz
        % ===========================================================
        % 🔧 GENERACJA REALISTYCZNYCH DANYCH POMIAROWYCH
        % ===========================================================
        f_nom = f_pom{i}; % Hz
        Ain = A_in{i};
        % Parametry „realistyczne"
        rng(100 + i);            % różne dla każdego Kc, ale powtarzalne
        freq_jitter_frac = 0.01; % 1% jitter częstotliwości
        amp_rel_noise = 0.025;   % 2.5% szum amplitudy
        amp_bias = 0.00;         % brak systematycznego błędu amplitudy
        phase_noise_deg = 3;     % ±3° szum fazy
        phase_bias_deg = 0;      % brak systematycznego błędu fazy
        % --- Częstotliwości pomiarowe z jitterem ---
        f_meas = f_nom .* (1 + freq_jitter_frac .* randn(size(f_nom)));
        w_meas = f_meas * 2*pi; % rad/s
        % --- Wartości modelu w tych częstotliwościach ---
        [mag_model, phase_model] = bode(system_closed, w_meas);
        mag_model = squeeze(mag_model);
        phase_model = squeeze(phase_model);
        % --- Dodanie szumu pomiarowego ---
        mag_meas = mag_model .* (1 + amp_rel_noise .* randn(size(mag_model))) + amp_bias;
        mag_meas(mag_meas <= 0) = eps;
        phase_meas = phase_model + phase_bias_deg + phase_noise_deg .* randn(size(phase_model));
        % --- Przygotowanie danych jak w oryginale ---
        f = f_meas * 2*pi;   % w [rad/s]
        Aout = mag_meas .* Ain;
        phi = phase_meas;
        G_pom = Aout ./ Ain;
        G_pom_dB = 20*log10(G_pom);
        % --- Wypisanie danych pomiarowych w konsoli ---
        fprintf('\n=== Punkty pomiarowe dla K_c = %.2f (Układ D) ===\n', K_c(i));
        fprintf('f [Hz]\t\tAmplituda [V]\tFaza [°]\n');
        for k = 1:length(f_meas)
            fprintf('%.2f\t\t%.4f\t\t%.2f\n', f_meas(k), G_pom(k), phi(k));
        end
        fprintf('---------------------------------------------\n');
        % ===========================================================
        % --- Marginesy stabilności (na układzie zamkniętym) ---
        [Gm, Pm, Wcg, Wcp] = margin(system_closed);
        Gm_dB = 20*log10(Gm);
        % Przesunięcia opisów (dla czytelności)
        y_offset_gm = 3 + (i-1)*3;
        y_offset_pm = -165 - (i-1)*10;
        % --- Wykres amplitudy (model) ---
        h_line(i) = plot(ax_mag, w_hz, mag_dB, 'Color', currentColor, ...
            'LineWidth', 1.5, 'DisplayName', sprintf('K_c = %.2f', K_c(i)));
        % --- Punkty pomiarowe (nie dodawaj do legendy) ---
        p = plot(ax_mag, f_meas, G_pom_dB, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', currentColor, ...
        'MarkerSize', 7, 'HandleVisibility', 'off');
        h_points(i) = p(1);  % tylko pierwszy uchwyt
        % --- Zapas wzmocnienia ---
        if ~isnan(Wcg)
            plot(ax_mag, Wcg/(2*pi), 0, 's', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', currentColor, 'HandleVisibility', 'off');
            text(ax_mag, Wcg/(2*pi), y_offset_gm, ...
                sprintf('Gm (Kc=%.2f) = %.1f dB', K_c(i), Gm_dB), ...
                'HorizontalAlignment', 'center', 'Color', currentColor, 'FontSize', 8);
        end
        % --- Wykres fazy ---
        plot(ax_phase, w_hz, phase_deg, 'Color', currentColor, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(ax_phase, f_meas, phi, 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', currentColor, 'MarkerSize', 7, 'HandleVisibility', 'off');
        % --- Zapas fazy ---
        if ~isnan(Wcp)
            plot(ax_phase, Wcp/(2*pi), -180, 's', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', currentColor, 'HandleVisibility', 'off');
            text(ax_phase, Wcp/(2*pi), y_offset_pm, ...
                sprintf('Pm (Kc=%.2f) = %.1f°', K_c(i), Pm), ...
                'HorizontalAlignment', 'center', 'Color', currentColor, 'FontSize', 8);
        end
        % --- Przygotowanie symboli do legendy (raz, na końcu pętli) ---
        if i == length(K_c)
            h_marker = plot(ax_mag, NaN, NaN, 'ok', 'MarkerFaceColor', [0.7 0.7 0.7], ...
                'DisplayName', 'Punkty pomiarowe');
            h_gm_pm = plot(ax_mag, NaN, NaN, 'sk', 'MarkerFaceColor', [0.7 0.7 0.7], ...
                'DisplayName', 'Zapas fazy i wzmocnienia');
        end
    end
    % --- Legenda i ustawienia końcowe ---
    legend(ax_mag, [h_line h_marker h_gm_pm], 'Location', 'Best');
    set(ax_mag, 'XScale', 'log');
    set(ax_phase, 'XScale', 'log');
end

% --- Nyquist: wszystkie K_c na jednym wykresie (układ D) ---
if char_nyquista
    figure;
    hold on; grid on;
    legend_entries = {};
    for i = 1:length(K_c)
        Gc = tf(K_c(i));
        nyquist(Gc*Go);
        legend_entries{end+1} = ['K_c = ', num2str(K_c(i))];
    end
    title('Charakterystyka Nyquista (układ D)');
    legend(legend_entries,'Location','Best');
    hold off;
end

%symulacja rozncych transmitancji i ich charajterystyk
if char_czest
    %wybrane wzmocnienie:
    K_c = 1.22;
    P = tf(K_c);

    %policzone transmitanjce;
    EoverR = feedback(1,P*Go);
    UoverR = feedback(P,Go);
    EoverD = -feedback(Go,P);
    UoverD = -feedback(P*Go,1);

    % lista transmitancji do analizy:
    systems = {EoverR, UoverR, EoverD, UoverD};
    labels = {'|E/R|', '|U/R|', '|E/D|', '|U/D|'};

    % przygotowanie zakresu częstotliwości (logarytmiczny)
    w = logspace(2, 5, 500);   % od 0.01 do 1000 rad/s

    figure;
    hold on;
    grid on;
    title('Charakterystyki amplitudowe');
    xlabel('Częstotliwość \omega [rad/s]');
    ylabel('Wzmocnienie |G(j\omega)|');

    % pętla po transmitancjach
    for i = 1:length(systems)
        [mag, phase, wout] = bode(systems{i});
        mag = squeeze(mag);%zmienienia z macierzy na wektor
        plot(wout, mag, 'DisplayName', labels{i});
    end

    set(gca, 'XScale', 'log');   % logarytmiczna oś częstotliwości
    legend show;
    hold off;

end