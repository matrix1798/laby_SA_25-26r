%% Uklad A: 

close all;
clear all;

K_c = [0.52 1.12 1.67];

T_i = 0.0013;
Gi = tf(1,[T_i 0]); %czlon calkujacy

w0 = 2560;
zeta = 0.37;
Go = tf(w0^2,[1 2*zeta*w0 w0^2]);

%odpowiedz ukladu
if false
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
if false
    figure;
    rlocus(Gi*Go);
    %title('Linie pieriwastkowe ukladu A');
    xlim([-3000 2000]);
    ylim([-5000 5000]);
    
  [k_graniczne, bieguny] = rlocfind(Go*Gc);
  disp(k_graniczne);
   % crosscut = find(real(r(:,1)) >= 0, 1,'first');
   % disp(crosscut);
end

% Char. Bodego: dla ukladu zamknietego (bo takie mamy pomiary), trzy
% charakterystyki
if false
    f_pom = { [1 10 100], [1 5 50 100], [1 10 100 1000] };     % Hz
    A_in  = { [1 1 1], [1 1 1 1], [1 1 1 1] };
    A_out = { [0.9 0.7 0.3], [0.95 0.8 0.4 0.2], [1 0.9 0.5 0.25] };
    phi_pom = { [-5 -45 -120], [-3 -30 -100 -150], [-10 -60 -130 -170] };  % °

    % --- Główna pętla ---
    for i = 1 : length(K_c)
        figure(i);
        
        % --- Tworzenie transmitancji ---
        Gc = tf(K_c(i));
        system_open = Gc * Gi * Go;         % układ otwarty (dla marginesów)
        system_closed = feedback(system_open, 1); % układ zamknięty (dla Bodego)
        
        % --- Dane do Bodego ---
        [mag, phase, wout] = bode(system_closed);
        mag = squeeze(mag);
        phase = squeeze(phase);
        w = wout / (2*pi); % Hz
        
        % --- Dane pomiarowe (dla danej iteracji) ---
        f = f_pom{i}*2*pi;
        Ain = A_in{i};
        Aout = A_out{i};
        phi = phi_pom{i};
        
        G_pom = Aout ./ Ain; % transmitancja z pomiaru
        
        % --- Marginesy ---
        [Gm, Pm, Wcg, Wcp] = margin(system_closed);
        Gm_dB = 20*log10(Gm);
        
        % --- Wykres amplitudowy ---
        subplot(2,1,1);
        semilogx(w, 20*log10(mag), 'b', 'LineWidth', 1.5); hold on;
        semilogx(f, 20*log10(G_pom), 'ro', 'MarkerFaceColor', 'r');
        grid on;
        ylabel('Amplituda [dB]');
        title(['Charakterystyka Bodego - K_c = ', num2str(K_c(i))]);
        legend('Model', 'Pomiar', 'Location', 'Best');
        
        % --- Zaznaczenie zapasu wzmocnienia ---
        if ~isnan(Wcg)
            semilogx(Wcg/(2*pi), 0, 'ks', 'MarkerFaceColor', 'y', 'DisplayName', 'Zapas wzmocnienia');
            text(Wcg/(2*pi), 3, sprintf('Gm = %.1f dB', Gm_dB), 'HorizontalAlignment', 'center');
        end
        
        % --- Wykres fazowy ---
        subplot(2,1,2);
        semilogx(w, phase, 'b', 'LineWidth', 1.5); hold on;
        semilogx(f, phi, 'ro', 'MarkerFaceColor', 'r');
        grid on;
        ylabel('Faza [°]');
        xlabel('Częstotliwość [Hz]');
        legend('Model', 'Pomiar', 'Location', 'Best');
        
        % --- Zaznaczenie zapasu fazy ---
        if ~isnan(Wcp)
            semilogx(Wcp/(2*pi), -180, 'ks', 'MarkerFaceColor', 'c', 'DisplayName', 'Zapas fazy');
            text(Wcp/(2*pi), -165, sprintf('Pm = %.1f°', Pm), 'HorizontalAlignment', 'center');
        end
    end
end

%bode bez zaznaczania punktów
if false 
    for i = 1 : 3
        figure(i); 
        Gc = tf(K_c(i)); 
         system_closed = feedback(Gc*Gi*Go,1);
     margin(system_closed); %zaznacza odrazu zapas wzmocnienia 
     legend(['Kc = ',num2str(K_c(i))]);
    end
end

% Wykres Nequista: dla ukladu otwartego z regulatorem P, 3 wykresy
if true
    for i=1 : 3
        figure(i);
        Gc = tf(K_c(i));
        nyquist(Gc*Go*Gi);
        xlim([-1.4,0.1]);
        ylim([-3,3]);
        grid on;
        legend(['K_c = ', num2str(K_c(i))]);
    end 
    %title('Charakterystyka Nequista ukladu A')
end

%% Uklad B

K_c = [0.47,1,2.1];

T_i = 0.0013;
Gi = tf(1,[T_i 0]);

T_x = 0.000342;
T_y = 0.0001;
Go = tf([-T_x 1],[T_y 1]);

%symulacja
if false
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

% Char. Bodego: dla ukladu zamknietego (bo takie mamy pomiary), trzy
% charakterystyki
if true
    f_pom = { [1 10 100], [1 5 50 100], [1 10 100 1000] };     % Hz
    A_in  = { [1 1 1], [1 1 1 1], [1 1 1 1] };
    A_out = { [0.9 0.7 0.3], [0.95 0.8 0.4 0.2], [1 0.9 0.5 0.25] };
    phi_pom = { [-5 -45 -120], [-3 -30 -100 -150], [-10 -60 -130 -170] };  % °

    % --- Główna pętla ---
    for i = 1 : length(K_c)
        figure(i);
        
        % --- Tworzenie transmitancji ---
        Gc = tf(K_c(i));
        system_open = Gc * Gi * Go;         % układ otwarty (dla marginesów)
        system_closed = feedback(system_open, 1); % układ zamknięty (dla Bodego)
        
        % --- Dane do Bodego ---
        [mag, phase, wout] = bode(system_closed);
        mag = squeeze(mag);
        phase = squeeze(phase);
        w = wout / (2*pi); % Hz
        
        % --- Dane pomiarowe (dla danej iteracji) ---
        f = f_pom{i} * 2 * pi;
        Ain = A_in{i};
        Aout = A_out{i};
        phi = phi_pom{i};
        
        G_pom = Aout ./ Ain; % transmitancja z pomiaru
        
        % --- Marginesy ---
        [Gm, Pm, Wcg, Wcp] = margin(system_closed);
        Gm_dB = 20*log10(Gm);
        
        % --- Wykres amplitudowy ---
        subplot(2,1,1);
        semilogx(w, 20*log10(mag), 'b', 'LineWidth', 1.5); hold on;
        semilogx(f, 20*log10(G_pom), 'ro', 'MarkerFaceColor', 'r');
        grid on;
        ylabel('Amplituda [dB]');
        title(['Charakterystyka Bodego - K_c = ', num2str(K_c(i))]);
        legend('Model', 'Pomiar', 'Location', 'Best');
        
        % --- Zaznaczenie zapasu wzmocnienia ---
        if ~isnan(Wcg)
            semilogx(Wcg/(2*pi), 0, 'ks', 'MarkerFaceColor', 'y', 'DisplayName', 'Zapas wzmocnienia');
            text(Wcg/(2*pi), 3, sprintf('Gm = %.1f dB', Gm_dB), 'HorizontalAlignment', 'center');
        end
        
        % --- Wykres fazowy ---
        subplot(2,1,2);
        semilogx(w, phase, 'b', 'LineWidth', 1.5); hold on;
        semilogx(f, phi, 'ro', 'MarkerFaceColor', 'r');
        grid on;
        ylabel('Faza [°]');
        xlabel('Częstotliwość [Hz]');
        legend('Model', 'Pomiar', 'Location', 'Best');
        
        % --- Zaznaczenie zapasu fazy ---
        if ~isnan(Wcp)
            semilogx(Wcp/(2*pi), -180, 'ks', 'MarkerFaceColor', 'c', 'DisplayName', 'Zapas fazy');
            text(Wcp/(2*pi), -165, sprintf('Pm = %.1f°', Pm), 'HorizontalAlignment', 'center');
        end
    end
end

%bode bez zaznaczania punktów
if false 
    for i = 1 : 3
        figure(i); 
        Gc = tf(K_c(i)); 
        system_closed = feedback(Gc*Gi*Go,1);
        margin(system_closed); %zaznacza odrazu zapas wzmocnienia 
        legend(['Kc = ',num2str(K_c(i))]);
    end
end

% Wykres Nequista: dla ukladu otwartego z regulatorem P, 3 wykresy
if false
    for i=1 : 3
        figure(i);
        Gc = tf(K_c(i));
        nyquist(Gc*Go*Gi);
        xlim([-1.4,0.1]);
        ylim([-3,3]);
        grid on;
        legend(['K_c = ', num2str(K_c(i))]);
    end 
    %title('Charakterystyka Nequista ukladu A')
end

%% uklad D

K_p = 1.47;
P = tf(K_p);

K = 1;
Tp = 0.005;
Go = tf(K,[Tp 1]);

object_open_loop = P * Go;

G_ry = feedback(object_open_loop,1);
G_dy = feedback(Go,P);

%symulacja
t = 0:0.001:0.1;

r = zeros(size(t));
r(t>=0.01) = 1;

f = 50;
A = 0.1;
d = A*sin(2*pi*f*t);

[Y_r,T] = lsim(G_ry,r,t);
[Y_d,~] = lsim(G_dy,d,t);

Y_out = Y_r + Y_d; %zgodnie z zasada superpozycji

figure;
plot(T,Y_out,'r','LineWidth',1);
hold on;
plot(T,r,'b');
title('Symulacja ukladu D z zakłóceniami')
xlabel('Czas');
ylabel('Amplituda');
grid on;
