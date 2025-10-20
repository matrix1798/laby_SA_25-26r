%% Uklad A: 

close all
clear all

K_p = 1.12;
G_p = tf(K_p); %regulator P

T_i = 0.0013;
Gc = tf(1,[T_i 0]); %czlon calkujacy

k = 1;
T = 0.000385; % w = 2*pi/T
w0 = 2560;
zeta = 0.37;
%Go = tf(k,[T^2,2*zeta*T,1]); %uklad 2 rzedu
Go = tf(w0^2,[1 2*zeta*w0 w0^2]);
system = Gc * Go * G_p; %uklad otwarty

object_closed = feedback(system,1); %feedback(G,H) : G/(1 + GH)

%odpowiedz ukladu
if false
    t = 0:0.0001:0.01;
    
    u = zeros(size(t));
    u(t>=0) = 2 ; %definicja skoku
    
    [Y,T] = lsim(object_closed,u,t);
    
    figure;
    plot(T,Y);
    grid on;
    title('Odpwoeidz na skok');
    xlabel('Time (s)');
    ylabel('Output Response');
    legend('Sygnal wyjsciowy');
end

% Linie pierwiastkowe:
if true
    figure;
    rlocus(Go*Gc);
    %title('Linie pieriwastkowe ukladu A');
    xlim([-3000 2000]);
    ylim([-5000 5000]);
    
  [k_graniczne, bieguny] = rlocfind(Go*Gc);
  disp(k_graniczne);
   % crosscut = find(real(r(:,1)) >= 0, 1,'first');
   % disp(crosscut);
end

% Char. Bodego:
if false
    figure;
    margin(Gc*Go); %zaznacza odrazu zapas wzmocnienia
    title('Char Bodego dla ukladu 2 rzedu')
end

% Wykres Nequista:
if false
    figure;
    nyquist(Go*Gc);
    axis equal;
    title('Charakterystyka Nequista ukladu A')
end

%% Uklad B

K_p = 3.8;
Gp = tf(K_p);

T_i = 0.0013;
Gi = tf(1,[T_i 0]);

T_x = 0.000342;
T_y = 0.0001;
Go = tf([-T_x 1],[T_y 1]);

object_open_loop = Gp * Gi * Go;

object_closed_loop = feedback(object_open_loop,1);

Gz = object_closed_loop;

%symulacja

t = 0:0.001:0.1;

u = zeros(size(t));
u(t>=0.01) = 1;

[output,T] = lsim(Gz,u,t);

figure;
plot(T, output);
grid on;
xlabel('time');
ylabel('amplitude');
legend('Odpowiedz ukladu');


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
