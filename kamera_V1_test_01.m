clc
clear all
close all
%% specyfikacja testu
%maksymalny czas pobierania danych
MaxTime = 20 ; 

% parametry filtru kalmana
% b�edy modelu
e_a_max = 1; % maksymalne przy�pieszenie kulki [cm / s^2]
e_s_max = 1; % maksymalna zmiana po�o�enia zwiazana z zewn�trznym wymuszeniem ruchu kulki [cm/s]

e_pomiar_max = 1; %maksymalny b��d pomiaru pi�eczki

%% preparation

i=2; % znacznik numeru pobranego obrazu

%% measurement (test)
% inicjalizacja czasu
time(1) = 0;
dt = 0.1;
bc(1,1) = randp(1, 0, 1); 
bc(1,2) = randp(1, 0, 1); 
while(time(i - 1) < MaxTime) %predefiniowany czas pomiaru
    dt = randp(1,0.1,0.2);
    time(i) = time(i-1) + dt;
    
    bc(i,1) = bc(i - 1, 1) + randp(1, -e_s_max, e_s_max) * dt;
    bc(i,2) = bc(i - 1, 2) + randp(1, -e_s_max, e_s_max) * dt;

    
    
    
    i = i + 1;
end
%% obr�bka danych
%pomiar
x_pomiar = bc(:,1);
y_pomiar = bc(:,2);
%pozycja pocz�tkowa
X(1) = x_pomiar(1);
Y(1) = y_pomiar(1);
%pr�dkosc pocz�tkowa
Vx(1) = 0;
Vy(1) = 0;
% inicjalizacja macierzy stanu
x_state = [X(1); Vx(1)];
y_state = [Y(1); Vy(1)];
% inicjalizacja macierzy kowariancji
Px = zeros(2,2);
Py = zeros(2,2);
Q = zeros(2,2);
R = e_pomiar_max^2; %mo�na przyjac za sta�e mo�na uzale�ni� od dt
B = [0; 0];         % macierz sta�a w czasie
C = [1 0];          % macierz sta�a w czasie
% model matematyczny uk��du
% x(i+1) = 1 * x(i) + dt * Vx(i) + dt * e_s(i)
% Vx(i+1)= 0 * x(i) + 1 *  Vx(i) + dt * e_a(i)
% x_pomiar(i+1) = x(i+1) + e_pomiar
% 
% Q = [e_s_max^2 * dt^2 0
%      0                e_a_max^2 * dt^2]
% R = e_pomiar_max^2 

% filtracja klamana

for i = 2 : length(time)
    dt = time(i) - time(i-1); % model nie musi by� z sta�ymi przyrostami czasu
    Q = [e_s_max^2 * dt^2 0
        0 e_a_max^2 * dt^2];
    A = [1 dt
        0 1];
    % mo�na rozpatrywa� fitlracje x i y osobno przyjmuj�c ten sam model
    [x_state, Px, K] = KalmanMS(A,B,C,R,Q,Px,0,x_state,x_pomiar(i));
    [y_state, Py] = KalmanMS(A,B,C,R,Q,Py,0,y_state,y_pomiar(i));
    
    X(i) = x_state(1);
    Vx(i) = x_state(2);
    Y(i) = y_state(1);
    Vy(i) = y_state(2);
    
    Vx_kalk(i) = (X(i) - X(i - 1)) / dt;
    Vy_kalk(i) = (Y(i) - Y(i - 1)) / dt;
    
    
end
figure(1);
plot(x_pomiar, y_pomiar, 'r'); grid on; hold on
plot(X,Y,'b');

figure(2);
subplot(2,1,1)
plot(time, Vx, 'r'); grid on; hold on
plot(time, Vx_kalk, 'b');
subplot(2,1,2)
plot(time, Vy, 'r'); grid on; hold on
plot(time, Vy_kalk, 'b');

