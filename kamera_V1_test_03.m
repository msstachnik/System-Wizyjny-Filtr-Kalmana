clc
clear all
close all
%% specyfikacja testu
%maksymalny czas pobierania danych
MaxTime = 10 ; 

% parametry filtru kalmana
% b³edy modelu
e_a_max = 1; % maksymalne przyœpieszenie kulki [cm / s^2]
e_s_max = 1; % maksymalna zmiana po³o¿enia zwiazana z zewnêtrznym wymuszeniem ruchu kulki [cm/s]

e_pomiar_max = 1; %maksymalny b³¹d pomiaru pi³eczki

%% preparation

i=2; % znacznik numeru pobranego obrazu

%% measurement (test)
% inicjalizacja czasu
time(1) = 0;
dt = 0.1;
bc(1,1) = 0; 
bc(1,2) = 0; 
while(time(i - 1) < MaxTime) %predefiniowany czas pomiaru
    dt = randp(1,0.1,0.2);
    time(i) = time(i-1) + dt;
    
    e = randp(1 , -e_a_max, e_a_max) * dt^2 + randp(1, -e_s_max, e_s_max) * dt;
    bc(i,1) =  sin(time(i) / 2) + e;
    e = randp(1 , -e_a_max, e_a_max) * dt^2 + randp(1, -e_s_max, e_s_max) * dt;
    bc(i,2) =  sin(time(i) / 1) + e;

    
    
    
    i = i + 1;
end
%% obróbka danych
%pomiar
x_pomiar = bc(:,1);
y_pomiar = bc(:,2);
%pozycja pocz¹tkowa
X(1) = x_pomiar(1);
Y(1) = y_pomiar(1);
%prêdkosc pocz¹tkowa
Vx(1) = 0;
Vy(1) = 0;
% inicjalizacja macierzy stanu
x_state = [X(1); Vx(1)];
y_state = [Y(1); Vy(1)];
% inicjalizacja macierzy kowariancji
Px = zeros(2,2);
Py = zeros(2,2);
Q = zeros(2,2);
R = e_pomiar_max^2; %mo¿na przyjac za sta³e mo¿na uzale¿niæ od dt
B = [0; 0];         % macierz sta³a w czasie
C = [1 0];          % macierz sta³a w czasie
% model matematyczny uk³¹du
% x(i+1) = 1 * x(i) + dt * Vx(i) + dt * e_s(i)
% Vx(i+1)= 0 * x(i) + 1 *  Vx(i) + dt * e_a(i)
% x_pomiar(i+1) = x(i+1) + e_pomiar
% 
% Q = [e_s_max^2 * dt^2 0
%      0                e_a_max^2 * dt^2]
% R = e_pomiar_max^2 

% filtracja klamana

for i = 2 : length(time)

    dt = time(i) - time(i-1); % model nie musi byæ z sta³ymi przyrostami czasu
    Q = [e_s_max^2 * dt^2 0
        0 e_a_max^2 * dt^2];
    A = [1 dt
        0 1];
    % mo¿na rozpatrywaæ fitlracje x i y osobno przyjmuj¹c ten sam model
    [x_state, Px, K] = KalmanMS(A,B,C,R,Q,Px,0,x_state,x_pomiar(i));
    [y_state, Py] = KalmanMS(A,B,C,R,Q,Py,0,y_state,y_pomiar(i));
    
    X(i) = x_state(1);
    Vx(i) = x_state(2);
    Y(i) = y_state(1);
    Vy(i) = y_state(2);
    
    Vx_kalk(i) = (X(i) - X(i - 1)) / dt;
    Vy_kalk(i) = (Y(i) - Y(i - 1)) / dt;
    
    Vx_pomiar(i) = (x_pomiar(i) - x_pomiar(i - 1)) / dt;
    Vy_pomiar(i) = (y_pomiar(i) - y_pomiar(i - 1)) / dt;

    
end

figure(1);
plot(x_pomiar, y_pomiar, 'r'); grid on; hold on
plot(X,Y,'b');
title('X Y');
legend('pomiar', 'filtr Kalmana')

figure(2);
subplot(2,1,1)
plot(time, Vx_pomiar, 'r'); grid on; hold on
plot(time, Vx, 'b');
plot(time, Vx_kalk, 'g'); 
legend('pomiar', 'filtr Kalmana', 'obliczenia z filtru Kalmana')
title('Vx')
subplot(2,1,2)
plot(time, Vy_pomiar, 'r'); grid on; hold on
plot(time, Vy, 'b');
plot(time, Vy_kalk, 'g'); 
legend('pomiar', 'filtr Kalmana', 'obliczenia z filtru Kalmana')
title('Vy')
