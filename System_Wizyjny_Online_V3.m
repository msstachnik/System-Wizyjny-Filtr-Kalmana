% System_Wizyjny_Online_V1 - pierwsza wersja systemu wizyjnego.
% najprostsza z rysowaniem przebieg�w po�o�enia oraz pr�dko�ci.
clc
clear all
close all
%% specyfikacja testu
MaxTime = 20 ; %maksymalny czas pobierania danych
display_figure = 'ON'; % je�li 'ON' to wy�wietlany jest pogl�d aktualnej pozycji
odleglosc_wierzcholkow=20; % odleglosc w cm mi�dzy wierzcho�kami kwadratu ma�ego - u�ywane w kalibracji
D_prog = 40; % maksymalna odloeg�osc pomi�dzy dwoma pobranymi punktami - warunek odsiewu szumu pomiarowego

%% parametry filtru kalmana
% e_a_sigma -  odchylenie standardowe szumu procesu okre�la jakie mo�e
% dzia�ac przy�pieszenie na k�lke - czyli zak��ci poruszanie si� k�lki
% e_pomiar_sigma - odchylenie standardowe szumu pomiaru okresla jakich
% spodziewany si� szum�w pomiarze po�ozenia k�lki
e_a_sigma = 98 / 2;          
e_pomiar_sigma = 2;

%% inicjalizacja filtru kalmana
% filtr kalmana rozbito na 2 filtry osobny dla wsp�ednej x i osobny dla
% wsp�ednej y
Px = zeros(2,2); % macierz kowariancji szumu procesu w filtrze kalmana
Py = zeros(2,2);
R = e_pomiar_sigma^2; %mo�na przyjac za sta�e mo�na uzalse�ni� od dt
B = [0; 0];         % macierz sta�a w czasie
C = [1 0];          % macierz sta�a w czasie
Vx(1) = 0;
Vy(1) = 0;
% model matematyczny uk��du
% x(i+1) = 1 * x(i) + dt * Vx(i) + dt^2/2 * e_a(i)
% Vx(i+1)= 0 * x(i) + 1 *  Vx(i) + dt * e_a(i)
% x_pomiar(i+1) = x(i+1) + e_pomiar
% 
% Q = [e_s_max^2 * dt^2 0
%      0                e_a_max^2 * dt^2]
% R = e_pomiar_max^2 %% inicjalizacja filtru kalmana
% filtr kalmana rozbito na 2 filtry osobny dla wsp�ednej x i osobny dla
% wsp�ednej y

%% definicja po��czenia

a=imaqhwinfo('winvideo');
vid = videoinput('winvideo', 1, 'YUY2_320x240');

set(vid, 'FramesPerTrigger', Inf);
set(vid, 'ReturnedColorspace', 'rgb')
% src = getselectedsource(vid);
% src.FocusMode = 'manual';
% src.ExposureMode = 'manual';
set(gcf, 'doublebuffer', 'on')

%vid = videoinput('winvideo',1,'YUY2_320x240');%'YUY2_320x240');%,'VideoResolution', '320x240');%%%);
%vidRes = get(cam, 'Resolution'); %res1=str2num(vidRes(1:3));
%res2=str2num(vidRes(5:end));%nBands = get(cam, 'NumberOfBands'); 
%hImage = image( zeros(res2,res1, 3),'Parent',handles.axes1); 
%axes(handles.axes1)%preview(cam,hImage)
%set(vid,'VideoResolution',[320 240]);

%vidRes = get(vid, 'VideoResolution');
%nBands = get(vid, 'NumberOfBands');
%hImage = image( zeros(vidRes(2), vidRes(1), nBands) );

%========================================
% ile_klatek=3;
% vid.FrameGrabInterval = ile_klatek; % okre�lenie co kt�r� klatka b�dzie brana do oblicze�
%=======================================

%% kalibracja
dl_boku = 32;
piks_to_cm=odleglosc_wierzcholkow /abs(dl_boku); % obliczenie d�ugo�ci jednego piksela
%% inicjalizacja po��czenia
start(vid); % rozpocz�cie pobierania obraz�w

k=0; % zmienna do numeracji kilejnych po�o�e� kulki
%data_old = getsnapshot(vid);
%figure;
i=0; % znacznik numeru pobranego obrazu

%inicjalizacja czasu
time(1) = 0;
tic
%% pobieranie danych
while(time(i) <MaxTime)
    i=i+1;
    
    
    vid.FramesAcquired
    data = getsnapshot(vid);
    time(i) = toc;
     
    if i == 2 % pobieranie rozmiar�w okna
        size_film = size(obraz_badany);
        wsp_x_max = size_film(2);
        wsp_y_max = size_film(1);
    end
    
    % wyciagniecie obiektow koloru czerwonego poprzez odejmowanie od siebie
    % obrazow
    obraz_badany = imsubtract(data(:,:,1), rgb2gray(data));
    %funkcja napisana do wykrywania pozycji kulki
    [wsp_x, wsp_y] = image2ball_pos(obraz_badany);
    
    if strcmp(display_figure, 'ON')
        imshow(data)
    end

    
    hold on
    
    %% j�sli wykryto pozycje k�lki to zapisa�
    if not(isnan(wsp_x + wsp_y)) 
        
        k = k + 1; %inkrementacja licznika pozycji k�lki
        time_ball(k) = time(i);

        X(k) = wsp_x * piks_to_cm;
        Y(k) = (wsp_y_max - wsp_y) * piks_to_cm;   

        %  pr�dko�� mo�na obliczy� dopiero przy minimum dw�ch pr�bkach
     %% Warunek akceptacji punktu startowego
        
        if k == 2
            Droga = sqrt((X(2) - X(1))^2 + (Y(2) - Y(1))^2);
            if Droga == 0 || Droga > Droga_prog %brak akceptacji po�o�enia pocz�tkowego
                k = 1;
                X(1) = X(2); % przypisanie bie��cej pozycji jako pozycji referencyjnej
                Y(1) = Y(2); % przypisanie bie��cej pozycji jako pozycji referencyjnej
                time_ball(1) = time_ball(2); 
            else % w przeciwnym przypadku akceptacja punktu startowego i mozna liczy� dalej
                
            end
            
        end
        
        %  pr�dko�� mo�na obliczy� dopiero przy minimum dw�ch pr�bkach
        if k >= 2
            %% Warunek akceptacji ka�dego nast�pnego punktu
            Droga = sqrt((X(k) - X(k-1))^2 + (Y(k) - Y(k-1))^2);

            if Droga == 0 || Droga > Droga_prog %brak akceptacji nowego punktu
                X = X(1:end-1); % usuni�cie tego punktu z tablicy danych
                Y = Y(1:end-1); % usuni�cie tego punktu z tablicy danych
                time_ball = time_ball(1:end-1); % usuni�cie tego punktu z tablicy danych
                k = k -1;
            else % je�li punkt jest zaakceptowany to normalne rpzetwarzanie
                
                Vx(k) = (X(k) - X(k-1)) / (time_ball(k) - time_ball(k-1));
                Vy(k) = (Y(k) - Y(k-1)) / (time_ball(k) - time_ball(k-1));

                % pokazanie na wizualizacji w kt�rym punkcie jeste�my
                if strcmp(display_figure, 'ON')
                    hold on
                    plot(wsp_x(1),wsp_y(1),'rs')
                    drawnow
                end

                %% filtracja kalmana online - mo�e zachodzi� dopiero dla K >= 2
                if k == 2 %inicjalizacja filtru i stan�w
                    % inicjalizacja poczatkowej pr�dkosci i po�ozenia w filtrze
                    % Klamana
                    X_k(1) = X(1);
                    Vx_k(1) = Vx(2) / 2;
                    Y_k(1) = Y(1);
                    Vy_k(1) = Vy(2) / 2;


                    % zastanawia� moze czemu w inicjalizacji u�yto pierwszej pr�bki
                    % po�ozenia i drugiej pr�bki pr�dko�ci. Teoretycznie
                    % powinno si� u�y� r�wnie� pierwszje pr�bki pr�dko�ci jednak
                    % pr�dko�� w pierwszej pr�bce jest r�wna 0 i ci�zko estymowac
                    % jej warto��. Dlatego przyj�to za bardziej reprezentatywn�
                    % drug� pr�bk� pr�dko�ci kt�ra w zasadzie si� wylicza na
                    % podstawie 2 pr�bek po�o�enia. pr�dko�c ta zosta�a
                    % podzielona przez 2 aby niec z�agodzi� jej wp�yw na
                    % pomiar
                    x_state = [X_k(1); Vx_k(1)];
                    y_state = [Y_k(1); Vy_k(1)];


                end
                %% w�asciwa filtracja kalmana
                if k >= 2
                   
                    % obliczanie rzeczywistego przyrostu czasu
                    dt_real = time_ball(k) - time_ball(k-1);

                    % oblicznenie macierzy kowariancji szumu procesu
                    Q = e_a_sigma^2 * [dt_real^4/4 dt_real^3/2; dt_real^3/2 dt_real^2];

                    % obliczenie dynamicznej ze wzgl�du na przyrost czasu macierzy
                    % stanu
                    A = [1 dt_real
                        0 1];

                    % filracja kalmana
                    [x_state, Px] = KalmanMS(A,B,C,R,Q,Px,0,x_state,X(k));
                    [y_state, Py] = KalmanMS(A,B,C,R,Q,Py,0,y_state,Y(k)); 

                    % wyciagni�cie z przestrzeni stanu pozycji i pr�dko�ci
                    X_k(k) = x_state(1);
                    Vx_k(k) = x_state(2);
                    Y_k(k) = y_state(1);
                    Vy_k(k) = y_state(2);


                end
            end
        end
    end
        

end
%% wizualiazacja pozycja XY
figure
plot(X,Y); grid on; hold on;
plot(X_k,Y_k,'g')
title('X Y')
xlabel('X [cm]');
ylabel('Y [cm]');
legend('pozycja zmierzona', 'pozcyja z filtru kalmana')
  
%% wizualiazacja pr�dko�ci

figure
subplot(3,1,1)
plot(time_ball,Vx); grid on; hold on;
plot(time_ball,Vx_k,'g')
title('Pr�dko��i w osi X')
xlabel('time [s]');
ylabel('V_X [cm]');
legend('pr�dko�� zmierzona', 'pr�dko�� z filtru kalmana')

subplot(3,1,2)
plot(time_ball,Vy); grid on; hold on;
plot(time_ball,Vy_k,'g')
title('Pr�dko��i w osi Y')
xlabel('time [s]');
ylabel('V_Y [cm]');
legend('pr�dko�� zmierzona', 'pr�dko�� z filtru kalmana')

V_abs = sqrt(Vx.^2 + Vy.^2);
V_abs_k = sqrt(Vx_k.^2 + Vy_k.^2);

subplot(3,1,3)
plot(time_ball,V_abs); grid on; hold on;
plot(time_ball,V_abs_k,'g')
title('Pr�dko��i absolutna')
xlabel('time [s]');
ylabel('V [cm]');
legend('pr�dko�� zmierzona', 'pr�dko�� z filtru kalmana')