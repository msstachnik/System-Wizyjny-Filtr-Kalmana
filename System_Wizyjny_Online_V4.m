% System_Wizyjny_Online_V4 - Zaimplementowany filtr Kalmana i optymalizacja
% okna pomiarowego
clc
clear all
close all
%% specyfikacja testu
MaxTime = 20 ; %maksymalny czas pobierania danych
display_figure = 'ON'; % jeœli 'ON' to wyœwietlany jest pogl¹d aktualnej pozycji
odleglosc_wierzcholkow=20; % odleglosc w cm miêdzy wierzcho³kami kwadratu ma³ego - u¿ywane w kalibracji
D_prog = 40; % maksymalna odloeg³osc pomiêdzy dwoma pobranymi punktami - warunek odsiewu szumu pomiarowego
srednica_pilki = 5; %srednica pilki w cm

% poziom sigma okreœla kompromis pomiêdzy optymalizacj¹ obliczeniow¹  a
% dok³adnoœci¹ pomiarów. Im wiêkszy poziom sigma tym wiêksze okno pomiarowe
poziom_sigma = 6;  % poziom sigma nie powininem byæ mniejszy ni¿ 1.5

%% parametry filtru kalmana
% e_a_sigma -  odchylenie standardowe szumu procesu okreœla jakie mo¿e
% dzia³ac przyœpieszenie na kólke - czyli zak³óci poruszanie siê kólki
% e_pomiar_sigma - odchylenie standardowe szumu pomiaru okresla jakich
% spodziewany siê szumów pomiarze po³ozenia kólki
e_a_sigma = 98 / 2;          
e_pomiar_sigma = 2;

%% inicjalizacja filtru kalmana
% filtr kalmana rozbito na 2 filtry osobny dla wspó³¿ednej x i osobny dla
% wspó³¿ednej y
% Px = zeros(2,2); % macierz kowariancji szumu procesu w filtrze kalmana
% Py = zeros(2,2);
R = e_pomiar_sigma^2; %mo¿na przyjac za sta³e mo¿na uzalse¿niæ od dt
B = [0; 0];         % macierz sta³a w czasie
C = [1 0];          % macierz sta³a w czasie
Vx(1) = 0;
Vy(1) = 0;
% model matematyczny uk³¹du
% x(i+1) = 1 * x(i) + dt * Vx(i) + dt^2/2 * e_a(i)
% Vx(i+1)= 0 * x(i) + 1 *  Vx(i) + dt * e_a(i)
% x_pomiar(i+1) = x(i+1) + e_pomiar
% 
% Q = [e_s_max^2 * dt^2 0
%      0                e_a_max^2 * dt^2]
% R = e_pomiar_max^2 %% inicjalizacja filtru kalmana
% filtr kalmana rozbito na 2 filtry osobny dla wspó³¿ednej x i osobny dla
% wspó³¿ednej y

%% definicja po³¹czenia

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
% vid.FrameGrabInterval = ile_klatek; % okreœlenie co któr¹ klatka bêdzie brana do obliczeñ
%=======================================

%% kalibracja
dl_boku = 32;
piks_to_cm=odleglosc_wierzcholkow /abs(dl_boku); % obliczenie d³ugoœci jednego piksela
%% inicjalizacja po³¹czenia
start(vid); % rozpoczêcie pobierania obrazów

k=0; % zmienna do numeracji kilejnych po³o¿eñ kulki
%data_old = getsnapshot(vid);
%figure;
i=0; % znacznik numeru pobranego obrazu

%inicjalizacja czasu
time(1) = 0;
tic
%% pobieranie danych
while(time(i) <MaxTime)
    %% optymalizacja okan pobierajacego dane
    
    
    
    i=i+1;
    
    
    vid.FramesAcquired
    data = getsnapshot(vid);
    time(i) = toc;
     
    if i == 2 % pobieranie rozmiarów okna
        size_film = size(data);
        wsp_x_max = size_film(2);
        wsp_y_max = size_film(1);
    end
    
    
    if k >= 2 % optymalizacja jest w³¹czona dopiero po zaakceptowaniu punktu pocz¹tkowego
        % i - obecna iteracja, (time(k) / dt) - ostatnia zaakceptowana iteracja pozycji
        dt_okna = (i - (time(k) / dt)) * dt;         
        
        % estymacja po³o¿enia X Y (nie licz¹c wariancji)
        X_est = X_k(k) + Vx_k(k) * dt_okna;
        Y_est = Y_k(k) + Vy_k(k) * dt_okna;
        
        % konwersja estymacji po³ozenia X Y po wspó³¿êdnych wsp_x wsp_y
        wsp_x_okna = round(X_est / piks_to_cm);
        wsp_y_okna = round(wsp_y_max - (Y_est / piks_to_cm));
        
        %sprawdzenia czy wspó³¿êdne nie znajduj¹ sie poza obszarem filmu
        wsp_x_okna = down_up_limit(wsp_x_okna, 1, wsp_x_max);
        wsp_y_okna = down_up_limit(wsp_y_okna, 1, wsp_y_max);

        % bok okna wynika z statystyki. odchylenia stadardowe b³edu pozycji
        % z macierzy stany wynosz¹ Px(1,1), prêdkosci Px(2,2), a b³¹d
        % pomiaru e_pomiar_sigma
        sigma_x = Px(1,1);
        sigma_Vx = Px(2,2) * dt_okna;
        sigma_y = Py(1,1);
        sigma_Vy = Py(2,2) * dt_okna;
        % rzeczywista signa u¿ywana do obliczania boku prostok¹ta na
        % podstawie poziomu sigma
        real_sigma = down_up_limit(poziom_sigma - 1.5, 0.1, inf);
        
        

        % obliczenie boków prostok¹ta - estymacja dok³¹dnosci pomiaru
        bok_okna_x = srednica_pilki + real_sigma * sqrt(sigma_x^2 + sigma_Vx^2 + e_pomiar_sigma^2);
        bok_okna_y = srednica_pilki + real_sigma * sqrt(sigma_y^2 + sigma_Vy^2 + e_pomiar_sigma^2);
        
        % konwersja na wspó³¿êdne
        wsp_bok_okna_x = round( bok_okna_x / piks_to_cm);
        wsp_bok_okna_y = round( bok_okna_y / piks_to_cm);
        
        % wyliczenie pozcyji skrajnych prostok¹ta
        wsp_xmin = down_up_limit(wsp_x_okna - wsp_bok_okna_x, 1, wsp_x_max);
        wsp_ymin = down_up_limit(wsp_y_okna - wsp_bok_okna_y, 1, wsp_y_max);
        wsp_xmax = down_up_limit(wsp_x_okna + wsp_bok_okna_x, 1, wsp_x_max);
        wsp_ymax = down_up_limit(wsp_y_okna + wsp_bok_okna_y, 1, wsp_y_max);
        
    else % jeœli nie ma akceptacji pierwszych 2 punktów to skanuj ca³a przestrzeñ
        wsp_xmin = 1;
        wsp_ymin = 1;
        wsp_xmax = wsp_x_max;
        wsp_ymax = wsp_y_max;
        
    
    end
    
    
    % wyciagniecie obiektow koloru czerwonego poprzez odejmowanie od siebie
    % obrazow
    data = data(wsp_ymin:wsp_ymax,wsp_xmin:wsp_xmax,:);
    obraz_badany = imsubtract(data(:,:,1), rgb2gray(data));
    %funkcja napisana do wykrywania pozycji kulki
    [wsp_x, wsp_y] = image2ball_pos(obraz_badany);
    
    if strcmp(display_figure, 'ON')
        imshow(data)
        hold on
    end

    
    
    
    %% jêsli wykryto pozycje kólki to zapisaæ
    if not(isnan(wsp_x + wsp_y)) 
        
        k = k + 1; %inkrementacja licznika pozycji kólki
        time_ball(k) = time(i);

        X(k) = wsp_x * piks_to_cm;
        Y(k) = (wsp_y_max - wsp_y) * piks_to_cm;   

        %  prêdkoœæ mo¿na obliczyæ dopiero przy minimum dwóch próbkach
     %% Warunek akceptacji punktu startowego
        
        if k == 2
            Droga = sqrt((X(2) - X(1))^2 + (Y(2) - Y(1))^2);
            if Droga == 0 || Droga > Droga_prog %brak akceptacji po³o¿enia pocz¹tkowego
                k = 1;
                X(1) = X(2); % przypisanie bie¿¹cej pozycji jako pozycji referencyjnej
                Y(1) = Y(2); % przypisanie bie¿¹cej pozycji jako pozycji referencyjnej
                time_ball(1) = time_ball(2); 
            else % w przeciwnym przypadku akceptacja punktu startowego i mozna liczyæ dalej
                
            end
            
        end
        
        %  prêdkoœæ mo¿na obliczyæ dopiero przy minimum dwóch próbkach
        if k >= 2
            %% Warunek akceptacji ka¿dego nastêpnego punktu
            Droga = sqrt((X(k) - X(k-1))^2 + (Y(k) - Y(k-1))^2);

            if Droga == 0 || Droga > Droga_prog %brak akceptacji nowego punktu
                X = X(1:end-1); % usuniêcie tego punktu z tablicy danych
                Y = Y(1:end-1); % usuniêcie tego punktu z tablicy danych
                time_ball = time_ball(1:end-1); % usuniêcie tego punktu z tablicy danych
                k = k -1;
            else % jeœli punkt jest zaakceptowany to normalne rpzetwarzanie
                
                Vx(k) = (X(k) - X(k-1)) / (time_ball(k) - time_ball(k-1));
                Vy(k) = (Y(k) - Y(k-1)) / (time_ball(k) - time_ball(k-1));

                % pokazanie na wizualizacji w którym punkcie jesteœmy
                if strcmp(display_figure, 'ON')
                    hold on
                    plot(wsp_x(1),wsp_y(1),'rs')
                    drawnow
                end

                %% filtracja kalmana online - mo¿e zachodziæ dopiero dla K >= 2
                if k == 2 %inicjalizacja filtru i stanów
                    % inicjalizacja poczatkowej prêdkosci i po³ozenia w filtrze
                    % Klamana
                    X_k(1) = X(1);
                    Vx_k(1) = Vx(2) / 2;
                    Y_k(1) = Y(1);
                    Vy_k(1) = Vy(2) / 2;


                    % zastanawiaæ moze czemu w inicjalizacji u¿yto pierwszej próbki
                    % po³ozenia i drugiej próbki prêdkoœci. Teoretycznie
                    % powinno siê u¿yæ równie¿ pierwszje próbki prêdkoœci jednak
                    % prêdkoœæ w pierwszej próbce jest równa 0 i ciêzko estymowac
                    % jej wartoœæ. Dlatego przyjêto za bardziej reprezentatywn¹
                    % drug¹ próbkê prêdkoœci która w zasadzie siê wylicza na
                    % podstawie 2 próbek po³o¿enia. prêdkoœc ta zosta³a
                    % podzielona przez 2 aby niec z³agodziæ jej wp³yw na
                    % pomiar
                    x_state = [X_k(1); Vx_k(1)];
                    y_state = [Y_k(1); Vy_k(1)];
                    
                    % obliczanie rzeczywistego przyrostu czasu
                    dt_real = time_ball(k) - time_ball(k-1);
                    
                    % oliczenie inicjalizacji macierzy kowariancji b³êdu
                    Px = diag([e_pomiar_sigma, e_a_sigma * dt_real]); % macierz kowariancji szumu procesu w filtrze kalmana
                    Py = Px;             

                end
                %% w³asciwa filtracja kalmana
                if k >= 2
                   
                    % obliczanie rzeczywistego przyrostu czasu
                    dt_real = time_ball(k) - time_ball(k-1);

                    % oblicznenie macierzy kowariancji szumu procesu
                    Q = e_a_sigma^2 * [dt_real^4/4 dt_real^3/2; dt_real^3/2 dt_real^2];

                    % obliczenie dynamicznej ze wzglêdu na przyrost czasu macierzy
                    % stanu
                    A = [1 dt_real
                        0 1];

                    % filracja kalmana
                    [x_state, Px] = KalmanMS(A,B,C,R,Q,Px,0,x_state,X(k));
                    [y_state, Py] = KalmanMS(A,B,C,R,Q,Py,0,y_state,Y(k)); 

                    % wyciagniêcie z przestrzeni stanu pozycji i prêdkoœci
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
  
%% wizualiazacja prêdkoœci

figure
subplot(3,1,1)
plot(time_ball,Vx); grid on; hold on;
plot(time_ball,Vx_k,'g')
title('Prêdkoœæi w osi X')
xlabel('time [s]');
ylabel('V_X [cm]');
legend('prêdkoœæ zmierzona', 'prêdkoœæ z filtru kalmana')

subplot(3,1,2)
plot(time_ball,Vy); grid on; hold on;
plot(time_ball,Vy_k,'g')
title('Prêdkoœæi w osi Y')
xlabel('time [s]');
ylabel('V_Y [cm]');
legend('prêdkoœæ zmierzona', 'prêdkoœæ z filtru kalmana')

V_abs = sqrt(Vx.^2 + Vy.^2);
V_abs_k = sqrt(Vx_k.^2 + Vy_k.^2);

subplot(3,1,3)
plot(time_ball,V_abs); grid on; hold on;
plot(time_ball,V_abs_k,'g')
title('Prêdkoœæi absolutna')
xlabel('time [s]');
ylabel('V [cm]');
legend('prêdkoœæ zmierzona', 'prêdkoœæ z filtru kalmana')