% System_Wizyjny_Offline_V5 - wersja z filtrem kalmana i odpornym
% pobieraniem danych oraz optymalizacj¹ okna na podstawie poziomu sigma
% procesu
% zmiana okreœlania orientacji pi³eczki + napisanie do tego zewnêtrznej
% funkcji


clc % 
close all % zamkniêcie wszystkich okien figure
clear all
%% okreœlnie parametrów
klatki=24; % sta³a filmu ile klatek na sekundê
dt = 1 / klatki; % okreslenie podstawy czasu
odleglosc_wierzcholkow=20; % odleglosc w cm miêdzy wierzcho³kami kwadratu ma³ego - u¿ywane w kalibracji
kalibracja_on_off = 'OFF'; % mo¿na napisaæ 'ON' wtedy wyœwietli siê okienko z kalibracj¹
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

% inicjalizacja macierzy kowariancji
Px = zeros(2,2); % macierz kowariancji szumu procesu w filtrze kalmana
Py = zeros(2,2);
R = e_pomiar_sigma^2; %mo¿na przyjac za sta³e mo¿na uzalse¿niæ od dt
B = [0; 0];         % macierz sta³a w czasie
C = [1 0];          % macierz sta³a w czasie
% model matematyczny uk³¹du
% x(i+1) = 1 * x(i) + dt * Vx(i) + dt^2/2 * e_a(i)
% Vx(i+1)= 0 * x(i) + 1 *  Vx(i) + dt * e_a(i)
% x_pomiar(i+1) = x(i+1) + e_pomiar
% 
% Q = [e_s_max^2 * dt^2 0
%      0                e_a_max^2 * dt^2]
% R = e_pomiar_max^2 

%% pobranie obrazu
[filename, pathname] = uigetfile('*.mp4','Wybierz plik z filmem'); %pobranie œæiezki i nazwy pliku
pathname_filename = [pathname, filename];   % sklejenie nazwy isciezki

xyloObj = VideoReader(pathname_filename);
nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;
film(1:nFrames) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'),'colormap',[]);

for i = 1 : nFrames
    film(i).cdata = read(xyloObj,i);
end
%% pobieranie rozmiarów okna
size_film = size(film(1).cdata);
wsp_x_max = size_film(2);
wsp_y_max = size_film(1);



%% kalibracja
if strcmp(kalibracja_on_off, 'ON')
    figure % otwarcie okna
    image(film(1).cdata) % wyœwietlenie obrazu na którym zaznaczamu mysza 2 punkty na jednym poziomie znajduj¹ce siê w dwóch rogach ma³ego kwadratu
    % dziêki temu bedziemy wiedziæ ile pikseli jest pomiêdzy wierzcho³kami
    % kwadratu co odpowiada 20 cm odlegloœci w rzeczywistoœci
    [x_calibration,y_calibration] = ginput(2); % pobranie wspó³rzednych dwóch punktów
    dl_boku=abs(x_calibration(2)-x_calibration(1));
    close(gcf) % zamkniêcie otwartego okna
else
    dl_boku = 32;
end
piks_to_cm=odleglosc_wierzcholkow /abs(dl_boku); % obliczenie d³ugoœci jednego piksela


%% Obróbka obrazu
k=0; % zmienna do numeracji kilejnych po³o¿eñ kulki

Vx(1) = 0;
Vy(1) = 0;

for i=2:nFrames % 2:nFrames % rozmiar 4 to iloœæ klatek, rozmiar ma 4 elemnty [a,b,c,d], a i b wys i szer, c- iloœæ wartw (3 dla RGB), d - iloœc klatek
    % zaczynami do 2 bo bierzemy klatke wczesniejsz¹ (2-1) i bierz¹c¹ 2
    %% optymalizacja okan pobierajacego dane
    
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
    
    
    
    %% pobieranie obrazu w zoptymalizowanym oknie + znajodwanie pozcyji kulki


    obraz_badany = film(i).cdata(wsp_ymin:wsp_ymax,wsp_xmin:wsp_xmax,:);
    obraz_referencyjny = film(1).cdata(wsp_ymin:wsp_ymax,wsp_xmin:wsp_xmax,:);
    obraz_roznica=medfilt2(rgb2gray(obraz_referencyjny))-medfilt2(rgb2gray(obraz_badany));
    imshow(obraz_badany)
    %funkcja napisana do wykrywania pozycji kulki
    [raw_x, raw_y] = image2ball_pos(obraz_roznica);
    
    % konwersja na pozycje absolutn¹ w pikselach - nale¿y dodaæ po³o¿enie
    % pocz¹tku zoptymalizowanego obrazu
    wsp_x=wsp_xmin + raw_x;
    wsp_y=wsp_ymin + raw_y;

    %% sprawdzenie czy pozycja jest 
    if not(isnan(wsp_x + wsp_y)) 
        k = k + 1; %inkrementacja licznika pozycji kólki
        time(k) = dt * i; %okreslenie czasu w którym dokonany jest pomiar
        
        % okreslenie po³o¿enia w zale¿noœci od piksela
        X(k) = wsp_x * piks_to_cm;
        Y(k) = (wsp_y_max - wsp_y) * piks_to_cm;
        %% Warunek akceptacji punktu startowego
        
        if k == 2
            Droga = sqrt((X(2) - X(1))^2 + (Y(2) - Y(1))^2);
            if Droga == 0 || Droga > D_prog  %brak akceptacji po³o¿enia pocz¹tkowego
                k = 1;
                X(1) = X(2); % przypisanie bie¿¹cej pozycji jako pozycji referencyjnej
                Y(1) = Y(2); % przypisanie bie¿¹cej pozycji jako pozycji referencyjnej
                time(1) = time(2); 
            else % w przeciwnym przypadku akceptacja punktu startowego i mozna liczyæ dalej
                
            end
            
        end
        
        %  prêdkoœæ mo¿na obliczyæ dopiero przy minimum dwóch próbkach
        if k >= 2
            %% Warunek akceptacji ka¿dego nastêpnego punktu
            Droga = sqrt((X(k) - X(k-1))^2 + (Y(k) - Y(k-1))^2);

            if Droga == 0 || Droga > D_prog %brak akceptacji nowego punktu
                X = X(1:end-1); % usuniêcie tego punktu z tablicy danych
                Y = Y(1:end-1); % usuniêcie tego punktu z tablicy danych
                time = time(1:end-1); % usuniêcie tego punktu z tablicy danych
                k = k -1;
            else % jeœli punkt jest zaakceptowany to normalne rpzetwarzanie
                
                Vx(k) = (X(k) - X(k-1)) / (time(k) - time(k-1));
                Vy(k) = (Y(k) - Y(k-1)) / (time(k) - time(k-1));

                % pokazanie na wizualizacji w którym punkcie jesteœmy
                hold on
                plot(wsp_x - wsp_xmin,wsp_y - wsp_ymin,'rs')
                drawnow

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


                end
                %% w³asciwa filtracja kalmana
                if k >= 2
                   
                    % obliczanie rzeczywistego przyrostu czasu
                    dt_real = time(k) - time(k-1);

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
                    
                    %% 


                end
            end
        end

            
        

    end

    hold off
    
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
plot(time,Vx); grid on; hold on;
plot(time,Vx_k,'g')
title('Prêdkoœæi w osi X')
xlabel('time [s]');
ylabel('V_X [cm]');
legend('prêdkoœæ zmierzona', 'prêdkoœæ z filtru kalmana')

subplot(3,1,2)
plot(time,Vy); grid on; hold on;
plot(time,Vy_k,'g')
title('Prêdkoœæi w osi Y')
xlabel('time [s]');
ylabel('V_Y [cm]');
legend('prêdkoœæ zmierzona', 'prêdkoœæ z filtru kalmana')

V_abs = sqrt(Vx.^2 + Vy.^2);
V_abs_k = sqrt(Vx_k.^2 + Vy_k.^2);

subplot(3,1,3)
plot(time,V_abs); grid on; hold on;
plot(time,V_abs_k,'g')
title('Prêdkoœæi absolutna')
xlabel('time [s]');
ylabel('V [cm]');
legend('prêdkoœæ zmierzona', 'prêdkoœæ z filtru kalmana')



