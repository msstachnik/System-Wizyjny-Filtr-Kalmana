% System_Wizyjny_Offline_V2 - wersja z filtrem kalmana


clc % 
close all % zamkni�cie wszystkich okien figure
clear all
%% okre�lnie parametr�w
klatki=24; % sta�a filmu ile klatek na sekund�
dt = 1 / klatki; % okreslenie podstawy czasu
odleglosc_wierzcholkow=20; % odleglosc w cm mi�dzy wierzcho�kami kwadratu ma�ego - u�ywane w kalibracji
kalibracja_on_off = 'OFF'; % mo�na napisa� 'ON' wtedy wy�wietli si� okienko z kalibracj�


%% parametry filtru kalmana
% e_a_sigma -  odchylenie standardowe szumu procesu okre�la jakie mo�e
% dzia�ac przy�pieszenie na k�lke - czyli zak��ci poruszanie si� k�lki
% e_pomiar_sigma - odchylenie standardowe szumu pomiaru okresla jakich
% spodziewany si� szum�w pomiarze po�ozenia k�lki
e_a_sigma = 98 / 2 ; % [cm/s^2]          
e_pomiar_sigma = 2; % [cm]

%% inicjalizacja filtru kalmana
% filtr kalmana rozbito na 2 filtry osobny dla wsp�ednej x i osobny dla
% wsp�ednej y

% inicjalizacja macierzy kowariancji
Px = zeros(2,2); % macierz kowariancji szumu procesu w filtrze kalmana
Py = zeros(2,2);
R = e_pomiar_sigma^2; %mo�na przyjac za sta�e mo�na uzalse�ni� od dt
B = [0; 0];         % macierz sta�a w czasie
C = [1 0];          % macierz sta�a w czasie
% model matematyczny uk��du
% x(i+1) = 1 * x(i) + dt * Vx(i) + dt^2/2 * e_a(i)
% Vx(i+1)= 0 * x(i) + 1 *  Vx(i) + dt * e_a(i)
% x_pomiar(i+1) = x(i+1) + e_pomiar
% 
% Q = [e_s_max^2 * dt^2 0
%      0                e_a_max^2 * dt^2]
% R = e_pomiar_max^2 

%% pobranie obrazu
[filename, pathname] = uigetfile('*.mp4','Wybierz plik z filmem'); %pobranie ��iezki i nazwy pliku
pathname_filename = [pathname, filename];   % sklejenie nazwy isciezki

xyloObj = VideoReader(pathname_filename);
nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;
film(1:nFrames) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'),'colormap',[]);

for i = 1 : nFrames
    film(i).cdata = read(xyloObj,i);
end
%% pobieranie rozmiar�w okna
size_film = size(film(1).cdata);
wsp_x_max = size_film(2);
wsp_y_max = size_film(1);



%% kalibracja
if strcmp(kalibracja_on_off, 'ON')
    figure % otwarcie okna
    image(film(1).cdata) % wy�wietlenie obrazu na kt�rym zaznaczamu mysza 2 punkty na jednym poziomie znajduj�ce si� w dw�ch rogach ma�ego kwadratu
    % dzi�ki temu bedziemy wiedzi� ile pikseli jest pomi�dzy wierzcho�kami
    % kwadratu co odpowiada 20 cm odleglo�ci w rzeczywisto�ci
    [x_calibration,y_calibration] = ginput(2); % pobranie wsp�rzednych dw�ch punkt�w
    dl_boku=abs(x_calibration(2)-x_calibration(1));
    close(gcf) % zamkni�cie otwartego okna
else
    dl_boku = 32;
end
piks_to_cm=odleglosc_wierzcholkow /abs(dl_boku); % obliczenie d�ugo�ci jednego piksela


%% Obr�bka obrazu
k=0; % zmienna do numeracji kilejnych po�o�e� kulki

Vx(1) = 0;
Vy(1) = 0;

  for i=5:nFrames % 2:nFrames % rozmiar 4 to ilo�� klatek, rozmiar ma 4 elemnty [a,b,c,d], a i b wys i szer, c- ilo�� wartw (3 dla RGB), d - ilo�c klatek
    % zaczynami do 2 bo bierzemy klatke wczesniejsz� (2-1) i bierz�c� 2

    %% pobieranie obrazu + konwertowanie go na macierz zancznik�w logicznych
    % kt�re nastepnie s� u�ywane w warunku sprawdzajacycm wykrycie obrazu
    
    obraz_roznica=medfilt2(rgb2gray(film(i-4).cdata))-medfilt2(rgb2gray(film(i).cdata));
    imshow(rgb2gray(film(i-2).cdata))
    wynik = obraz_roznica(5:end-5,5:end-5);%rgb2gray(obraz_roznica); % przepisanie obrazu r�nicowego czarno bia�ego do zmiennej wynik
    %%%%%%%%%%%%%%%%%%
    znacznik_logiczny = wynik>50;
    %%imshow(aa)
    
    
    
    %%%%%%%%%%%%%%%%%%
    % szukanie wsp�rzednych kulki, kt�ra w obrazie bw jest bia�a - szukamy
    % miejsca gdzie jest najwi�cej pikseli bia�ych kt�re nale�a do kulki
    %poziom=sum(aa); %sum((wynik)) sumowanie ilo�ci pikseli (bo piksele maj� warto�ci 0 i 1) w poziomie
    %wsp_x=find(poziom==max(poziom));% szukanie miejsca indeksu gdzie w wektorze w kt�rym sa sumy dla ka�ej kolumny wynik wyst�puj� maximum
    % wsp_x=mean( wsp_x); opcja gdyby by�y wskazane 2 lub 3 maksima to
    % wybierzemy u�rednion� wsp�rzedn�

    %pion=sum(aa');%sum(wynik') sumowanie ilo�ci pikseli (bo piksele maj� warto�ci 0 i 1) w pionie (naprawde w poziomie ale macierz jest obr�cona - transponowana  co daje ten sam efekt)
    %wsp_y=find(pion==max(pion));  % szukanie miejsca indeksu gdzie w wektorze w kt�rym sa sumy dla ka�ej kolumny wynik wyst�puj� maximum
    % wsp_y=mean( wsp_y); opcja gdyby by�y wskazane 2 lub 3 maksima to
    % wybierzemy u�rednion� wsp�rzedn�
    %% wyliczanie pozycji pi�eczni + g��wny warunek determinujacy czy wykryto czy te� nie pi�eczke
    if sum(sum(znacznik_logiczny))<100 && sum(sum(znacznik_logiczny))>5  %max(max(wynik))>50 &&sum(poziom)>500 && sum(pion)>500 % za�o�ono �e powierznia kulki ma mniej ni� 700 pix (w rzeczywisto�ci dla rozdzielczo�ci 800x600 ma ok.440 pix , je�eli na obrazie suma pikseli bia�ych jest mniejsza ni� 700
        k = k + 1; %inkrementacja licznika pozycji k�lki
        time(k) = dt * i; %okreslenie czasu w kt�rym dokonany jest pomiar
        
        % szukanie wsp�rzednych kulki, kt�ra w obrazie bw jest bia�a - szukamy
        % miejsca gdzie jest najwi�cej pikseli bia�ych kt�re nale�a do kulki
        poziom=sum((znacznik_logiczny)); %sum((wynik)) sumowanie ilo�ci pikseli (bo piksele maj� warto�ci 0 i 1) w poziomie
        wsp_x=find(poziom==max(poziom));% szukanie miejsca indeksu gdzie w wektorze w kt�rym sa sumy dla ka�ej kolumny wynik wyst�puj� maximum
        wsp_x=mean( wsp_x);
        
        pion=sum(znacznik_logiczny,2);%sum(wynik') sumowanie ilo�ci pikseli (bo piksele maj� warto�ci 0 i 1) w pionie (naprawde w poziomie ale macierz jest obr�cona - transponowana  co daje ten sam efekt)
        wsp_y=find(pion==max(pion));  % szukanie miejsca indeksu gdzie w wektorze w kt�rym sa sumy dla ka�ej kolumny wynik wyst�puj� maximum
        wsp_y=mean( wsp_y);
        
        % okreslenie po�o�enia w zale�no�ci od piksela
        X(k) = wsp_x * piks_to_cm;
        Y(k) = (wsp_y_max - wsp_y) * piks_to_cm;
        
        %  pr�dko�� mo�na obliczy� dopiero przy minimum dw�ch pr�bkach
        if k >= 2
            Vx(k) = (X(k) - X(k-1)) / (time(k) - time(k-1));
            Vy(k) = (Y(k) - Y(k-1)) / (time(k) - time(k-1));
        end


        
        % pokazanie na wizualizacji w kt�rym punkcie jeste�my
        hold on
        plot(wsp_x(1),wsp_y(1),'rs')
        drawnow

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
            % podstawie 2 pr�bek po�o�enia
            x_state = [X_k(1); Vx_k(1)];
            y_state = [Y_k(1); Vy_k(1)];


        end

        if k >= 2
            % w�asciwa filtracja kalmana

            % obliczanie rzeczywistego przyrostu czasu
            dt_real = time(k) - time(k-1);

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
  
%% wizualiazacja pr�dko�ci

figure
subplot(3,1,1)
plot(time,Vx); grid on; hold on;
plot(time,Vx_k,'g')
title('Pr�dko��i w osi X')
xlabel('time [s]');
ylabel('V_X [cm]');
legend('pr�dko�� zmierzona', 'pr�dko�� z filtru kalmana')

subplot(3,1,2)
plot(time,Vy); grid on; hold on;
plot(time,Vy_k,'g')
title('Pr�dko��i w osi Y')
xlabel('time [s]');
ylabel('V_Y [cm]');
legend('pr�dko�� zmierzona', 'pr�dko�� z filtru kalmana')

V_abs = sqrt(Vx.^2 + Vy.^2);
V_abs_k = sqrt(Vx_k.^2 + Vy_k.^2);

subplot(3,1,3)
plot(time,V_abs); grid on; hold on;
plot(time,V_abs_k,'g')
title('Pr�dko��i absolutna')
xlabel('time [s]');
ylabel('V [cm]');
legend('pr�dko�� zmierzona', 'pr�dko�� z filtru kalmana')


%   hist(time(7:end),20)
%    median(time)
% 
% % wyliczenie warto�ci �redniej
% ile_klatek=k;
% V_sr=droga/(ile_klatek*1/klatki); 
% 
% disp('==================================')
% disp('Warto�c �rednia')
% tekst=['Czas ruchu=',num2str(ile_klatek*1/klatki),' sekund'];
% disp(tekst)
% tekst=['�rednia predko��= ',num2str(V_sr),' cm/sek'];
% disp(tekst)
% tekst=['Przebyty dustans= ',num2str(droga), ' cm'];
% disp(tekst)
% close gcf
% 
% 
% figure
% 
% plotyy(1:length(druk(:,4)),druk(:,4),1:length(druk(:,5)),druk(:,5)) % wyrysowanie danych na podw�jnej osi OY
% 
% hold on
% line([1 k],[V_sr V_sr],'color','r') % narysowanie warto�ci sredniejpredko�ci
% legend('Pr�dko�� chwilowa','Pr�dko�� �rednia','Droga','Location','Best')
% 
% 
% %++++++++++++++++++++++++++++++++++++++++++++++++++++
% % predko��
% 
% t=(druk(:,1)-druk(1,1))*1/( xyloObj.FrameRate);
% v=medfilt1(druk(:,4),5); % filtracja medianowa pr�dko�ci dla 5 warto�ci
% 
% 
% %figure
% f = fit(t(2:end),v(2:end),'poly1');
% %plot(f,t(2:end),v(2:end));
% a=f.p1;
% b=f.p2;
% wynik=roots([a,b]);
% 
% 
% x=0:abs(wynik);
% y=a*x+b;
% 
% figure
% plot(t(2:end),v(2:end),'.r');
% hold on
% plot(x,y)
% title('Estymacja pr�dko�ci')
% 
% % =========================================================
% % po�ozenie
% xx=druk(:,2);%2);
% yy=druk(:,3);%3);3
% 
% figure
% f = fit(xx,yy,'poly1');
% plot(f,xx,yy);
% a=f.p1;
% b=f.p2;
% %wynik=roots([a,b]);
% title('Wykres po�ozenia')
% 
% 
% x=0:1000;
% y=a*x+b;
% 
% figure
% plot(xx,yy,'.r')
% hold on
% plot(x,y)
% title('Estymacja po�ozenia')
% 

