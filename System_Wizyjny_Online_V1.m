% System_Wizyjny_Online_V1 - pierwsza wersja systemu wizyjnego.
% najprostsza z rysowaniem przebieg�w po�o�enia oraz pr�dko�ci.
clc
clear all
close all
%% specyfikacja testu
MaxTime = 20 ; %maksymalny czas pobierania danych
display_figure = 'ON'; % je�li 'ON' to wy�wietlany jest pogl�d aktualnej pozycji
odleglosc_wierzcholkow=20; % odleglosc w cm mi�dzy wierzcho�kami kwadratu ma�ego - u�ywane w kalibracji


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
     
    % wyciagniecie obiektow koloru czerwonego poprzez odejmowanie od siebie
    % obrazow
    obraz_badany = imsubtract(data(:,:,1), rgb2gray(data));
    %funkcja napisana do wykrywania pozycji kulki
    [wsp_x, wsp_y] = image2ball_pos(obraz_badany);
    
    if strcmp(display_figure, 'ON')
        imshow(data)
    end

    
    hold on
    
    %oznaczanie
    if not(isnan(wsp_x + wsp_y)) 
        
        k = k + 1; %inkrementacja licznika pozycji k�lki
        time_ball(k) = time(i);

        X(k) = wsp_x * piks_to_cm;
        Y(k) = (wsp_y_max - wsp_y) * piks_to_cm;   

        %  pr�dko�� mo�na obliczy� dopiero przy minimum dw�ch pr�bkach
        if k >= 2
            Vx(k) = (X(k) - X(k-1)) / (time_ball(k) - time_ball(k-1));
            Vy(k) = (Y(k) - Y(k-1)) / (time_ball(k) - time_ball(k-1));
        end
        
        
        if strcmp(display_figure, 'ON')
            % pokazanie na wizualizacji w kt�rym punkcie jeste�my
            hold on
            plot(wsp_x(1),wsp_y(1),'rs')
            drawnow
        end
    end
        

end
%% wizualiazacja pozycja XY
figure
plot(X,Y);
title('X Y')
xlabel('X [cm]');
ylabel('Y [cm]');
  
%% wizualiazacja pr�dko�ci

figure
subplot(3,1,1)
plot(time_ball,Vx)
title('Pr�dko��i w osi X')
xlabel('time [s]');
ylabel('V_X [cm]');


subplot(3,1,2)
plot(time_ball,Vy)
title('Pr�dko��i w osi Y')
xlabel('time [s]');
ylabel('V_Y [cm]');

V_abs = sqrt(Vx.^2 + Vy.^2);


subplot(3,1,3)
plot(time_ball,V_abs)
title('Pr�dko��i absolutna')
xlabel('time [s]');
ylabel('V [cm]');