% System_Wizyjny_Online_V1 - pierwsza wersja systemu wizyjnego.
% najprostsza z rysowaniem przebiegów po³o¿enia oraz prêdkoœci.
clc
clear all
close all
%% specyfikacja testu
MaxTime = 20 ; %maksymalny czas pobierania danych
display_figure = 'ON'; % jeœli 'ON' to wyœwietlany jest pogl¹d aktualnej pozycji
odleglosc_wierzcholkow=20; % odleglosc w cm miêdzy wierzcho³kami kwadratu ma³ego - u¿ywane w kalibracji


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
    i=i+1;
    
    
    vid.FramesAcquired
    data = getsnapshot(vid);
    time(i) = toc;
    
    if i == 2 % pobieranie rozmiarów okna - wsytarczy tylko raz
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
    
    %% jêsli wykryto pozycje kólki to zapisaæ
    if not(isnan(wsp_x + wsp_y)) 
        
        k = k + 1; %inkrementacja licznika pozycji kólki
        time_ball(k) = time(i);

        X(k) = wsp_x * piks_to_cm;
        Y(k) = (wsp_y_max - wsp_y) * piks_to_cm;   

        %  prêdkoœæ mo¿na obliczyæ dopiero przy minimum dwóch próbkach
        if k >= 2
            Vx(k) = (X(k) - X(k-1)) / (time_ball(k) - time_ball(k-1));
            Vy(k) = (Y(k) - Y(k-1)) / (time_ball(k) - time_ball(k-1));
        end
        
        
        if strcmp(display_figure, 'ON')
            % pokazanie na wizualizacji w którym punkcie jesteœmy
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
  
%% wizualiazacja prêdkoœci

figure
subplot(3,1,1)
plot(time_ball,Vx)
title('Prêdkoœæi w osi X')
xlabel('time [s]');
ylabel('V_X [cm]');


subplot(3,1,2)
plot(time_ball,Vy)
title('Prêdkoœæi w osi Y')
xlabel('time [s]');
ylabel('V_Y [cm]');

V_abs = sqrt(Vx.^2 + Vy.^2);


subplot(3,1,3)
plot(time_ball,V_abs)
title('Prêdkoœæi absolutna')
xlabel('time [s]');
ylabel('V [cm]');