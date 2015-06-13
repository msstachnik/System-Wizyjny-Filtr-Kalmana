% System_Wizyjny_Online_V1 - pierwsza wersja systemu wizyjnego.
% najprostsza z rysowaniem przebiegów po³o¿enia oraz prêdkoœci.
clc
clear all
close all
%% specyfikacja testu
MaxTime = 20 ; %maksymalny czas pobierania danych
display_figure = 'ON'; % jeœli 'ON' to wyœwietlany jest pogl¹d aktualnej pozycji


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
tic
%% pobieranie danych
while(vid.FramesAcquired<=100)
    i=i+1;
    time(i) = toc;
    
    vid.FramesAcquired
    data = getsnapshot(vid);
    
     
    % wyciagniecie obiektow koloru czerwonego poprzez odejmowanie od siebie
    % obrazow
    obraz_badany = imsubtract(data(:,:,1), rgb2gray(data));
    %funkcja napisana do wykrywania pozycji kulki
    [raw_x, raw_y] = image2ball_pos(obraz_badany);
    
    if strcmp(display_figure, 'ON')
        imshow(data)
    end

    
    hold on
    
    %oznaczanie
    if not(isnan(wsp_x + wsp_y)) 
        
        k = k + 1; %inkrementacja licznika pozycji kólki
        time(k) = time(i);

        X(k) = wsp_x * piks_to_cm;
        Y(k) = (wsp_y_max - wsp_y) * piks_to_cm;   

        %  prêdkoœæ mo¿na obliczyæ dopiero przy minimum dwóch próbkach
        if k >= 2
            Vx(k) = (X(k) - X(k-1)) / (time(k) - time(k-1));
            Vy(k) = (Y(k) - Y(k-1)) / (time(k) - time(k-1));
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
plot(time,Vx)
title('Prêdkoœæi w osi X')
xlabel('time [s]');
ylabel('V_X [cm]');


subplot(3,1,2)
plot(time,Vy)
title('Prêdkoœæi w osi Y')
xlabel('time [s]');
ylabel('V_Y [cm]');

V_abs = sqrt(Vx.^2 + Vy.^2);


subplot(3,1,3)
plot(time,V_abs)
title('Prêdkoœæi absolutna')
xlabel('time [s]');
ylabel('V [cm]');