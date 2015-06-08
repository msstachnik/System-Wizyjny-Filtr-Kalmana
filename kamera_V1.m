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

a=imaqhwinfo('winvideo');


%camera_name = char(a.AdaptorName);
%%camera_info = imaqhwinfo(camera_name);
%cc=camera_info.DeviceInfo;
%camera_id = char(cc.DeviceID);
%resolution = char(cc.SupportedFormats);

vid = videoinput('winvideo', 1, 'YUY2_320x240');
%vid = videoinput('winvideo', 1, 'YUY2_320x240');

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
ile_klatek=3;
vid.FrameGrabInterval = ile_klatek; % okre�lenie co kt�r� klatka b�dzie brana do oblicze�
%=======================================
start(vid); % rozpocz�cie pobierania obraz�w


% Display the video data in your GUI.
%preview(handles.vid, hImage);

%==========================================================================
% lokalizacja obiektu
%==========================================================================

klatki=30; % sta�a filmu ile klatek na sekund�
odleglosc_wierzcholkow=20; % odleglosc w cm mi�dzy wierzcho�kami kwadratu ma�ego
    %dlugosc_boku_kwadratu=35; % odleg�o�c w pix zmierzona na filmie 320x240
stala=35/odleglosc_wierzcholkow; % obliczenie sta�ej ile pikseli przypada na 1 cm
droga=0; % ustalenie warto�ci poczatkowej przebytej drogi
znak=1; % znacznik dla oblicze�
k=0; % zmienna do numeracji kilejnych po�o�e� kulki
%data_old = getsnapshot(vid);
%figure;
i=1; % znacznik numeru pobranego obrazu

%% measurement
% inicjalizacja czasu
tic 
time(i) = 0; 
while(time(i) < MaxTime) %predefiniowany czas pomiaru
    time(i) = toc;
    vid.FramesAcquired
    data = getsnapshot(vid);
     
    % wyciagniecie obiektow koloru czerwonego poprzez odejmowanie od siebie
    % obrazow
    diff_im = imsubtract(data(:,:,1), rgb2gray(data));
    %filtracja medianowa w celu usuni�cia ewentualnych zak��ce� 
    diff_im = medfilt2(diff_im, [3 3]);
    % postac binarna z progiem 0.18
    diff_im = im2bw(diff_im,0.18);
    
    % domkniecie obiekty mniejsze niz 10pix
    diff_im = bwareaopen(diff_im,10);
    
    % etykietowanie obiektow
    bw = bwlabel(diff_im, 8);
    
    
    % parametry polozenie 
    stats = regionprops(bw, 'BoundingBox', 'Centroid');
    
    
    imshow(data)
    
    hold on
    
    %oznaczanie
    if length(stats)>0
        for object = 1:length(1)
            bb = stats(object).BoundingBox;
            bc(i,:) = stats(object).Centroid;
            rectangle('Position',bb,'EdgeColor','g','LineWidth',2)
            plot(bc(i,1),bc(i,2), '-m+')
            a=text(bc(i,1)+15,bc(i,2), strcat('X: ', num2str(round(bc(i,1))), '    Y: ', num2str(round(bc(i,2)))));
            set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'red');
        end
        hold off
        i=i+1;
    end        
% Czy trzeba to oblicza� ju� teraz?        
%         if i>2
%             s=sqrt(abs((bc(i-1,1)-bc(i,1))^2+(bc(i-1,2)-bc(i,2))^2)); % pokonanna drog� przez wykryty obiekt
% 
%            % if s>0 && s<odleglosc_wierzcholkow/5 % jesli droga pokonana przez kulk� jest wieksza ni� 1/4 d� bokuto traktujemy jako zak��cenie i nie uwzgl�dniamy
%                 %plot(wsp_x(1),wsp_y(1),'rs')
%                 %k=k+1;
% 
%                 s=s/stala;
%                 droga=droga+s;
%                 V(i)=s/(ile_klatek/klatki);
%                 %wsp zapamietane % przepisanie znalezionych wsp�rzendych do zapamietania (w celu w nast�pnym kroku uzycia ich jako poprzednie po�o�enie) 
%                 %wsp_x_poprzednie=wsp_x;
%                 %wsp_y_poprzednie=wsp_y;
% 
%                 %druk(k,:)=[i,wsp_x(1),wsp_y(1),round(V),round(droga),round(s)]; % przygotowanie danych do wy�wietleniai stworzenie macierzy z danymi
%                 %disp(druk(k,:)) % wy�wietlenie berzacych wynik�w
%            % end
% 
%         %end
%         k=k+1;
%         end

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

%% filtracja klamana

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
    
    
end
%% Wizualizacja


% wyliczenie warto�ci �redniej
liczba_klatek=k;
V_sr=droga/(liczba_klatek*ile_klatek/klatki); 

disp('==================================')
disp('Warto�c �rednia')
%tekst=['Czas ruchu=',num2str(k*ile_klatek/klatki),' sekund'];
%disp(tekst)
tekst=['�rednia predko��= ',num2str(V_sr),' cm/sek'];
disp(tekst)
tekst=['Przebyty dystans= ',num2str(droga), ' cm'];
disp(tekst)


stop(vid);

% usuniecie z pamieci pobranych klatek.
flushdata(vid);


figure

plot(bc(:,1),bc(:,2),'r')
title('X-Y')
figure
plot(1:i-1,V,'b')%(1:length(druk(:,4)),druk(:,4),1:length(druk(:,5)),druk(:,5)) % wyrysowanie danych na podw�jnej osi OY
title('pr�dko�c')
%legend('X-Y','V')
%hold on
%line([1 k],[V_sr V_sr],'color','r') % narysowanie warto�ci sredniejpredko�ci
%legend('Pr�dko�� chwilowa','Pr�dko�� �rednia','Droga','Location','Best')