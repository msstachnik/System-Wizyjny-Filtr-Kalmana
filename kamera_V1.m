clc
clear all
close all
%% specyfikacja testu
%maksymalny czas pobierania danych
MaxTime = 20 ; 


%%

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
vid.FrameGrabInterval = ile_klatek; % okreœlenie co któr¹ klatka bêdzie brana do obliczeñ
%=======================================
start(vid); % rozpoczêcie pobierania obrazów


% Display the video data in your GUI.
%preview(handles.vid, hImage);

%==========================================================================
% lokalizacja obiektu
%==========================================================================

klatki=30; % sta³a filmu ile klatek na sekundê
odleglosc_wierzcholkow=20; % odleglosc w cm miêdzy wierzcho³kami kwadratu ma³ego
    %dlugosc_boku_kwadratu=35; % odleg³oœc w pix zmierzona na filmie 320x240
stala=35/odleglosc_wierzcholkow; % obliczenie sta³ej ile pikseli przypada na 1 cm
droga=0; % ustalenie wartoœci poczatkowej przebytej drogi
znak=1; % znacznik dla obliczeñ
k=0; % zmienna do numeracji kilejnych po³o¿eñ kulki
%data_old = getsnapshot(vid);
%figure;
i=1; % znacznik numeru pobranego obrazu

while(vid.FramesAcquired<=100)
    vid.FramesAcquired
        data = getsnapshot(vid);
     
    % wyciagniecie obiektow koloru czerwonego poprzez odejmowanie od siebie
    % obrazow
    diff_im = imsubtract(data(:,:,1), rgb2gray(data));
    %filtracja medianowa w celu usuniêcia ewentualnych zak³óceñ 
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
        
        if i>2
            s=sqrt(abs((bc(i-1,1)-bc(i,1))^2+(bc(i-1,2)-bc(i,2))^2)); % pokonanna drog¹ przez wykryty obiekt

           % if s>0 && s<odleglosc_wierzcholkow/5 % jesli droga pokonana przez kulkê jest wieksza ni¿ 1/4 d³ bokuto traktujemy jako zak³ócenie i nie uwzglêdniamy
                %plot(wsp_x(1),wsp_y(1),'rs')
                %k=k+1;

                s=s/stala;
                droga=droga+s;
                V(i)=s/(ile_klatek/klatki);
                %wsp zapamietane % przepisanie znalezionych wspó³rzendych do zapamietania (w celu w nastêpnym kroku uzycia ich jako poprzednie po³o¿enie) 
                %wsp_x_poprzednie=wsp_x;
                %wsp_y_poprzednie=wsp_y;

                %druk(k,:)=[i,wsp_x(1),wsp_y(1),round(V),round(droga),round(s)]; % przygotowanie danych do wyœwietleniai stworzenie macierzy z danymi
                %disp(druk(k,:)) % wyœwietlenie berzacych wyników
           % end

        %end
        k=k+1;
        end
        i=i+1;
    end
end

% wyliczenie wartoœci œredniej
liczba_klatek=k;
V_sr=droga/(liczba_klatek*ile_klatek/klatki); 

disp('==================================')
disp('Wartoœc œrednia')
%tekst=['Czas ruchu=',num2str(k*ile_klatek/klatki),' sekund'];
%disp(tekst)
tekst=['Œrednia predkoœæ= ',num2str(V_sr),' cm/sek'];
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
plot(1:i-1,V,'b')%(1:length(druk(:,4)),druk(:,4),1:length(druk(:,5)),druk(:,5)) % wyrysowanie danych na podwójnej osi OY
title('prêdkoœc')
%legend('X-Y','V')
%hold on
%line([1 k],[V_sr V_sr],'color','r') % narysowanie wartoœci sredniejpredkoœci
%legend('Prêdkoœæ chwilowa','Prêdkoœæ œrednia','Droga','Location','Best')