clc % 
close all % zamkniêcie wszystkich okien figure
clear all
%% okreœlnie parametrów
klatki=30; % sta³a filmu ile klatek na sekundê
dt = 1 / klatki; % okreslenie podstawy czasu
odleglosc_wierzcholkow=20; % odleglosc w cm miêdzy wierzcho³kami kwadratu ma³ego


%% pobranie obrazu
[filename, pathname] = uigetfile('*.mp4','Wybierz plik z filmem'); %pobranie œæiezki i nazwy pliku
pathname_filename = [pathname, filename];   % sklejenie nazwy isciezki

xyloObj = VideoReader(pathname_filename);
nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;
film(1:nFrames) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'),'colormap',[]);

for k = 1 : nFrames
    film(k).cdata = read(xyloObj,k);
end


wsp_x_poprzednie=0; % ustalenie wartoœci poprzedniego po³o¿enia dla pierwszego wykrycia kulki  
wsp_y_poprzednie=0;
droga=0; % ustalenie wartoœci poczatkowej przebytej drogi

figure % otwarcie okna
image(film(round(nFrames/2)).cdata) % wyœwietlenie obrazu na którym zaznaczamu mysza 2 punkty na jednym poziomie znajduj¹ce siê w dwóch rogach ma³ego kwadratu
% dziêki temu bedziemy wiedziæ ile pikseli jest pomiêdzy wierzcho³kami
% kwadratu co odpowiada 20 cm odlegloœci w rzeczywistoœci
[x,y] = ginput(2); % pobranie wspó³rzednych dwóch punktów
dl_boku=abs(x(2)-x(1));
stala=abs(x(2)-x(1))/odleglosc_wierzcholkow; % obliczenie sta³ej ile pikseli przypada na 1 cm
close(gcf) % zamkniêcie otwartego okna

znak=1; % znacznik do wywo³ania wafunki WHILE
k=0; % zmienna do numeracji kilejnych po³o¿eñ kulki

  for i=5:nFrames % 2:nFrames % rozmiar 4 to iloœæ klatek, rozmiar ma 4 elemnty [a,b,c,d], a i b wys i szer, c- iloœæ wartw (3 dla RGB), d - iloœc klatek
    % zaczynami do 2 bo bierzemy klatke wczesniejsz¹ (2-1) i bierz¹c¹ 2
   %if i = 5 
       % odejmij od go³ej mapy
   %else
    
    data = film(i-4).cdata;
    n =50;
    
    obraz_roznica=medfilt2(rgb2gray(film(i-4).cdata))-medfilt2(rgb2gray(film(i).cdata));

    tic
    data = data(1:n,1:n,:);
    diff_im_1 = imsubtract(data(:,:,1), rgb2gray(data));
    %filtracja medianowa w celu usuniêcia ewentualnych zak³óceñ 
    diff_im_1 = medfilt2(diff_im_1, [3 3]);
    % postac binarna z progiem 0.18
    diff_im_1 = im2bw(diff_im_1,0.18);
    
    % domkniecie obiekty mniejsze niz 10pix
    diff_im_1 = bwareaopen(diff_im_1,10);
    
    % etykietowanie obiektow
    bw = bwlabel(diff_im_1, 8);
%     stats = regionprops(diff_im_1, 'BoundingBox', 'Centroid');
    stats = regionprops(bw, 'BoundingBox', 'Centroid');
    time(i) = toc;
    imshow(rgb2gray(film(i-2).cdata))
    wynik=obraz_roznica(5:end-5,5:end-5);%rgb2gray(obraz_roznica); % przepisanie obrazu ró¿nicowego czarno bia³ego do zmiennej wynik
    hold on % wstrzymanie wykresu, mo¿na dorysoswywac inne elementy (znacznik drogi kulki)
    %%%%%%%%%%%%%%%%%%
    aa=wynik>50;
    %%imshow(aa)
    
    
    
    %%%%%%%%%%%%%%%%%%
    % szukanie wspó³rzednych kulki, która w obrazie bw jest bia³a - szukamy
    % miejsca gdzie jest najwiêcej pikseli bia³ych które nale¿a do kulki
    %poziom=sum(aa); %sum((wynik)) sumowanie iloœci pikseli (bo piksele maj¹ wartoœci 0 i 1) w poziomie
    %wsp_x=find(poziom==max(poziom));% szukanie miejsca indeksu gdzie w wektorze w którym sa sumy dla ka¿ej kolumny wynik wystêpujê maximum
    % wsp_x=mean( wsp_x); opcja gdyby by³y wskazane 2 lub 3 maksima to
    % wybierzemy uœrednion¹ wspó³rzedn¹

    %pion=sum(aa');%sum(wynik') sumowanie iloœci pikseli (bo piksele maj¹ wartoœci 0 i 1) w pionie (naprawde w poziomie ale macierz jest obrócona - transponowana  co daje ten sam efekt)
    %wsp_y=find(pion==max(pion));  % szukanie miejsca indeksu gdzie w wektorze w którym sa sumy dla ka¿ej kolumny wynik wystêpujê maximum
    % wsp_y=mean( wsp_y); opcja gdyby by³y wskazane 2 lub 3 maksima to
    % wybierzemy uœrednion¹ wspó³rzedn¹
    
    if sum(sum(aa))<100 && sum(sum(aa))>5  %max(max(wynik))>50 &&sum(poziom)>500 && sum(pion)>500 % za³o¿ono ¿e powierznia kulki ma mniej ni¿ 700 pix (w rzeczywistoœci dla rozdzielczoœci 800x600 ma ok.440 pix , je¿eli na obrazie suma pikseli bia³ych jest mniejsza ni¿ 700
         % zmiana wartoœci znacznika
        % % narysowanie znacznika wskazuj¹cego kulke
        
        
         % szukanie wspó³rzednych kulki, która w obrazie bw jest bia³a - szukamy
    % miejsca gdzie jest najwiêcej pikseli bia³ych które nale¿a do kulki
    poziom=sum((aa)); %sum((wynik)) sumowanie iloœci pikseli (bo piksele maj¹ wartoœci 0 i 1) w poziomie
    wsp_x=find(poziom==max(poziom));% szukanie miejsca indeksu gdzie w wektorze w którym sa sumy dla ka¿ej kolumny wynik wystêpujê maximum
    % wsp_x=mean( wsp_x); opcja gdyby by³y wskazane 2 lub 3 maksima to
    % wybierzemy uœrednion¹ wspó³rzedn¹

    pion=sum(aa');%sum(wynik') sumowanie iloœci pikseli (bo piksele maj¹ wartoœci 0 i 1) w pionie (naprawde w poziomie ale macierz jest obrócona - transponowana  co daje ten sam efekt)
    wsp_y=find(pion==max(pion));  % szukanie miejsca indeksu gdzie w wektorze w którym sa sumy dla ka¿ej kolumny wynik wystêpujê maximum
    % wsp_y=mean( wsp_y); opcja gdyby by³y wskazane 2 lub 3 maksima to
    % wybierzemy uœrednion¹ wspó³rzedn¹
        
        
        
        while znak==1 % jelsi spe³niony warunek znak=1
            
            
            wsp_x_poprzednie=wsp_x(1); % przepisanie znalezionych wspó³rzendych do zapamietania (w celu w nastêpnym kroku uzycia ich jako poprzednie po³o¿enie) , 
            %jest to dla pierwszych znalzionych wsó³rzednych
            wsp_y_poprzednie=wsp_y(1);
            znak=0;
        end

        if i>=2 % sprawdzamy 
            
            s=sqrt((wsp_x_poprzednie(1)-wsp_x(1))^2+(wsp_y_poprzednie(1)-wsp_y(1))^2);
                        wsp_x_poprzednie=wsp_x(1);
            wsp_y_poprzednie=wsp_y(1);
            
            if s>0 && s<dl_boku/5 % jesli droga pokonana przez kulkê jest wieksza ni¿ 1/4 d³ bokuto traktujemy jako zak³ócenie i nie uwzglêdniamy
                plot(wsp_x(1),wsp_y(1),'rs')
                k=k+1;
            %line([wsp_x_poprzednie(1), wsp_x(1)],[wsp_y_poprzednie(1), wsp_y(1)],'LineWidth',5,'Color','r');    
            s=s/stala;
            droga=droga+s;
            V=s/(1/klatki);
            %wsp zapamietane % przepisanie znalezionych wspó³rzendych do zapamietania (w celu w nastêpnym kroku uzycia ich jako poprzednie po³o¿enie) 
            wsp_x_poprzednie=wsp_x(1);
            wsp_y_poprzednie=wsp_y(1);
            
                  druk(k,:)=[i,wsp_x(1),wsp_y(1),round(V),round(droga),round(s)]; % przygotowanie danych do wyœwietleniai stworzenie macierzy z danymi
                  disp(druk(k,:)) % wyœwietlenie berzacych wyników
            end
            
            
        end

       % if i>=2
        %wsp zapamietane % przepisanie znalezionych wspó³rzendych do zapamietania (w celu w nastêpnym kroku uzycia ich jako poprzednie po³o¿enie) 
       % wsp_x_poprzednie=wsp_x;
       % wsp_y_poprzednie=wsp_y;

      %  end

       % druk(k,:)=[i,wsp_x(1),wsp_y(1),round(V),round(droga),round(s)]; % przygotowanie danych do wyœwietleniai stworzenie macierzy z danymi
        %tekst=[num2str(i),'      ',num2str(wsp_x(1)),'     ',num2str(wsp_y(1)), '      ',num2str(round(V)), '      ',num2str(round(droga)), '           ',num2str(round(s))];
        %disp(druk(k,:)) % wyœwietlenie berzacych wyników
        drawnow % polecenie odswiezenia wykresu
    end

    hold off
    
  end

  hist(time(7:end),20)
   median(time)
% 
% % wyliczenie wartoœci œredniej
% ile_klatek=k;
% V_sr=droga/(ile_klatek*1/klatki); 
% 
% disp('==================================')
% disp('Wartoœc œrednia')
% tekst=['Czas ruchu=',num2str(ile_klatek*1/klatki),' sekund'];
% disp(tekst)
% tekst=['Œrednia predkoœæ= ',num2str(V_sr),' cm/sek'];
% disp(tekst)
% tekst=['Przebyty dustans= ',num2str(droga), ' cm'];
% disp(tekst)
% close gcf
% 
% 
% figure
% 
% plotyy(1:length(druk(:,4)),druk(:,4),1:length(druk(:,5)),druk(:,5)) % wyrysowanie danych na podwójnej osi OY
% 
% hold on
% line([1 k],[V_sr V_sr],'color','r') % narysowanie wartoœci sredniejpredkoœci
% legend('Prêdkoœæ chwilowa','Prêdkoœæ œrednia','Droga','Location','Best')
% 
% 
% %++++++++++++++++++++++++++++++++++++++++++++++++++++
% % predkoœæ
% 
% t=(druk(:,1)-druk(1,1))*1/( xyloObj.FrameRate);
% v=medfilt1(druk(:,4),5); % filtracja medianowa prêdkoœci dla 5 wartoœci
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
% title('Estymacja prêdkoœci')
% 
% % =========================================================
% % po³ozenie
% xx=druk(:,2);%2);
% yy=druk(:,3);%3);3
% 
% figure
% f = fit(xx,yy,'poly1');
% plot(f,xx,yy);
% a=f.p1;
% b=f.p2;
% %wynik=roots([a,b]);
% title('Wykres po³ozenia')
% 
% 
% x=0:1000;
% y=a*x+b;
% 
% figure
% plot(xx,yy,'.r')
% hold on
% plot(x,y)
% title('Estymacja po³ozenia')
% 

