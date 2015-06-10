figure 
clc % szyszczenie okna polece
close all % zamkni�cie wszystkich okien figure
clear all
klatki=30; % sta�a filmu ile klatek na sekund�
odleglosc_wierzcholkow=20; % odleglosc w cm mi�dzy wierzcho�kami kwadratu ma�ego

[filename, pathname] = uigetfile('*.mp4','Wybierz plik z filmem'); %pobranie ��iezki i nazwy pliku
pathname_filename = [pathname, filename];   % sklejenie nazwy isciezki


xyloObj = VideoReader(pathname_filename);
nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;
film(1:nFrames) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'),'colormap',[]);

for k = 1 : nFrames
    film(k).cdata = read(xyloObj,k);
end




%load(pathname_filename)
%load film % wczystanie pliku z filmem



%rozmiar=size(film); % sprawdzenie rozmiaru filmu
disp('Nr kl     X     Y    V[cm/s] S[cm] S_t') % wy�wietlenie nag��wka do wy�wietlania wynik�w

  wsp_x_poprzednie=0; % ustalenie warto�ci poprzedniego po�o�enia dla pierwszego wykrycia kulki  
  wsp_y_poprzednie=0;
  droga=0; % ustalenie warto�ci poczatkowej przebytej drogi
  
  figure % otwarcie okna
  image(film(round(nFrames/2)).cdata) % wy�wietlenie obrazu na kt�rym zaznaczamu mysza 2 punkty na jednym poziomie znajduj�ce si� w dw�ch rogach ma�ego kwadratu
  % dzi�ki temu bedziemy wiedzi� ile pikseli jest pomi�dzy wierzcho�kami
  % kwadratu co odpowiada 20 cm odleglo�ci w rzeczywisto�ci
  [x,y] = ginput(2); % pobranie wsp�rzednych dw�ch punkt�w
  dl_boku=abs(x(2)-x(1));
  stala=abs(x(2)-x(1))/odleglosc_wierzcholkow; % obliczenie sta�ej ile pikseli przypada na 1 cm
  close(gcf) % zamkni�cie otwartego okna
  
  znak=1; % znacznik do wywo�ania wafunki WHILE
  k=0; % zmienna do numeracji kilejnych po�o�e� kulki

  for i=5:nFrames % 2:nFrames % rozmiar 4 to ilo�� klatek, rozmiar ma 4 elemnty [a,b,c,d], a i b wys i szer, c- ilo�� wartw (3 dla RGB), d - ilo�c klatek
    % zaczynami do 2 bo bierzemy klatke wczesniejsz� (2-1) i bierz�c� 2
   %if i = 5 
       % odejmij od go�ej mapy
   %else
       
    obraz_roznica=medfilt2(rgb2gray(film(i-4).cdata))-medfilt2(rgb2gray(film(i).cdata));
    imshow(rgb2gray(film(i-2).cdata))
    wynik=obraz_roznica(5:end-5,5:end-5);%rgb2gray(obraz_roznica); % przepisanie obrazu r�nicowego czarno bia�ego do zmiennej wynik
    hold on % wstrzymanie wykresu, mo�na dorysoswywac inne elementy (znacznik drogi kulki)
    %%%%%%%%%%%%%%%%%%
    aa=wynik>50;
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
    
    if sum(sum(aa))<100 && sum(sum(aa))>5  %max(max(wynik))>50 &&sum(poziom)>500 && sum(pion)>500 % za�o�ono �e powierznia kulki ma mniej ni� 700 pix (w rzeczywisto�ci dla rozdzielczo�ci 800x600 ma ok.440 pix , je�eli na obrazie suma pikseli bia�ych jest mniejsza ni� 700
         % zmiana warto�ci znacznika
        % % narysowanie znacznika wskazuj�cego kulke
        
        
         % szukanie wsp�rzednych kulki, kt�ra w obrazie bw jest bia�a - szukamy
    % miejsca gdzie jest najwi�cej pikseli bia�ych kt�re nale�a do kulki
    poziom=sum((aa)); %sum((wynik)) sumowanie ilo�ci pikseli (bo piksele maj� warto�ci 0 i 1) w poziomie
    wsp_x=find(poziom==max(poziom));% szukanie miejsca indeksu gdzie w wektorze w kt�rym sa sumy dla ka�ej kolumny wynik wyst�puj� maximum
    % wsp_x=mean( wsp_x); opcja gdyby by�y wskazane 2 lub 3 maksima to
    % wybierzemy u�rednion� wsp�rzedn�

    pion=sum(aa');%sum(wynik') sumowanie ilo�ci pikseli (bo piksele maj� warto�ci 0 i 1) w pionie (naprawde w poziomie ale macierz jest obr�cona - transponowana  co daje ten sam efekt)
    wsp_y=find(pion==max(pion));  % szukanie miejsca indeksu gdzie w wektorze w kt�rym sa sumy dla ka�ej kolumny wynik wyst�puj� maximum
    % wsp_y=mean( wsp_y); opcja gdyby by�y wskazane 2 lub 3 maksima to
    % wybierzemy u�rednion� wsp�rzedn�
        
        
        
        while znak==1 % jelsi spe�niony warunek znak=1
            
            
            wsp_x_poprzednie=wsp_x(1); % przepisanie znalezionych wsp�rzendych do zapamietania (w celu w nast�pnym kroku uzycia ich jako poprzednie po�o�enie) , 
            %jest to dla pierwszych znalzionych ws�rzednych
            wsp_y_poprzednie=wsp_y(1);
            znak=0;
        end

        if i>=2 % sprawdzamy 
            
            s=sqrt((wsp_x_poprzednie(1)-wsp_x(1))^2+(wsp_y_poprzednie(1)-wsp_y(1))^2);
                        wsp_x_poprzednie=wsp_x(1);
            wsp_y_poprzednie=wsp_y(1);
            
            if s>0 && s<dl_boku/5 % jesli droga pokonana przez kulk� jest wieksza ni� 1/4 d� bokuto traktujemy jako zak��cenie i nie uwzgl�dniamy
                plot(wsp_x(1),wsp_y(1),'rs')
                k=k+1;
            %line([wsp_x_poprzednie(1), wsp_x(1)],[wsp_y_poprzednie(1), wsp_y(1)],'LineWidth',5,'Color','r');    
            s=s/stala;
            droga=droga+s;
            V=s/(1/klatki);
            %wsp zapamietane % przepisanie znalezionych wsp�rzendych do zapamietania (w celu w nast�pnym kroku uzycia ich jako poprzednie po�o�enie) 
            wsp_x_poprzednie=wsp_x(1);
            wsp_y_poprzednie=wsp_y(1);
            
                  druk(k,:)=[i,wsp_x(1),wsp_y(1),round(V),round(droga),round(s)]; % przygotowanie danych do wy�wietleniai stworzenie macierzy z danymi
                disp(druk(k,:)) % wy�wietlenie berzacych wynik�w
            end
            
            
        end

       % if i>=2
        %wsp zapamietane % przepisanie znalezionych wsp�rzendych do zapamietania (w celu w nast�pnym kroku uzycia ich jako poprzednie po�o�enie) 
       % wsp_x_poprzednie=wsp_x;
       % wsp_y_poprzednie=wsp_y;

      %  end

       % druk(k,:)=[i,wsp_x(1),wsp_y(1),round(V),round(droga),round(s)]; % przygotowanie danych do wy�wietleniai stworzenie macierzy z danymi
        %tekst=[num2str(i),'      ',num2str(wsp_x(1)),'     ',num2str(wsp_y(1)), '      ',num2str(round(V)), '      ',num2str(round(droga)), '           ',num2str(round(s))];
        %disp(druk(k,:)) % wy�wietlenie berzacych wynik�w
        drawnow % polecenie odswiezenia wykresu
    end
    
    hold off
    
end

% wyliczenie warto�ci �redniej
ile_klatek=k;
V_sr=droga/(ile_klatek*1/klatki); 

disp('==================================')
disp('Warto�c �rednia')
tekst=['Czas ruchu=',num2str(ile_klatek*1/klatki),' sekund'];
disp(tekst)
tekst=['�rednia predko��= ',num2str(V_sr),' cm/sek'];
disp(tekst)
tekst=['Przebyty dustans= ',num2str(droga), ' cm'];
disp(tekst)
close gcf


figure

plotyy(1:length(druk(:,4)),druk(:,4),1:length(druk(:,5)),druk(:,5)) % wyrysowanie danych na podw�jnej osi OY

hold on
line([1 k],[V_sr V_sr],'color','r') % narysowanie warto�ci sredniejpredko�ci
legend('Pr�dko�� chwilowa','Pr�dko�� �rednia','Droga','Location','Best')


%++++++++++++++++++++++++++++++++++++++++++++++++++++
% predko��

t=(druk(:,1)-druk(1,1))*1/( xyloObj.FrameRate);
v=medfilt1(druk(:,4),5); % filtracja medianowa pr�dko�ci dla 5 warto�ci


%figure
f = fit(t(2:end),v(2:end),'poly1');
%plot(f,t(2:end),v(2:end));
a=f.p1;
b=f.p2;
wynik=roots([a,b]);


x=0:abs(wynik);
y=a*x+b;

figure
plot(t(2:end),v(2:end),'.r');
hold on
plot(x,y)
title('Estymacja pr�dko�ci')

% =========================================================
% po�ozenie
xx=druk(:,2);%2);
yy=druk(:,3);%3);3

figure
f = fit(xx,yy,'poly1');
plot(f,xx,yy);
a=f.p1;
b=f.p2;
%wynik=roots([a,b]);
title('Wykres po�ozenia')


x=0:1000;
y=a*x+b;

figure
plot(xx,yy,'.r')
hold on
plot(x,y)
title('Estymacja po�ozenia')


