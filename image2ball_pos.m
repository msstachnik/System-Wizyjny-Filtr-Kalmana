function [x,y] = image2ball_pos(image)
% [x,y] = image2ball_pos(image)
% funkcja do odczytu pozycji pi³eczki na podstawie obrazu
image_filtered = medfilt2(image, [3 3]);
% postac binarna z progiem 0.18
image_bin = im2bw(image_filtered,0.18);

% % domkniecie obiekty mniejsze niz 10pix - u¿ywaæ jeœli obraz jest niestabilny
image_bin_closed = bwareaopen(image_bin,10);
% etykowanie obiektów
% figure(2)
% imshow(image_bin_closed)
% figure(1)
stats = regionprops(image_bin_closed, 'Area', 'Centroid');
liczba_obiektow = length(stats);
if liczba_obiektow == 0 % jeœli nie ma obiektów zwraca not a number
    x = NaN;
    y = NaN;
elseif liczba_obiektow == 1 % jeœli jeden przypisuje wartoœæ centraln¹
    x = stats.Centroid(1);
    y = stats.Centroid(2);
elseif liczba_obiektow > 1 % w innym przypadku wyznaczyæ œrodek masy
    % predefiniowanie
    dx_A(liczba_obiektow) = 0;
    dy_A(liczba_obiektow) = 0;
    % obliczenie iliczynów pozycja srodka masy x powierzchnia
    for i = 1 : liczba_obiektow
        dx_A(i) = stats(i).Centroid(1) * stats(i).Area;
        dy_A(i) = stats(i).Centroid(2) * stats(i).Area;
    end
    % obliczenie pozycji
    x = sum(dx_A) / sum([stats.Area]);
    y = sum(dy_A) / sum([stats.Area]);
else % w specyficznych przypadkach zróæ NaN i warning
    warning('wrong input of image2ball_pos function')
    x = NaN;
    y = NaN;
end

