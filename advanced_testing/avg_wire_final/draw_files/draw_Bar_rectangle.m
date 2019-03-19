function Layer1 = draw_Bar_rectangle(L,W,P,winkel)
% L=350;
% W=100;
% P=500;
% winkel=0;

% first half of drawing
% x und y vertauscht denken

vecx1=linspace(-P/2,P/2,501);
vecy1=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vecx1,vecy1);

winkel_rad=winkel;

X=x*cosd(winkel_rad)+y*sind(winkel_rad);
Y=-x*sind(winkel_rad)+y*cosd(winkel_rad);

% left bar
Layer1 = zeros(size(Y));
Layer1(X>-L/2 & X<0 & Y<(W/2-P/2)) = 1;
% upper bar
Layer2 = zeros(size(Y));
Layer2(Y>-L/2 & Y<0 & X<(W/2-P/2)) = 1;
% mirror down
Layer3 = flipud(Layer1+Layer2);
% mirror right
Layer4 = fliplr(Layer2+Layer3+Layer1);

Layer1 = Layer1 + Layer2 + Layer3 + Layer4;
% figure(7);imagesc(vecx1,vecy1,Layer1);axis square