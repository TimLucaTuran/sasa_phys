function Layer1 = draw_Bar_sqspecial(L,W,P,winkel)
% L=300;
% W=80;
% P=450;
% winkel=0;

% first half of drawing
% x uny vertauscht denken

vecx1=linspace(-P/2,P/2,501);
vecy1=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vecx1,vecy1);

winkel_rad=winkel;

X=x*cosd(winkel_rad)+y*sind(winkel_rad);
Y=-x*sind(winkel_rad)+y*cosd(winkel_rad);

% central half-bar
Layer1=zeros(size(Y));
Layer1(abs(X)>L/2)=1;
Layer1(abs(Y)>W/2)=1;
Layer1=1-Layer1;
% upper-right corner bar
Layer2=ones(size(y));
Layer2(X>-P/2+L/2)=0;
Layer2(Y<P/2-W/2)=0;
% lower-right corner bar
Layer3 = flipud(Layer2);
% lower-right corner bar
Layer4 = fliplr(Layer2+Layer3);

Layer1 = Layer1 + Layer2 + Layer3 + Layer4;
% figure(7);imagesc(Layer1);axis image;