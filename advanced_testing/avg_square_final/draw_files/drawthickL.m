function Layer1 = drawthickL(L,W,P,winkel)

% OLD VERSION: Use draw_L.m instead!

% L=240;
% W=40;
% P=300;
% winkel=135;

% linear diaboloid coordinates

vec=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vec,vec);


radang=winkel;
X=cosd(radang)*x+sind(radang)*y;
Y=-sind(radang)*x+cosd(radang)*y;

Layer1=ones(size(X));
Layer1(abs(Y)>L/3)=0;
Layer1(abs(X-L/2+W)>W)=0;

Layer2=ones(size(X));
Layer2(abs(X)>L/2)=0;
Layer2(abs(Y+L/3-W)>W)=0;

Layer1=Layer1+Layer2;
Layer1(Layer1>1)=1;

% figure;imagesc(Layer1);axis image;

% Layer1=fliplr(Layer1);


% figure;imagesc(Layer1);axis image;

% end