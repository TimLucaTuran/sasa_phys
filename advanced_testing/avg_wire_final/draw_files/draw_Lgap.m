function Layer1 = draw_Lgap(L,W,D,P,winkel)

% L=250;
% W=150;
% D=30;
% L=L-D;
% P=300;
% winkel=0;

% linear diaboloid coordinates
% x und y vertauscht denken...

vec=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vec,vec);


radang=winkel;
X=cosd(radang)*x+sind(radang)*y;
Y=-sind(radang)*x+cosd(radang)*y;

% Ruecken
Layer1=ones(size(X));
Layer1(abs(X+D/2)>L/2)=0;
Layer1(X>D/2)=0;
Layer1(abs(Y+W/2-D/2)>D/2)=0;
% Fuss
Layer2=ones(size(X));
Layer2(abs(Y)>W/2)=0;
Layer2(abs(X-L/2)>D/2)=0;
% Layer2=zeros(size(X));

Layer1=Layer1+Layer2;
Layer1(Layer1>1)=1;

% figure(7);clf;imagesc(vec,vec,Layer1);axis image;

% Layer1=fliplr(Layer1);


% figure;imagesc(Layer1);axis image;

% end