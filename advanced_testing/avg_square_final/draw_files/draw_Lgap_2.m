function Layer1 = draw_Lgap_2(L,W,D,P,winkel)

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
