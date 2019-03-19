function LayerOUT = draw_L(PAR_STRUCT)

L = PAR_STRUCT.p1;
W = PAR_STRUCT.p2;
D = PAR_STRUCT.p3;
P = PAR_STRUCT.p7;
ang = PAR_STRUCT.p4;

vec=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vec,vec);

X=cosd(ang)*x+sind(ang)*y;
Y=-sind(ang)*x+cosd(ang)*y;

% Back
Layer1=ones(size(X));
Layer1(abs(X)>L/2)=0;
Layer1(abs(Y+W/2-D/2)>D/2)=0;
% Foot
Layer2=ones(size(X));
Layer2(abs(Y)>W/2)=0;
Layer2(abs(X-L/2+D/2)>D/2)=0;

LayerOUT=Layer1+Layer2;
LayerOUT(LayerOUT>1)=1;
