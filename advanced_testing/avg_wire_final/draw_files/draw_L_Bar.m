function LayerOUT = draw_L_Bar(PAR_STRUCT)

L = PAR_STRUCT.p1;
W = PAR_STRUCT.p2;
D = PAR_STRUCT.p3;
P = PAR_STRUCT.p7;
ang = PAR_STRUCT.p4;

% translation params in nm
T3I_nm = 0;
T3J_nm = -75;
T4I_nm = 0;
T4J_nm = T3J_nm + 250;
% transform to pixels
pixel = P/500;
T3I = floor(T3I_nm/pixel);
T3J = floor(T3J_nm/pixel);
T4I = floor(T4I_nm/pixel);
T4J = floor(T4J_nm/pixel);

% x und y vertauscht denken...

vec=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vec,vec);

X=cosd(ang)*x+sind(ang)*y;
Y=-sind(ang)*x+cosd(ang)*y;
%% draw L
% Back
Layer1=ones(size(X));
Layer1(abs(X)>L/2)=0;
Layer1(abs(Y+W/2-D/2)>D/2)=0;
% Foot
Layer2=ones(size(X));
Layer2(abs(Y)>W/2)=0;
Layer2(abs(X-L/2+D/2)>D/2)=0;

Layer3=Layer1+Layer2;
Layer3(Layer3>1)=1;

% translate
[LI3,LJ3] = find(Layer3==1);
LIND    = sub2ind(size(Layer3),LI3,LJ3);
LIND_T  = sub2ind(size(Layer3),LI3+T3I,LJ3+T3J);
Layer3(LIND)   = 0;
Layer3(LIND_T) = 1;

%% draw Bar
% Use L Back and translate
Layer4 = Layer1;
[LI4,LJ4] = find(Layer4==1);
LIND    = sub2ind(size(Layer4),LI4,LJ4);
LIND_T  = sub2ind(size(Layer4),LI4+T4I,LJ4+T4J);
Layer4(LIND)   = 0;
Layer4(LIND_T) = 1;

LayerOUT = Layer3 + Layer4;