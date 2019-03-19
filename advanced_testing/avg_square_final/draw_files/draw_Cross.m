function Layer=draw_Cross(PAR_STRUCT)

% PX  = 300;
% PY  = PX;
% L   = PX/1.6;
% W   = 0.24*L;
% ang = 0;

L   = PAR_STRUCT.L;
W   = PAR_STRUCT.W;
ang = PAR_STRUCT.A;
PX  = PAR_STRUCT.PX;
PY  = PAR_STRUCT.PY;

vecx  = linspace(-PX/2,PX/2,501);
vecy  = linspace(-PY/2,PY/2,501);
[x,y] = ndgrid(vecx,vecy);

X = x*cosd(ang)+y*sind(ang);
Y = -x*sind(ang)+y*cosd(ang);

Layer1=ones(size(X));
Layer1(abs(X)>L/2)=0;
Layer1(abs(Y)>W/2)=0;
Layer2=ones(size(X));
Layer2(abs(Y)>L/2)=0;
Layer2(abs(X)>W/2)=0;

Layer=Layer1+Layer2;
Layer(Layer>1)=1;

% figure(1);imagesc(vecx,vecy,Layer);axis image;

