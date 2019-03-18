function LayerOUT = draw_antiBar(PAR_STRUCT)

L   = PAR_STRUCT.L;
W   = PAR_STRUCT.W;
ang = PAR_STRUCT.A;
PX  = PAR_STRUCT.PX;
PY  = PAR_STRUCT.PY;

NX = 501;
NY = NX;
vecx = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy = PY*linspace(1-1/NY,1/NY,NY) - PY/2;
[x,y] = ndgrid(vecx,vecy);

X = x*cosd(ang)+y*sind(ang);
Y = -x*sind(ang)+y*cosd(ang);

Layer1 = zeros(size(Y));
Layer1(abs(X)>W/2) = 1;
Layer1(abs(Y)>L/2) = 1;

LayerOUT = Layer1;

% figure(1);
% imagesc(vecx,vecy,LayerOUT)
% colorbar