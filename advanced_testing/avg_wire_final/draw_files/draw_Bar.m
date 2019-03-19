function LayerOUT = draw_Bar(PAR_STRUCT)

L   = PAR_STRUCT.L;
W   = PAR_STRUCT.W;
ang = PAR_STRUCT.A;
PX  = PAR_STRUCT.PX;
PY  = PAR_STRUCT.PY;

NX = 501;
NY = 501;

vecx  = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy  = PY*linspace(1-1/NY,1/NY,NY) - PY/2;
[x,y] = ndgrid(vecx,vecy);

X = x*cosd(ang)+y*sind(ang);
Y = -x*sind(ang)+y*cosd(ang);

Layer1 = ones(size(Y));
Layer1(abs(X)>W/2) = 0;
Layer1(abs(Y)>L/2) = 0;

LayerOUT = Layer1;