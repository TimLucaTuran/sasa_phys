function Layer=draw_Square(PAR_STRUCT)

PX   = PAR_STRUCT.PX;
PY   = PAR_STRUCT.PY;
W    = PAR_STRUCT.W;
ang  = PAR_STRUCT.A;

NX = 501;
NY = 501;

vecx  = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy  = PY*linspace(1-1/NY,1/NY,NY) - PY/2;
[x,y] = ndgrid(vecx,vecy);

X = x*cosd(ang)+y*sind(ang);
Y = -x*sind(ang)+y*cosd(ang);

Layer = zeros(size(X));
Layer(abs(X)>W/2) = 1;
Layer(abs(Y)>W/2) = 1;
Layer = 1 - Layer;
