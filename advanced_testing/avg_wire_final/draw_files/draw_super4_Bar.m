function LayerOUT = draw_super4_Bar(PAR_STRUCT)
% Creates a super cell of wires/bars containing 9 unit cells. This only makes
% sense in combination with a second periodic layer, where the period ratio
% from this to the other is 3/2.

L   = PAR_STRUCT.L_Bar;
W   = PAR_STRUCT.W_Bar;
ang = PAR_STRUCT.A;
PX  = PAR_STRUCT.PX;
PY  = PAR_STRUCT.PY;

PX   = PX/2;
PY   = PY/2;

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

Layer1 = repmat(Layer1,[1,2]);
Layer1 = repmat(Layer1,[2,1]);

LayerOUT = Layer1;