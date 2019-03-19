function Layer = draw_super9_Square(PAR_STRUCT)
% Creates a super cell of squares containing 9 unit cells. This only makes
% sense in combination with a second periodic layer, where the period ratio
% from this to the other is 2/3.

PX   = PAR_STRUCT.PX;
PY   = PAR_STRUCT.PY;
W    = PAR_STRUCT.W_Sq;
ang_ = PAR_STRUCT.A;

PX   = PX/3;
PY   = PY/3;

NX = 501;
NY = 501;

vecx  = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy  = PY*linspace(1-1/NY,1/NY,NY) - PY/2;
[x,y] = ndgrid(vecx,vecy);

X = x * cosd(ang_) + y * sind(ang_);
Y = -x * sind(ang_) + y * cosd(ang_);

Layer = zeros(size(X));

Layer( abs(X) > W/2 ) = 1;
Layer( abs(Y) > W/2 ) = 1;

Layer = 1 - Layer;
Layer = repmat(Layer,[1,3]);
Layer = repmat(Layer,[3,1]);