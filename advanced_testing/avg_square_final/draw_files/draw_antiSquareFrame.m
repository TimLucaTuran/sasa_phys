function Layer=draw_antiSquareFrame(PAR_STRUCT)
% Draw a square frame of over all width W and frame girth G.


P    = PAR_STRUCT.P;
W    = PAR_STRUCT.W;
G    = PAR_STRUCT.G;
ang_ = PAR_STRUCT.A;

NX = 501;
NY = 501;

vecx  = P * linspace( 1-1/NX, 1/NX,NX ) - P/2;
vecy  = P * linspace( 1-1/NY, 1/NY,NY ) - P/2;
[x,y] = ndgrid( vecx, vecy );

X = x * cosd(ang_) + y * sind(ang_);
Y = -x * sind(ang_) + y * cosd(ang_);

% draw main square
Layer1 = zeros(size(X));
Layer1(abs(X) > W/2) = 1;
Layer1(abs(Y) > W/2) = 1;

% remove center to create frame
Layer2 = ones(size(X));
Layer2(abs(X) > W/2 - G) = 0;
Layer2(abs(Y) > W/2 - G) = 0;

Layer = Layer1 + Layer2;

% figure(1);
% imagesc(vecx,vecy,Layer)
% colorbar
% axis image