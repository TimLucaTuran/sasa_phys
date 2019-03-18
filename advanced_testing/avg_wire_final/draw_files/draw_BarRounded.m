function LayerOUT = draw_BarRounded(PAR_STRUCT)
% Draws a matrix of 1s and 0s in the shape of a rod/bar/wire with rounded
% cornders. The x-dimension is defined by the width W and the y-dimension
% is defined by the length L. The cornes are each rounded by a quarter
% circle of radius RAD. The periods in x and y which define the unit cell
% are given by PX and PY. The rotation angle in mathematical direction is
% set by ang in degrees.
%
% Usage: LayerOUT = draw_Bar_rounded(PAR_STRUCT)
%
% Input via the struct PAR_STRUCT:
% L   = PAR_STRUCT.L;
% W   = PAR_STRUCT.W;
% ang = PAR_STRUCT.A;
% PX  = PAR_STRUCT.PX;
% PY  = PAR_STRUCT.PY;
% RAD = PAR_STRUCT.radius;

L   = PAR_STRUCT.L;
W   = PAR_STRUCT.W;
ang = PAR_STRUCT.A;
PX  = PAR_STRUCT.PX;
PY  = PAR_STRUCT.PY;
RAD = PAR_STRUCT.radius;

NX = 501;
NY = 501;

vecx  = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy  = PY*linspace(1-1/NY,1/NY,NY) - PY/2;
[x,y] = ndgrid(vecx,vecy);

X =  x * cosd(ang) + y * sind(ang);
Y = -x * sind(ang) + y * cosd(ang);

Layer1 = ones( size(Y) );
Layer1( abs(X) > (W - 2 * RAD) / 2 ) = 0;
Layer1( abs(Y) > L / 2 ) = 0;

Layer2 = ones( size(Y) );
Layer2( abs(X) > W / 2 ) = 0;
Layer2( abs(Y) > (L - 2 * RAD) / 2 ) = 0;

X_shift = W/2 - RAD;
Y_shift = L/2 - RAD;

Layer3 = zeros( size(Y) );
Layer3( (abs(X + X_shift).^2 + abs(Y + Y_shift).^2) < RAD^2 ) = 1;
Layer3( (abs(X - X_shift).^2 + abs(Y + Y_shift).^2) < RAD^2 ) = 1;
Layer3( (abs(X + X_shift).^2 + abs(Y - Y_shift).^2) < RAD^2 ) = 1;
Layer3( (abs(X - X_shift).^2 + abs(Y - Y_shift).^2) < RAD^2 ) = 1;

LayerOUT = Layer1 + Layer2 + Layer3;
LayerOUT( LayerOUT > 1 ) = 1;

% figure(11)
% imagesc(LayerOUT)
% axis image