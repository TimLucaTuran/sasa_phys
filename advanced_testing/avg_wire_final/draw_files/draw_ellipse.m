function Layer_out = draw_ellipse(PAR_STRUCT)
% Draws a matrix of 1s and 0s in the shape of an ellipse. The x-dimension 
% is defined as the Minor axis 'minor' and the y-dimension
% is defined asthe major axis 'major'. The periods in x and y which define
% the unit cell are given by PX and PY. The rotation angle in mathematical
% direction is set by ang in degrees.
%
% Usage: LayerOUT = draw_ellipse(PAR_STRUCT)
%
% Input via the struct PAR_STRUCT:
% PX    = PAR_STRUCT.PX;
% PY    = PAR_STRUCT.PY;
% Major = PAR_STRUCT.Major;
% Minor = PAR_STRUCT.Minor;
% ang_  = PAR_STRUCT.A;

PX    = PAR_STRUCT.PX;
PY    = PAR_STRUCT.PY;
Major = PAR_STRUCT.Major / 2;
Minor = PAR_STRUCT.Minor / 2;
ang_  = PAR_STRUCT.A;

% Test parameters
% P     = 300;
% Major = 120;
% Minor = 70;
% ang_  = 0;

vecx  = linspace(-PX/2,PX/2,501);
vecy  = linspace(-PY/2,PY/2,501);
[x,y] = ndgrid(vecx,vecy);

X =  x * cosd(ang_) + y * sind(ang_);
Y = -x * sind(ang_) + y * cosd(ang_);

Layer_out = zeros(size(X));

Layer_out( abs(X).^2 / Minor^2 + abs(Y).^2 / Major^2 <= 1 ) = 1;

% figure(1)
% imagesc(vecx, vecy, Layer_out)
% axis image
% colorbar