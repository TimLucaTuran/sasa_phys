function Layer1 = draw_split_circ_sym3010(PAR_STRUCT)
% draw a v lying on on side; open to the left

P   = PAR_STRUCT.p7;
D   = PAR_STRUCT.p1;
W   = PAR_STRUCT.p2;
ANG = PAR_STRUCT.p4;
R1  = D/2;
R2  = R1 + W;

ANG1 = 30*pi/180;
ANG2 = 10*pi/180;

vec = linspace(-P/2,P/2,501);
[x,y] = ndgrid(vec,vec);

radang = ANG;
X = cosd(radang)*x+sind(radang)*y;
Y = -sind(radang)*x+cosd(radang)*y;

% draw ring #1
Layer1 = zeros(size(x));
Layer1( (X.^2+Y.^2) > R1^2 & (X.^2+Y.^2) < R2^2 ) = 1;

% cut out #1
Layer1( Y./X > tan(-1*ANG1/2) & Y./X < tan(1*ANG1/2) ) = 0;
Layer1( X > 0 ) = 0;

% draw ring #2
Layer2 = zeros(size(x));
Layer2( (X.^2+Y.^2) > R1^2 & (X.^2+Y.^2) < R2^2 ) = 1;

% cut out #2
Layer2( Y./X > tan(-1*ANG2/2) & Y./X < tan(1*ANG2/2) ) = 0;
Layer2( X < 0 ) = 0;

Layer1 = Layer1 + Layer2;
Layer1(Layer1 > 1) = 1;

