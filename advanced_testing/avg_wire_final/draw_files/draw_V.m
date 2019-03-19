function Layer1 = draw_V(L,W,D,P,winkel)
% draw a v lying on on side; open to the left

% opening angle of V
AngIn = 40;

vec=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vec,vec);


radang=winkel;
X=cosd(radang)*x+sind(radang)*y;
Y=-sind(radang)*x+cosd(radang)*y;

% foot
Layer1 = ones(size(X));
Layer1(abs(Y)>W/2)=0;
Layer1(abs(X-L/2)>D/2)=0;

% head
radang = AngIn*pi/180;
shift1 = (W/2-D)*(tan(radang).^2);
shift2 = (W/2-D)*tan(radang);
Xh = cosd(radang)*(X+shift1)-sind(radang)*(Y-shift2);
Yh = sind(radang)*(X+shift1)+cosd(radang)*(Y-shift2);

% move head vertically
% maybe readjust for future purposes
% coord_ = D/3;
% Xh = Xh - coord_/cos(radang);
% Yh = Yh + coord_/sin(radang);
% 
% move head horizontally
% coord_2 = D/4;
% Xh = Xh + coord_2/sin(radang);
% Yh = Yh + coord_2/cos(radang);

% add layers together
Layer2 = ones(size(Xh));
Layer2(abs(Yh)>W/2)=0;
Layer2(abs(Xh-L/2)>D/2)=0;

Layer1=Layer1+Layer2;
Layer1(Layer1>1)=1;
