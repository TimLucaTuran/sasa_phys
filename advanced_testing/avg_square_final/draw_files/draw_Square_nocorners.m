function LayerOUT = draw_Square_nocorners(PAR_STRUCT)

P   = PAR_STRUCT.p7;
ang = PAR_STRUCT.p4;
W   = PAR_STRUCT.p2;

vecx  = linspace(-P/2,P/2,501);
vecy  = linspace(-P/2,P/2,501);
[x,y] = ndgrid(vecx,vecy);

X  =x*cosd(ang)+y*sind(ang);
Y = -x*sind(ang)+y*cosd(ang);

%% draw square
Layer1 = zeros(size(X));
% Layer1(abs(X)<W/2 & abs(Y)<W/2) = 1;
Layer1(abs(X)>W/2 & abs(Y)>W/2) = 10;
%% remove corners
% first: define square to be removed
Layer2_ = zeros(size(X));
Layer2_(abs(X)<W/8 & abs(Y)<W/8) = 1;
% translation params in nm
shift1 = 3*W/8;
Layer2_T1 = my_translate(Layer2_,P/500,shift1,shift1);
Layer2_T2 = my_translate(Layer2_,P/500,-1*shift1,shift1);
Layer2_T3 = my_translate(Layer2_,P/500,shift1,-1*shift1);
Layer2_T4 = my_translate(Layer2_,P/500,-1*shift1,-1*shift1);
% creat overlay layer
Layer2 = Layer2_T1 + Layer2_T2 + Layer2_T3 + Layer2_T4;
% substract overlay layer
Layer3 = Layer1 - Layer2;

%% add circle to make rounded corners
% draw circle
circle_ = zeros(size(x));
circle_(abs(X).^2+abs(Y).^2<abs(W).^2/16) = 1;
% move circle
shift2 = W/4;
circle_T1 = my_translate(circle_,P/500,shift2,shift2);
circle_T2 = my_translate(circle_,P/500,-1*shift2,shift2);
circle_T3 = my_translate(circle_,P/500,shift2,-1*shift2);
circle_T4 = my_translate(circle_,P/500,-1*shift2,-1*shift2);
% creat overlay layer
Layer4 = circle_T1 + circle_T2 + circle_T3 + circle_T4;
% add overlay layer
LayerOUT = Layer3 + Layer4;
LayerOUT(LayerOUT ~= 10 & LayerOUT > 0 ) = 1;
LayerOUT(LayerOUT == 10) = 0;
LayerOUT(LayerOUT < 0) = 0;