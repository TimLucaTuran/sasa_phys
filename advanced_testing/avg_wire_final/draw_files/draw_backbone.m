function LayerOUT = draw_backbone(PAR_STRUCT)

% n_G = 1.46;
% PY  = 10;%nm
% PX  = 400;
% W_G = 100;%nm
PX = PAR_STRUCT.p6;
PY = PAR_STRUCT.p7;
W1 = PAR_STRUCT.p1;

pixel_num = 1001;

vecx  = linspace(-PX/2,PX/2,pixel_num);
vecy  = linspace(-PY/2,PY/2,pixel_num);
[x,y] = ndgrid(vecx,vecy);

Layer1 = ones(size(x));
Layer1(abs(x)>W1/2) = 0;

LayerOUT = Layer1;

% figure(4);
% imagesc(vecy,vecx,LayerOUT)
% colorbar