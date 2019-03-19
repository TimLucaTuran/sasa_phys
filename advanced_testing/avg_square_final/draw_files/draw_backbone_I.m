function LayerOUT = draw_backbone_I(n_TiO2,PX,PY,W_G,W_TiO2)

% n_TiO2 = 2.35;
% PX     = 10;%nm
% PY     = 400;
% W_G    = 100;%nm
% W_TiO2 = 20;%nm

W = W_G + 2*W_TiO2;

NX = 1001;
NY = 11;

vecx = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy = PY*linspace(1-1/NY,1/NY,NY) - PY/2;
[x,y] = ndgrid(vecx,vecy);

Layer1 = ones(size(x));
Layer1(abs(x)>W/2) = 0;

LayerOUT = Layer1*n_TiO2;
LayerOUT(LayerOUT==0) = 1;

% figure(10);
% imagesc(LayerOUT)
% colorbar