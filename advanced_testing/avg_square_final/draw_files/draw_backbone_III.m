function LayerOUT = draw_backbone_III(n_TiO2,n_G,PX,PY,W_G)

% n_TiO2 = 2.35;
% n_G    = 1.46;
% P      = 400;%nm
% W_G    = 100;%nm
% W_TiO2 = 20;%nm

W1 = W_G;

NX = 1001;
NY = 11;

vecx = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy = PY*linspace(1-1/NY,1/NY,NY) - PY/2;

[x,y] = ndgrid(vecx,vecy);

Layer1 = ones(size(x))*n_G;
Layer1(abs(x)>W1/2) = n_TiO2;

LayerOUT = Layer1;

% figure(13);
% imagesc(LayerOUT)
% colorbar