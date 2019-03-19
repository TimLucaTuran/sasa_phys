function LayerOUT = draw_backbone_II(n_TiO2,n_G,PX,PY,W_G,W_TiO2)

% n_TiO2 = 2.35;
% n_G    = 1.46;
% P      = 400;%nm
% W_G    = 100;%nm
% W_TiO2 = 20;%nm

W1 = W_G;
W2 = W_G + 2*W_TiO2;

NX = 1001;
NY = 11;

vecx = PX*linspace(1-1/NX,1/NX,NX) - PX/2;
vecy = PY*linspace(1-1/NY,1/NY,NY) - PY/2;

[x,y] = ndgrid(vecx,vecy);

Layer1 = ones(size(x));
Layer1(abs(x)>W1/2) = 0;

Layer2 = ones(size(x));
Layer2(abs(x)>W2/2) = 0;

LayerOUT = Layer1 + Layer2;

LayerOUT(LayerOUT==1) = n_TiO2;
LayerOUT(LayerOUT==0) = 1;
LayerOUT(LayerOUT==2) = n_G;

% figure(11);
% imagesc(LayerOUT)
% colorbar