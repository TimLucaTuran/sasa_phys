function Layer=draw_antiSquare(PAR_STRUCT)

P    = PAR_STRUCT.P;
W    = PAR_STRUCT.W;
ang_ = PAR_STRUCT.A;

vecx=linspace(-P/2,P/2,501);
vecy=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vecx,vecy);

X=x*cosd(ang_)+y*sind(ang_);
Y=-x*sind(ang_)+y*cosd(ang_);

Layer=zeros(size(X));
Layer(abs(X)>W/2)=1;
Layer(abs(Y)>W/2)=1;

% figure(1);
% imagesc(vecx,vecy,Layer)
% colorbar