function Layer=draw_Square_scaled(PAR_STRUCT)

P     = PAR_STRUCT.p7;
scale = PAR_STRUCT.p2;
W     = P*scale;
ang_  = PAR_STRUCT.p4;

vecx=linspace(-P/2,P/2,501);
vecy=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vecx,vecy);

X=x*cosd(ang_)+y*sind(ang_);
Y=-x*sind(ang_)+y*cosd(ang_);

Layer=zeros(size(X));
Layer(abs(X)>W/2)=1;
Layer(abs(Y)>W/2)=1;
Layer=1-Layer;
