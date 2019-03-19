function Layer1 = draw_plain(PAR_STRUCT)
% draw a v lying on on side; open to the left

P     = PAR_STRUCT.p7;

% opening angle of V

vec=linspace(-P/2,P/2,501);
[x,y]=ndgrid(vec,vec);

Layer1 = ones(size(x));

