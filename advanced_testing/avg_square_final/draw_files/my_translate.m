function TLayer = my_translate(Layer_,pixel,TI_nm,TJ_nm)
% Translates layers (layer matrices) for draw_ files.
% TI_nm and TJ_nm give the translation length in nm. I... row and J
% ...column
% A scale is needed to convert nm into pixels: pixel = scale/Period.
% CAUTION: This method works only matrices with values 0 and 1.


% transform nm to pixels
TI = round(TI_nm/pixel);
TJ = round(TJ_nm/pixel);

TLayer = Layer_;
[LI,LJ] = find(TLayer==1);
LIND    = sub2ind(size(TLayer),LI,LJ);
LIND_T  = sub2ind(size(TLayer),LI+TI,LJ+TJ);
TLayer(LIND)   = 0;
TLayer(LIND_T) = 1;