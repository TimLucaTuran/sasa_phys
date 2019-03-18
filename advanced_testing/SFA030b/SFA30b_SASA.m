%% load SMATs
clear
load1 = 'Square_single_13_CVD';
load2 = 'Bar_single_2';

if ispc
    file_path_1 = ['X:\Research\Programming\MatLab\Calculation\func_MS_design\',load1,'\',load1,'_Daten_gesamt.mat'];
    file_path_2 = ['X:\Research\Programming\MatLab\Calculation\func_MS_design\',load2,'\',load2,'_Daten_gesamt.mat'];
elseif ismac
    file_path_1 = ['/Volumes/sperrhake/Research/Programming/MatLab/Calculation/func_MS_design/',load1,'/',load1,'_Daten_gesamt.mat'];
    file_path_2 = ['/Volumes/sperrhake/Research/Programming/MatLab/Calculation/func_MS_design/',load2,'/',load2,'_Daten_gesamt.mat'];
end
load(file_path_1)
MAT1 = SMAT_;
clear SMAT_;
load(file_path_2)
MAT2 = SMAT_;
%% perform SASA

% save name of mat-file
STORAGE_NAME = 'SFA30b_SASA_dat';

NP2 = 2; NW = 9;

% Smat of gold squares with 500 nm CVD glass on top
MAT1_ = MAT1;
% Smat of gold wires embedded in glass
MAT2_ = squeeze(MAT2(NP2,NW,:,:,:));

figh1 = 10;

if ispc
    PATH = 'X:\Research\Programming\MatLab\Calculation\calc\SASA';
    SAVE_PATH = ...
        ['X:\Research\Programming\MatLab\Calculation\func_MS_design\SFA030b\',STORAGE_NAME];
elseif ismac
    PATH = '/Volumes/sperrhake/Research/Programming/MatLab/Calculation/calc/SASA';
    SAVE_PATH = ...
        ['/Volumes/sperrhake/Research/Programming/MatLab/Calculation/func_MS_design/SFA030b/',STORAGE_NAME];
end

% measured parameters

h_wire    = 41.08;           % Au wire height in nm
h_sq      = 58.2;            % Au square height in nm 
dist_FIB  = 520.4 - h_wire;  % distance from FIB measurement in nm
dist_clad = 575.2 - h_sq;    % height of cladding layer CVDed on top

% define SASA parameters

MATFILE_ = {MAT1_,MAT2_};               % input single layer smats
ANG_     = [0 0];                       % layer rotation angles
MIRROR_  = [0 0];                       % layer mirroring on/off
step_dist= dist_FIB/100;
min_dist = 75*step_dist;
max_dist = 125*step_dist;
dist_vec = (min_dist:step_dist:max_dist)/1000;% layer separation
Hnormal  = {0, h_sq, dist_vec, h_wire, 0};  % system height profile

% input into SASA parameter struct

SASA_PAR.lambda_vec = 1./linspace(1/1.7,1/0.47,64); % wavelength in nm
SASA_PAR.n_subs     = 1.42;                         % substrate index
SASA_PAR.n_clad     = 1;                            % cladding index
SASA_PAR.Hnormal    = Hnormal;                      % height profile
SASA_PAR.dist_vec   = dist_vec;                     % distance(s) in \mum

%% perform SASA
addpath(PATH)
SMAT_SASA = SASA_for2(MATFILE_,ANG_,MIRROR_,SASA_PAR,SAVE_PATH);
SMAT_SASA = squeeze(SMAT_SASA);
rmpath(PATH)

%% imagesc showing s-matrix components over distance and frequency

addpath ../../calc/utilities

lambda_vec = SASA_PAR.lambda_vec;
FREQ_VEC   = 3e8./(lambda_vec/1e6)/(1e12);

% set coefficients for the s-matrix subplot
COEF = cell(4,4);
COEF{1,1} = '\tau^f_{xx}'; COEF{1,2} = '\tau^f_{xy}'; COEF{1,3} = '\rho^b_{xx}'; COEF{1,4} = '\rho^b_{xy}';
COEF{2,1} = '\tau^f_{yx}'; COEF{2,2} = '\tau^f_{yy}'; COEF{2,3} = '\rho^b_{yx}'; COEF{2,4} = '\rho^b_{yy}';
COEF{3,1} = '\rho^f_{xx}'; COEF{3,2} = '\rho^f_{xy}'; COEF{3,3} = '\tau^b_{xx}'; COEF{3,4} = '\tau^b_{xy}';
COEF{4,1} = '\rho^f_{yx}'; COEF{4,2} = '\rho^f_{yy}'; COEF{4,3} = '\tau^b_{yx}'; COEF{4,4} = '\tau^b_{yy}';

sp_vec = [1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16];
S_IND = [1 6 3 8];

% Determine the factor Re(ksz)/kcz (definition of clad and subs is 
% reveresed in FOMO script).
n_fac = SASA_PAR.n_subs/SASA_PAR.n_clad;
% multiply trans_fac to transmission blocks (for constant ref. indices)
ABS_SMAT_SASA              = abs( SMAT_SASA ).^2;
ABS_SMAT_SASA(:,:,1:2,1:2) = n_fac * ABS_SMAT_SASA(:,:,1:2,1:2);
ABS_SMAT_SASA(:,:,3:4,3:4) = n_fac * ABS_SMAT_SASA(:,:,3:4,3:4);
% 
% % throw out singularity
FREQ_VEC(12)            = [];
ABS_SMAT_SASA(:,12,:,:) = [];

IN_DATA.x       = FREQ_VEC;
IN_DATA.y       = dist_vec*1000;
IN_DATA.smat    = ABS_SMAT_SASA;
LabelStr.x      = 'frequency in THz';
LabelStr.y      = 'layer distance in nm';
LabelStr.CB     = COEF;
LabelStr.fontsz = 12;

figure(figh1);clf
omni_imagesc_fullSMAT(IN_DATA,LabelStr)

