% match square-wire stack
% Compares CES-measurment of the stack with SASA-results
clear

addpath SFA030b
addpath ../advanced_testiing/materials
load('SFA030b.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load FMM data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_name_1 = 'avg_square_final';
data_name_2 = 'avg_wire_final';

addpath(data_name_1)
addpath(data_name_2)
addpath ../advanced_testing/materials

load([data_name_2, '_Daten_gesamt.mat'])
SMAT_Wire = SMAT_;

load([data_name_1,'_Daten_gesamt.mat'])
SMAT_Square = SMAT_;

% wavelength in \mum
lambda_FMM    = linspace(0.47, 1.2, 512);

%% refractive index with dispersion (fused silica)

% FMM-based data
n_Futu = zeros(1, length(lambda_FMM));
n_Sili = zeros(1, length(lambda_FMM));

for loop = 1:length(lambda_FMM)
    
    % Futurrex index
    n_Futu(loop) = real(refractive_ind(lambda_FMM(loop) * 1000, 10, 2));
    
    % fused silica index
    n_Sili(loop) = refractive_ind(lambda_FMM(loop), 2, 1);

end

% vacuum or air index
n_vac = ones(1, length(lambda_FMM));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SASA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stack from right to left, yo!!!

addpath ../calc/SASA_source/
addpath ../calc/utilities

% SASA parameters (distances in \mum)

h_spacer      = 450 / 1000; % spacer height
h_cladding    = 585 / 1000; % cladding thickness

meta_layer_u = struct(...
    'Smatrix', SMAT_Square, ...
    'embedding', {n_Futu, n_Sili} ...
    );

meta_layer_l = struct(...
    'Smatrix', SMAT_Wire, ...
    'embedding', {n_Sili, n_Futu} ...
    );

stack_medium = {meta_layer_l, n_Futu, meta_layer_u, n_Sili};

stack_height = {0, h_spacer, 0, h_cladding};

n_embed      = struct('cladding', n_Sili, 'substrate', n_vac);

SASA_PARS = struct(...
    'Lambda', lambda_FMM,...
    'medium', stack_medium,...
    'height', stack_height,...
    'embedding', n_embed);

% run SASA
SMAT_SASA = SASA(SASA_PARS);

SMAT_SASA_INT = abs(SMAT_SASA).^2;
SMAT_SASA_ARG = angle(SMAT_SASA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load CES-data of squares %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dosis and field choice
Names = {'D100','D110','D121','D133'};

lambda_CES_amp = Samples.([Names{1}, 'I']).Amp.wl * 1e6;
lambda_CES_arg = Samples.([Names{1}, 'I']).Arg.wl * 1e6;

% initialize
INT_CES_X  = cell(1, length(Names));
INT_CES_Y  = cell(1, length(Names));
ARG_CES_X  = cell(1, length(Names));
ARG_CES_Y  = cell(1, length(Names));

CES_amp_correct = 0.96;

% glass height which adds phase (to substract)
h_wire  = 45/1000;
h_patch = (55 + 4)/1000;

some_height_1 = h_wire + h_spacer;
some_height_2 = h_patch + h_cladding;

n_Futu_CES = real(refractive_ind(lambda_CES_arg * 1000,10,2));
n_Sili_CES = refractive_ind(lambda_CES_arg, 2, 1);

CES_phase_correct_1 = zeros(length(lambda_CES_arg),1);
for loop = 1:length(lambda_CES_arg)
    
    CES_phase_correct_1(loop) = n_Futu_CES(loop) * 2*pi * ...
        some_height_1 / lambda_CES_arg(loop);
end

CES_phase_correct_2 = zeros(length(lambda_CES_arg),1);
for loop = 1:length(lambda_CES_arg)
    
    CES_phase_correct_2(loop) = n_Sili_CES(loop) * 2*pi * ...
        some_height_2 / lambda_CES_arg(loop);
end

CES_phase_correct = CES_phase_correct_1 + CES_phase_correct_2;

for name_loop = 1:length(Names)

    % extract intensity (transmittance)
    INT_CES_X{name_loop} = Samples.([Names{name_loop}, 'I']).Amp.PX * CES_amp_correct;
    INT_CES_Y{name_loop} = Samples.([Names{name_loop}, 'I']).Amp.PY * CES_amp_correct;
    
    % extract argument (phase)
    ARG_CES_X{name_loop} = Samples.([Names{name_loop}, 'I']).Arg.PX + CES_phase_correct;
    ARG_CES_Y{name_loop} = Samples.([Names{name_loop}, 'I']).Arg.PY + CES_phase_correct;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot transmittance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pick a field of the waver
field_N = 3;

figure(506)
SASA_INT_X = squeeze(SMAT_SASA_INT(:,1,1)) .* real(n_Sili).';
SASA_INT_Y = squeeze(SMAT_SASA_INT(:,2,2)) .* real(n_Sili).';

p_hand = plot( lambda_CES_amp * 1000, INT_CES_X{field_N},'-b',...
    lambda_CES_amp * 1000, INT_CES_Y{field_N},':b',...
    lambda_FMM * 1000, SASA_INT_X, 'r--',...
    lambda_FMM * 1000, SASA_INT_Y, 'r');

LabelStr.x = 'wavelength in nm';
LabelStr.y = 'transmittance';
%omni_plot(p_hand,17,LabelStr)
title('Square-wire stack')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pick a field of the waver
field_N = 3;

figure(507)
SASA_ARG_X = unwrap(squeeze(SMAT_SASA_ARG(:,1,1)));
SASA_ARG_Y = unwrap(squeeze(SMAT_SASA_ARG(:,2,2)));

p_hand = plot( lambda_CES_arg * 1000, ARG_CES_X{field_N} - 2*pi,'-b',...
    lambda_CES_arg * 1000, ARG_CES_Y{field_N} - 2*pi,':b',...
    lambda_FMM * 1000, SASA_ARG_X + 6*pi, 'r--',...
    lambda_FMM * 1000, SASA_ARG_Y + 6*pi, 'r');

LabelStr.x = 'wavelength in nm';
LabelStr.y = 'transmittance';
%omni_plot(p_hand,17,LabelStr)
title('Square-wire stack')

save('match_square_wire_stack_final.mat')
