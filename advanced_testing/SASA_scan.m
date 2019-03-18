% match square-wire stack
% SASA scan over layer distance

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load FMM data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_name_1 = 'avg_square_final';
data_name_2 = 'avg_wire_final';

addpath(data_name_1, data_name_2)
addpath ../calc/materials
addpath ../calc/SASA_source/
addpath ../calc/utilities

load([data_name_1, '_Daten_gesamt.mat'])
SMAT_Wire = SMAT_;

load([data_name_2,'_Daten_gesamt.mat'])
SMAT_Square = SMAT_;

% parameter space:

% --- squares ---
% Groove height in nm
Height_squares = [55 58 60] / 1000; % metal layer thickness in nm

% --- wires ---
% glass data: 1 - Sellmeier SiO2, 2 - Futurrex
Glass_data = [1 2];
% gold data: 1 - Johnson, 2 - Kay, 3 - olmon
Au_data = [1 2 3];
% averaging threshold: all, largest, smallest
Thresh_ = [0 1 9];
% Groove height in nm
Height_wires = [35 41 45] / 1000; % metal layer thickness in nm

% wavelength in \mum
lambda_FMM    = linspace(0.47, 1.2, 512);

%% refractive index with dispersion (fused silica)

% FMM-based data
n_Sili = zeros(1, length(lambda_FMM));
n_Futu = zeros(1, length(lambda_FMM));

for loop = 1:length(lambda_FMM)
    
    % Futurrex index
    n_Futu(loop) = refractive_ind(lambda_FMM(loop) * 1000, 10, 2);
    
    % fused silica index (Sellmeier)
    n_Sili(loop) = refractive_ind(lambda_FMM(loop), 2, 1);
    
end

% vacuum or air index
n_vac = ones(1, length(lambda_FMM));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SASA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stack from right to left, yo!!!

% spacer height for scan
h_spacer      = (400:1000) / 1000;

% SASA parameters (distances in \mum)
h_cladding    = 585 / 1000; % cladding thickness

meta_layer_1 = struct(...
    'Smatrix', SMAT_Square, ...
    'embedding', {n_Futu, n_Sili} ...
    );

meta_layer_2 = struct(...
    'Smatrix', SMAT_Wire, ...
    'embedding', {n_Sili, n_Futu} ...
    );

stack_medium = {meta_layer_2, n_Futu, meta_layer_1, n_Sili};

stack_height = {0, h_spacer, 0, h_cladding};


n_embed      = struct('cladding', n_Sili, 'substrate', n_vac);

SASA_PARS = struct(...
    'Lambda', lambda_FMM,...
    'medium', stack_medium,...
    'height', stack_height,...
    'embedding', n_embed);

% run SASA
tic
SMAT_SASA = SASA(SASA_PARS);
toc

SMAT_SASA_INT = abs(SMAT_SASA).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot transmittance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_Sili_mat = repmat(n_Sili, [length(h_spacer) 1]);

SASA_INT_X = squeeze(SMAT_SASA_INT(:,:,1,1)) .* real(n_Sili_mat);
SASA_INT_Y = squeeze(SMAT_SASA_INT(:,:,2,2)) .* real(n_Sili_mat);

figure(110)
imagesc(lambda_FMM * 1000, h_spacer * 1000, SASA_INT_X)
colorbar

figure(111)
imagesc(lambda_FMM * 1000, h_spacer * 1000, SASA_INT_Y)
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot transmittance with frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_Sili_mat = repmat(n_Sili, [length(h_spacer) 1]);

SASA_INT_X = squeeze(SMAT_SASA_INT(:,:,1,1)) .* real(n_Sili_mat);
SASA_INT_Y = squeeze(SMAT_SASA_INT(:,:,2,2)) .* real(n_Sili_mat);

lambda_FMM_nonlin = 1./linspace(1/1.2,1/0.47,512);
freq_vec = 300./lambda_FMM_nonlin;

SASA_INT_X_freq = SASA_INT_X;
SASA_INT_Y_freq = SASA_INT_Y;

for loop1 = 1:length(h_spacer)
    SASA_INT_X_freq(loop1,:) = ...
        interp1(lambda_FMM, squeeze(SASA_INT_X(loop1, :)), lambda_FMM_nonlin);
    SASA_INT_Y_freq(loop1,:) = ...
        interp1(lambda_FMM, squeeze(SASA_INT_Y(loop1, :)), lambda_FMM_nonlin);
end

figure(112)
imagesc(freq_vec, h_spacer * 1000, SASA_INT_X_freq)
colorbar

figure(113)
imagesc(freq_vec, h_spacer * 1000, SASA_INT_Y_freq)
colorbar


