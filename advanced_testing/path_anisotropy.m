% evaluate the level of anisotropy of each Feynman path
% 
% run square_wire_stack_GEO.m before using this script
% OR load 

% in the following, we interpret the Feynman paths as metasurfaces in their
% own right and analyze their symmetry

%% anisotropy at 0th order
% get characteristic optical parameter using opt_char_par.m
opt_char_pars_0 = opt_char_par(SMAT_SASA_GEO_CELL{1});

%
elli_0_x = opt_char_pars_0.transmission.ellipticity.front.x;
elli_0_y = opt_char_pars_0.transmission.ellipticity.front.y;

figure(10)
plot(lambda_FMM, elli_0_x, lambda_FMM, elli_0_y)
ylim([-0.4, 0.6])

%% anisotropy at 1st order
% get characteristic optical parameter using opt_char_par.m
opt_char_pars_1 = opt_char_par(HIGHER_ORDER_CELL{1});

%% 
elli_1_x = opt_char_pars_1.transmission.ellipticity.front.x;
elli_1_y = opt_char_pars_1.transmission.ellipticity.front.y;

figure(11)
plot(lambda_FMM, elli_1_x, lambda_FMM, elli_1_y)

%% anisotropy at 2nd order
% get characteristic optical parameter using opt_char_par.m
opt_char_pars_2 = opt_char_par(HIGHER_ORDER_CELL{2});

%% 
elli_2_x = opt_char_pars_2.transmission.ellipticity.front.x;
elli_2_y = opt_char_pars_2.transmission.ellipticity.front.y;

figure(12)
plot(lambda_FMM, elli_2_x, lambda_FMM, elli_2_y)

%% anisoptropy for A LOT of orders: see how it evolves
clear HIGHER_ORDER_CELL SMAT_SASA_GEO_CELL
load higher_orders.mat

opt_char_par_full = opt_char_par(SMAT_SASA);
elli_x_full = opt_char_par_full.transmission.ellipticity.front.x;
elli_y_full = opt_char_par_full.transmission.ellipticity.front.y;

figure(13)
hold on

elli_arrx = zeros(length(lambda_FMM), length(SMAT_SASA_GEO_CELL));
elli_arry = zeros(length(lambda_FMM), length(SMAT_SASA_GEO_CELL));

for loop = 1:length(SMAT_SASA_GEO_CELL)
    opt_char_pars_HIGH = opt_char_par(SMAT_SASA_GEO_CELL{loop});
    elli_high_x = opt_char_pars_HIGH.transmission.ellipticity.front.x;
    elli_high_y = opt_char_pars_HIGH.transmission.ellipticity.front.y;
    
    plt = plot(lambda_FMM, elli_high_x, 'b',...
    lambda_FMM, elli_high_y, 'g');
    
    elli_diff_x(:,loop) = abs(elli_high_x - elli_x_full);
    elli_diff_y(:,loop) = abs(elli_high_y - elli_y_full);
    
    elli_arrx(:,loop) = elli_high_x;
    elli_arry(:,loop) = elli_high_y;
end

hold off

%%

figure(15)
plot(lambda_FMM, elli_arry(:,1), lambda_FMM, elli_0_y, lambda_FMM, elli_y_full)

%% compare to full SASA results

figure(14)
imagesc(lambda_FMM, 1:length(SMAT_SASA_GEO_CELL), log10(elli_diff_y).')
colorbar

%% save ellipticity data (workspace) for processing in python

save('SASA_GEO_elli_data.mat')

