% perform a geometric expansion (Bremmer series) of the SASA stack
% calculated in 'match_square_wire_stack_final.m'

match_square_wire_stack_final

close all
%% run geometric SASA

% expand up to order GEO_ORDER and save complete (trunkated) series in cell
GEO_ORDER = 2;
SMAT_SASA_GEO_CELL = cell(1,GEO_ORDER + 1);
for loop = 0:GEO_ORDER
    
    tic
    SMAT_SASA_GEO_CELL{loop + 1} = SASA_GEO(SASA_PARS, loop);
    toc
    
end

disp('next')

% substract series of order (GEO_ORDER - 1) from the one of order GEO_ORDER
% to get only the GEO_ORDERth term
HIGHER_ORDER_CELL = cell(1,GEO_ORDER);
for loop = 1:GEO_ORDER
    tic
    
    HIGHER_ORDER_CELL{loop} = ...
        SMAT_SASA_GEO_CELL{loop + 1} - SMAT_SASA_GEO_CELL{loop};
    
    toc
end

%% amplitude and phases of each order coefficient comparing to full series

% plot parameters
fighand = 24;

NX      = 1;
NY      = 1;

stepsz = 12;

pic_points = 1:stepsz:length(lambda_FMM);

fsz = 17;
LabelStr.x = 'wavelength in \mum';
LabelStr.y = 'amplitude';

% plot all in one (0th to nth)
figure(fighand)
clf
% ======================= Amplitude subplot ===============================
% subplot(1,2,1)

% zeroth order
ph1 = plot(lambda_FMM, ...
    squeeze(abs(SMAT_SASA_GEO_CELL{1}(:,NX,NY))));
omni_plot(ph1,fsz,LabelStr,true)

hold on

% each higher order
for loop = 1:GEO_ORDER
    ph1 = plot(lambda_FMM, ...
        squeeze(abs(HIGHER_ORDER_CELL{loop}(:,NX,NY))));
    omni_plot(ph1,fsz,LabelStr,true)
end

% trunkated geometric series
ph1 = plot(lambda_FMM(pic_points), ...
    abs( squeeze( SMAT_SASA_GEO_CELL{end}(pic_points,NX,NY) ) ), 'o');
omni_plot(ph1,fsz,LabelStr,true)

% full SASA
ph1 = plot(lambda_FMM, abs( squeeze( SMAT_SASA(:,NX,NY) ) ));
omni_plot(ph1,fsz,LabelStr)

% set limits that fit ALL :)
ylim([min( squeeze( abs( HIGHER_ORDER_CELL{loop}(:,NX,NY) ) ) ) - 0.01 ...
    max( abs( squeeze( SMAT_SASA(:,NX,NY) ) ) )] + 0.01)

hold off

legend('leading order','NLO','NNLO','truncated series','full stack', ...
    'location', 'east')

% ========================= Phase subplot =================================
% subplot(1,2,2)

figure(fighand + 1)
LabelStr.y = 'phase';

% zeroth order
ph2 = plot(lambda_FMM, ...
    unwrap(squeeze(angle(SMAT_SASA_GEO_CELL{1}(:,NX,NY)))));
omni_plot(ph2,fsz,LabelStr,true);

hold on

% each higher order
for loop = 1:GEO_ORDER
    ph2 = plot(...
        lambda_FMM, ...
        unwrap( squeeze( angle( HIGHER_ORDER_CELL{loop}(:,NX,NY) ) ) )...
        );
    omni_plot(ph2,fsz,LabelStr)
end

% trunkated geometric series
ph2 = plot(lambda_FMM(pic_points), unwrap( angle( squeeze( ...
    SMAT_SASA_GEO_CELL{end}(pic_points,NX,NY) ...
    ) ) ), 'o');
omni_plot(ph2,fsz,LabelStr,true)

% full SASA
ph2 = plot(lambda_FMM, unwrap( angle( squeeze( SMAT_SASA(:,NX,NY) ) ) ) );
omni_plot(ph2,fsz,LabelStr,true)

hold off

legend('leading order','NLO','NNLO','truncated series','full stack', ...
    'location', 'best')

%% amplitude and phase of the series truncated at different points
% plot parameters
fighand = 34;

NX      = 1;
NY      = 1;

stepsz = 12;

pic_points = 1:stepsz:length(lambda_FMM);

fsz = 17;
LabelStr.x = 'wavelength in \mum';
LabelStr.y = 'amplitude';

% plot all in one (0th to nth)
figure(fighand)
clf
% ======================= Amplitude subplot ===============================

hold on

% trunkated geometric series at different order
for loop = 1:GEO_ORDER + 1
    ph1 = plot(lambda_FMM, ...
        abs( squeeze( SMAT_SASA_GEO_CELL{loop}(:,NX,NY) ) ));
    omni_plot(ph1,fsz,LabelStr,true)
end

% full SASA
ph1 = plot(lambda_FMM, abs( squeeze( SMAT_SASA(:,NX,NY) ) ));
omni_plot(ph1,fsz,LabelStr)

% set limits that fit ALL :)
% ylim([min( squeeze( abs( HIGHER_ORDER_CELL{loop}(:,NX,NY) ) ) ) - 0.01 ...
%     max( abs( squeeze( SMAT_SASA(:,NX,NY) ) ) )] + 0.01)

hold off

legend('leading order','NLO','NNLO','full stack', ...
    'location', 'east')

% ========================= Phase subplot =================================
figure(fighand + 1)
LabelStr.y = 'phase';
clf
hold on

% trunkated geometric series at different orders
for loop = 1:GEO_ORDER + 1
    ph2 = plot(...
        lambda_FMM, ...
        unwrap( angle( squeeze( ...
    SMAT_SASA_GEO_CELL{loop}(:,NX,NY) ...
                ) ...
            ) ...
        ) ...
    );
    omni_plot(ph2,fsz,LabelStr)
end

% full SASA
ph2 = plot(lambda_FMM, unwrap( angle( squeeze( SMAT_SASA(:,NX,NY) ) ) ) );
omni_plot(ph2,fsz,LabelStr,true)

hold off

legend('leading order','NLO','NNLO','full stack', ...
    'location', 'best')

%% save plot data for python processing

save('SASA_geo_data.mat', 'SMAT_SASA', 'SMAT_SASA_GEO_CELL', ...
    'HIGHER_ORDER_CELL', 'lambda_FMM')

