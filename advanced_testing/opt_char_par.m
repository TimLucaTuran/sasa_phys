function PAR_OUT = opt_char_par(SMAT_IN)

% give CH, OA, and other characterization parameters

% transform to circular basis
PATH1 ='../calc/utilities';
addpath(PATH1)

SMAT_C = trafo_cart2circ(SMAT_IN);

% denote s-amtrix elements

DIM_IN = ndims(SMAT_IN);
II = cell(4,4);
for loop1 = 1:4
    for loop2 = 1:4
        II{loop1,loop2} = smat_phys_ind_arb_sz(DIM_IN,loop1,loop2);
    end
end

% front transmission
TFXX = SMAT_IN( II{1,1}{:} );
TFXY = SMAT_IN( II{1,2}{:} );
TFYX = SMAT_IN( II{2,1}{:} );
TFYY = SMAT_IN( II{2,2}{:} );

% back transmissiom
TBXX = SMAT_IN( II{3,3}{:} );
TBXY = SMAT_IN( II{3,4}{:} );
TBYX = SMAT_IN( II{4,3}{:} );
TBYY = SMAT_IN( II{4,4}{:} );

% back reflection
RBXX = SMAT_IN( II{1,3}{:} );
RBXY = SMAT_IN( II{1,4}{:} );
RBYX = SMAT_IN( II{2,3}{:} );
RBYY = SMAT_IN( II{2,4}{:} );

% front reflection
RFXX = SMAT_IN( II{3,1}{:} );
RFXY = SMAT_IN( II{3,2}{:} );
RFYX = SMAT_IN( II{4,1}{:} );
RFYY = SMAT_IN( II{4,2}{:} );

% circular base
% T-front block
TFRR = SMAT_C( II{1,1}{:} );
TFRL = SMAT_C( II{1,2}{:} );
TFLR = SMAT_C( II{2,1}{:} );
TFLL = SMAT_C( II{2,2}{:} );

% T-back block
TBRR = SMAT_C( II{3,3}{:} );
TBRL = SMAT_C( II{3,4}{:} );
TBLR = SMAT_C( II{4,3}{:} );
TBLL = SMAT_C( II{4,4}{:} );

% R-back block
RBRR = SMAT_C( II{1,3}{:} );
RBRL = SMAT_C( II{1,4}{:} );
RBLR = SMAT_C( II{2,3}{:} );
RBLL = SMAT_C( II{2,4}{:} );

% R-front block
RFRR = SMAT_C( II{3,1}{:} );
RFRL = SMAT_C( II{3,2}{:} );
RFLR = SMAT_C( II{4,1}{:} );
RFLL = SMAT_C( II{4,2}{:} );

%% calculte characteristics for TRANSMISSION
% reciprocal media:

% asymmetric transmission
% asymmetric transmission for x-pol. (TM) 
DELTA_Y = abs(TFYX).^2 - abs(TFXY).^2;

% asymmetric transmission for y-pol. (TE)
DELTA_X = -1*DELTA_Y;

% asymmetric transmission for R-pol. (TM) 
DELTA_L = abs(TFLR).^2 - abs(TFRL).^2;

% asymmetric transmission for L-pol. (TE)
DELTA_R = -1*DELTA_L;

%% stokes parameters
% Stokes parameters front

% x-polarized state
PHIDIFF_x = angle(TFXX) - angle(TFYX); % phase differnce
S0_F_x = abs(TFXX).^2 + abs(TFYX).^2;
S1_F_x = abs(TFXX).^2 - abs(TFYX).^2;
S2_F_x = 2*abs(TFXX).*abs(TFYX).*cos(PHIDIFF_x);
S3_F_x = 2*abs(TFXX).*abs(TFYX).*sin(PHIDIFF_x);

% y-polarized state
PHIDIFF_y = angle(TFXY) - angle(TFYY); % phase differnce
S0_F_y = abs(TFYY).^2 + abs(TFXY).^2;
S1_F_y = abs(TFYY).^2 - abs(TFXY).^2;
S2_F_y = 2*abs(TFYY).*abs(TFXY).*cos(PHIDIFF_y);
S3_F_y = 2*abs(TFYY).*abs(TFXY).*sin(PHIDIFF_y);

%stokes parameters back

% x-polarized state
PHIDIFB_x = angle(TBXX) - angle(TBYX); % phase differnce
S0_B_x = abs(TBXX).^2 + abs(TBYX).^2;
S1_B_x = abs(TBXX).^2 - abs(TBYX).^2;
S2_B_x = 2*abs(TBXX).*abs(TBYX).*cos(PHIDIFB_x);
S3_B_x = 2*abs(TBXX).*abs(TBYX).*sin(PHIDIFB_x);

% y-polarized state
PHIDIFB_y = angle(TBXY) - angle(TBYY); % phase differnce
S0_B_y = abs(TBYY).^2 + abs(TBXY).^2;
S1_B_y = abs(TBYY).^2 - abs(TBXY).^2;
S2_B_y = 2*abs(TBYY).*abs(TBXY).*cos(PHIDIFB_y);
S3_B_y = 2*abs(TBYY).*abs(TBXY).*sin(PHIDIFB_y);

%% ellipticity

% front
% x-state
ELLI_F_x = S3_F_x./S0_F_x;

% y-state
ELLI_F_y = S3_F_y./S0_F_y;

% back
% x-state
ELLI_B_x = S3_B_x./S0_B_x;

% y-state
ELLI_B_y = S3_B_y./S0_B_y;

%% orientation of the polarization ellipse, its rotation angle

% front
% x-state
PSI_F_x = atan( S2_F_x ./ S0_F_x ) / 2;

% y-state
PSI_F_y = atan( S2_F_y ./ S0_F_y ) / 2;

% back
% x-state
PSI_B_x = atan( S2_B_x ./ S0_B_x ) / 2;

% y-state
PSI_B_y = atan( S2_B_y ./ S0_B_y ) / 2;

%% other parameters
% chiral dichroism
THETA_ = atan( (TFLR - TFRL)./(TFRR + TFLR + TFLL + TFRL) ); 

% circular birefringence front and back
ZETA_CIRC = angle(TFRR./TFLL);
ZETA_CIRCB = angle(TBRR./TBLL);

% optical rotation front and back
PHI_ = angle(TFRR + TFLR) - angle(TFLL + TFRL);
PHI_B = angle(TBRR + TBLR) - angle(TBLL + TBRL);

%% put everything 2gether into struct for transmission

PAR_ASYM = struct( ...
    'DeltaX', DELTA_X, ...
    'DeltaY', DELTA_Y, ...
    'DeltaR', DELTA_R, ...
    'DeltaL', DELTA_L);

PAR_CIRC_BI = struct( ...
    'front', ZETA_CIRC,...
    'back',  ZETA_CIRCB);

PAR_OPTROT = struct( ...
    'front',     PHI_, ...
    'back',      PHI_B, ...
    'DiffFront', PHIDIFF_y, ...
    'DiffBack',  PHIDIFB_y);

PAR_STOKES_F_x = struct( ...
    'S0', S0_F_x,...
    'S1', S1_F_x,...
    'S2', S2_F_x,...
    'S3', S3_F_x);

PAR_STOKES_F_y = struct( ...
    'S0', S0_F_y,...
    'S1', S1_F_y,...
    'S2', S2_F_y,...
    'S3', S3_F_y);

PAR_STOKES_F = struct('x', PAR_STOKES_F_x, 'y', PAR_STOKES_F_y);

PAR_STOKES_B_x = struct( ...
    'S0', S0_B_x,...
    'S1', S1_B_x,...
    'S2', S2_B_x,...
    'S3', S3_B_x);

PAR_STOKES_B_y = struct( ...
    'S0', S0_B_y,...
    'S1', S1_B_y,...
    'S2', S2_B_y,...
    'S3', S3_B_y);

PAR_STOKES_B = struct('x', PAR_STOKES_B_x, 'y', PAR_STOKES_B_y);

PAR_STOKES = struct('front', PAR_STOKES_F, 'back', PAR_STOKES_B);

PAR_ELLI_F = struct('x', ELLI_F_x, 'y', ELLI_F_y);
PAR_ELLI_B = struct('x', ELLI_B_x, 'y', ELLI_B_y);

PAR_ELLI = struct( ...
    'front', PAR_ELLI_F, ...
    'back',  PAR_ELLI_B);

PAR_ROTANG_F = struct('x', PSI_F_x, 'y', PSI_F_y);
PAR_ROTANG_B = struct('x', PSI_B_x, 'y', PSI_B_y);

PAR_ROTANG = struct(...
    'front', PAR_ROTANG_F,...
    'back', PAR_ROTANG_B);

PAR_TRANS = struct(...
    'asymmetric',        PAR_ASYM, ...
    'ChiralDichroism',   THETA_, ...
    'CircBirefringence', PAR_CIRC_BI, ...
    'OpticalRot',        PAR_OPTROT, ...
    'stokes',            PAR_STOKES, ...
    'ellipticity',       PAR_ELLI, ...
    'RotationAngle',     PAR_ROTANG);

%% calculte characteristics for REFLECTION
% asymmetric transmission for x-pol. (TM) 
DELTAX_R = abs(RFXX).^2 + abs(RFYX).^2 - abs(RBXX).^2 - abs(RBYX).^2;

% asymmetric transmission for y-pol. (TE)
DELTAY_R = abs(RFYY).^2 + abs(RFXY).^2 - abs(RBYY).^2 - abs(RBXY).^2;

% asymmetric transmission for R-pol. (TM) 
DELTAR_R = abs(RFRR).^2 + abs(RFRL).^2 - abs(RBRR).^2 - abs(RBRL).^2;

% asymmetric transmission for L-pol. (TE)
DELTAL_R = abs(RFLL).^2 + abs(RFLR).^2 - abs(RBLL).^2 - abs(RBLR).^2;

% chiral dichroism
THETA_R = atan( (RFRR + RFLR - RFLL - RFRL)./(RFRR + RFLR + RFLL + RFRL) ); 

% circular birefringence front and back
ZETA_CIRC_R = angle(RFRR./RFLL); 
ZETA_CIRCB_R = angle(RBRR./RBLL);

% optical rotation front and back
PHI_R = angle(RFRR + RFLR) - angle(RFLL + RFRL);
PHI_B_R = angle(RBRR + RBLR) - angle(RBLL + RBRL);

% Stokes parameters front
PHIDIFF_R = angle(RFYY) - angle(RFXY); % phase differnce
S0_F_R = abs(RFYY).^2 + abs(RFXY).^2;
S1_F_R = abs(RFYY).^2 - abs(RFXY).^2;
S2_F_R = 2*abs(RFYY).*abs(RFXY).*cos(PHIDIFF_R);
S3_F_R = 2*abs(RFYY).*abs(RFXY).*sin(PHIDIFF_R);

%stokes parameters back
PHIDIFB_R = angle(RBYY) - angle(RBXY); % phase differnce
S0_B_R = abs(RBYY).^2 + abs(RBXY).^2;
S1_B_R = abs(RBYY).^2 - abs(RBXY).^2;
S2_B_R = 2*abs(RBYY).*abs(RBXY).*cos(PHIDIFB_R);
S3_B_R = 2*abs(RBYY).*abs(RBXY).*sin(PHIDIFB_R);

% normalized ellipticity front
CHI_F_R = S3_F_R./S0_F_R;

% rotation angle for TE and TM front
CHIF_TE_R = atan(real(RFYX./RFXX));
CHIF_TM_R = atan(real(RFXY./RFYY));

% normalized ellipticity back
CHI_B_R = S3_B_R./S0_B_R;

% rotation angle for TE and TM back
CHIB_TE_R = atan(real(RBYX./RBXX));
CHIB_TM_R = atan(real(RBXY./RBYY));

%% put everything 2gether into struct for reflection

warning('Parameters in backward direction have yet to be updated...')

PAR_ASYM_R = struct( ...
    'DeltaX', DELTAX_R, ...
    'DeltaY', DELTAY_R, ...
    'DeltaR', DELTAR_R, ...
    'DeltaL', DELTAL_R);

PAR_CIRC_BI_R = struct( ...
    'front', ZETA_CIRC_R,...
    'back',  ZETA_CIRCB_R);

PAR_OPTROT_R = struct( ...
    'front',     PHI_R, ...
    'back',      PHI_B_R, ...
    'DiffFront', PHIDIFF_R, ...
    'DiffBack',  PHIDIFB_R);
  
PAR_STOKES_R_F = struct( ...
    'S0', S0_F_R,...
    'S1', S1_F_R,...
    'S2', S2_F_R,...
    'S3', S3_F_R);

PAR_STOKES_R_B = struct( ...
    'S0', S0_B_R,...
    'S1', S1_B_R,...
    'S2', S2_B_R,...
    'S3', S3_B_R);

PAR_STOKES_R = struct('front', PAR_STOKES_R_F, 'back', PAR_STOKES_R_B);

PAR_ELLI_R = struct(...
    'front', CHI_F_R,...
    'back',  CHI_B_R);

PAR_ROTANG_R = struct(...
    'frontTE', CHIF_TE_R,...
    'frontTM', CHIF_TM_R,...
    'backTE',  CHIB_TE_R,...
    'backTM',  CHIB_TM_R);

PAR_REFL = struct(...
    'asymmetric',        PAR_ASYM_R, ...
    'ChiralDichroism',   THETA_R, ...
    'CircBirefringence', PAR_CIRC_BI_R, ...
    'OpticalRot',        PAR_OPTROT_R, ...
    'stokes',            PAR_STOKES_R, ...
    'ellipticity',       PAR_ELLI_R, ...
    'RotationAngle',     PAR_ROTANG_R);

%% create final output struct
PAR_OUT = struct('transmission', PAR_TRANS, 'reflection', PAR_REFL);

