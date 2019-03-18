function INDEX = refractive_ind(LAMBDA,MATERIAL,FIT_INTERP)
% Usage: INDEX = refractive_ind(LAMBDA,MATERIAL,FIT_INTERP)
%
% Find refractive index INDEX for a given wavelength LAMBDA.
% Choose materials with MATERIALS: 
% 1 - saphire, 
% 2 - fused silica, 
% 3 - silicon,
% 4 - titania, 
% 5 - P205 (VIS), 
% 6 - P205 (amorphous, IR), 
% 7 - Si3N4, 
% 8 - MgF2 (amorphous),
% 9 - futurex (by Dennis, measured for ~200nm),
% 10 - futurrex (by FUTURREX)
% 11 - Ti from Johnson, Christy Data
% 12 - Ti from Rakic et. al. (1998) Drude-Lorentz-model
% 13 - Ti from Rakic et. al. (1998) Brendel-Bormann-model
%
% Choose with FIT_INTERP, whether to use a fit function (1; e.g. Sellmeier eq.)
% or interplolate(2; cubic). LAMBDA has the unit \mum
% The usual threshold is between 0.5 and 1.5 \mum (a bit lower or higher
% shouldn't hurt).
% 
% NOTE: only 1, 2, and 7 have fit functions available. The rest of the data
% (including the three mentioned) are based on measured data.

x = LAMBDA;
switch(FIT_INTERP)
    
    case 1 % use fit function (Sellmeier eq.)
        switch(MATERIAL)
            
            case 1 % saphire index (from refractiveindex.info)
                INDEX = sqrt(1+1.4313493./(1-power(0.0726631./x,2)) + ...
                    0.65054713./(1-power(0.1193242./x,2))...
                    + 5.3414021./(1-power(18.028251./x,2)));
                
            case 2 % silica index (from refractiveindex.info)
                INDEX = sqrt(1+0.6961663./(1-power(0.0684043./x,2)) + ...
                    0.4079426./(1-power(0.1162414./x,2))...
                    + 0.8974794./(1-power(9.896161./x,2)));
                
            case 3 % silicon
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
                
            case 4 % TiO2
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
                
            case 5 % P2O5
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
                
            case 6 % P2O5 (amorphous)
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
                
            case 7 % Si3N4 from 0.310 to 5.504 \mum (from refractiveindex.info)
                INDEX = sqrt(1 + 3.0249./(1 - (0.1353406./x).^2) + ...
                    40314./(1 - (1239.842./x).^2));
                
            case 8 % amorphous MgF2 measured 2012 by ??? (V700)
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
            case 9 % futurex (by Dennis, measured for ~200nm)
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
            case 10 % futurrex (by FUTURREX)
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
            case 11 % Ti from Johnson, Christy Data
                msg = ['No fit function available for material ',...
                    num2str(MATERIAL)];
                error(msg);
                     
            case 12 % Ti from Rakic et. al. (1998) Drude-Lorentz-model
                EPS_Ti = Ti_LD_Rakic(LAMBDA);
                INDEX  = sqrt(EPS_Ti);
                
            case 13 % Ti from Rakic et. al. (1998) Brendel-Bormann-model
                EPS_Ti = Ti_BB_Rakic(LAMBDA);
                INDEX  = sqrt(EPS_Ti);
                
        end
        
    case 2 % use Piecewise Cubic Hermite Interpolating Polynomial
        % the variable IN_FILE is loaded with the material mat-files
        switch(MATERIAL)
            
            case 1 % saphire index (from refractiveindex.info)
                load('Malitson-o.mat')
                INDEX = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
            case 2 % silica index (from refractiveindex.info)
                load('Malitson.mat')
                INDEX = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
                
            case 3 % silicon (from refractiveindex.info)
                % without absorption
                load('Vuye-20C.mat')
                INDEX = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
            case 4 % titania meassured in the Kley group (07.2016)
                load('n-k_TiO2_Goerke.mat')
                RE_IND = pchip(IN_FILE(:,1)/1000, IN_FILE(:,2),LAMBDA);
                IM_IND = pchip(IN_FILE(:,1)/1000, IN_FILE(:,3),LAMBDA);
                
                INDEX  = RE_IND + 1i * IM_IND;
                
            case 5 % T2O5 in the range of 0.35 to 1.8 \mum (from refractiveindex.info)
                load('Gao_Ta2O5_n.mat')
                RE_IND = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
                clear IN_FILE
                load('Gao_Ta2O5_k.mat')
                IM_IND = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
                INDEX  = RE_IND + 1i * IM_IND;
                
            case 6 % amorphous T2O5 in the range of 0.5 to 1000 \mum (from refractiveindex.info)
                load('Bright-amorphous_Ta2O5_n.mat')
                RE_IND = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
                clear IN_FILE
                load('Bright-amorphous_Ta2O5_k.mat')
                IM_IND = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
                INDEX  = RE_IND + 1i * IM_IND;
                
            case 7 % Si3N4 from 0.310 to 5.504 \mum (from refractiveindex.info)
                load('Luke_SI3N4.mat')
                INDEX = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
            case 8 % amorphous MgF2 measured 2012 by ??? (V700)
                load('MGFV700')
                INDEX = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);
                
            case 9 % futurex (by Dennis, measured for ~200nm)
                load('Futurrex_IC1-200_nk.mat')
                
                RE_IND = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);               
                IM_IND = pchip(IN_FILE(:,1), IN_FILE(:,3),LAMBDA);
                
                INDEX  = RE_IND + 1i * IM_IND;
            case 10 % futurex (by FUTUREX)
                load('n-k_Futurrex.mat')
                
                RE_IND = pchip(IN_FILE(:,1), IN_FILE(:,2),LAMBDA);               
                IM_IND = pchip(IN_FILE(:,1), IN_FILE(:,3),LAMBDA);
                
                INDEX  = RE_IND + 1i * IM_IND;
                
            case 11 % Ti from Johnson, Christy Data
                DATA_IN = load('Ti_Johnson.mat');
                
                RE_IND = pchip(DATA_IN.wl, DATA_IN.n,LAMBDA);
                IM_IND = pchip(DATA_IN.wl, DATA_IN.k,LAMBDA);
                
                INDEX = RE_IND + 1i * IM_IND;
                
            case 12 % Ti from Rakic et. al. (1998) Drude-Lorentz-model
                msg = ['No data available for material ',...
                    num2str(MATERIAL)];
                error(msg);
                
            case 13 % Ti from Rakic et. al. (1998) Brendel-Bormann-model
                msg = ['No data available for material ',...
                    num2str(MATERIAL)];
                error(msg);
                
        end
        
end