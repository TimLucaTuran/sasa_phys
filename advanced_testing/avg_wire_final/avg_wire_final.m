function avg_wire_final(PARGRID_processID, PARGRID_totalNumberOfProcesses, PARGRID_storageID)

% preferences: add/change parameters for the simulation
% job disembly: create grid of parameters VARIED in the loop (spectrum)
%               --> !after that: read them out of PARGRID_CurrentLoopVec!

PATH_1 = '/lustre/vol/sperrhake/programming/calc/FMM_source';
PATH_2 = '/lustre/vol/sperrhake/programming/calc/draw_files';
PATH_3 = '/lustre/vol/sperrhake/programming/calc/materials';
addpath(PATH_1)
addpath(PATH_2)
addpath(PATH_3)

myMaxNumCompThreads(2);

% ========================== preferences ================================ %
% storage path
pathF = ['/lustre/vol/sperrhake/programming/data/',mfilename];

STORAGE_FileName = [pathF,'/', mfilename,'_Daten.mat'];

% IMPORTANT: Which geometry should be drawn
Avg_Wire = load('avg_wire_shape_4_R_smallest.mat');

% fourier orders
Order_ = 22;

% output diffraction orders
% (for only the incomming 0th order mode)
orderSMAT_bool = false;

% angle of incidence
A_theta = 0;

% parameters that can be varied in the loop
% periode in nm
Period_ = 300;

% Groove height in nm
Height_ = 45; % metal layer thickness in nm

% rotation angle of unit cell in degree
Ang_ = 0;

% wavelength in \mu
Lambda_ = linspace(0.47,1.2,512);

% save parameters in cell array for use after the simulation
PARAM_CELL = {Order_, orderSMAT_bool, A_theta, Period_, ...
   Height_, Ang_, Lambda_};

% ============================ job disembly ============================= %
V1 = Lambda_;

VAR_CELL = {V1};

PARGRID_loopID = 0;

[PARGRID_loopID, PARGRID_DataStorageFileName, PARGRID_lastFileName, ...
	PARGRID_CurrentLoopVec, PARGRID_OriginalLoopVec] = ...
    PARGRID_initializeCurrentStep( PARGRID_loopID, PARGRID_processID, ...
    PARGRID_totalNumberOfProcesses, VAR_CELL, ...
    STORAGE_FileName);

VAR_CELL   = PARGRID_CurrentLoopVec;
Lambda_    = VAR_CELL{1};

% ========================== computation core =========================== %

for loop = 1 : length(VAR_CELL{1})
    tic 
    
    % extract parameters for loop (these which are/could be varied)
    lambda_   = Lambda_(loop);

    % set angel of incidence
    A = [0 A_theta 0];
    
    % embedding material
    % Futurrex
    n_embed = real( refractive_ind(lambda_ * 1000,10,2) );
    
    % refractive indices of substrate and cladding (embedding media) 
    N=[n_embed, n_embed];
    
    % unit cell dimension x/y in {\mu}m
    D0=[Period_ Period_]/1000;
    
    % height of the layer(s)

    H =[0 Height_]/1000;
    H0 = cumsum(H);
    
    % orders in x and y direction
    OX = Order_;
    OY = Order_;
    
    % draw unit cells
    
    % gold permittivity from measured data by Kay
    % Kay
    eps_Au = Aufitek_Meep(lambda_);
    n_Au   = sqrt(eps_Au);

    % switch thresholding type
    Layer_mask = Avg_Wire.Layer_thresh_spec;

    STRUC{1} = (n_Au - n_embed) * Layer_mask + n_embed;
        
    % solve the base wave ...
    [fetat,fetar,fTxmn,fTymn,fTzmn,fRxmn,fRymn,fRzmn,eigVecs_E,eigVecs_H,...
        eigVals, coeff_AB,R,T,SMatrixTotal] = ...
        fmm3d_ml_Generalized_Rauspurzel(STRUC, H0, D0, A, N, lambda_,...
        OX, OY, 1, false);
    
    mlkm = 2*OX+1;
    nlkm = 2*OY+1;
    
    % save all outgoing higher order modes
    if orderSMAT_bool == true
    % take only the incomming 0th order mode
        sz=(2*OX+1)*(2*OY+1);
               
        % pick a polarization and a s-matrix block      
        %    xx xy yx yy
        pol_m=[0 1 0 1];
        pol_n=[0 0 1 1];
        %      Tf Rb Rf Tb
        block_m=[0 1 0 1];
        block_n=[0 0 1 1];
        
        % get linear indices and combine with block and pol position 
        % (here could be a problem if OX and OY differ)
        % take only the incomming 0th order mode
        
        row_in=sub2ind([mlkm,nlkm],OX+1,OY+1);
        orderSMat = zeros(2*OX+1,2*OY+1,4,4);
        for loop2 = 1:4
            col_in=(1:sz)+pol_m(loop2)*sz;
            for loop3 = 1:4                
            	SCOL = col_in+block_m(loop3)*2*sz;
            	SROW = row_in+pol_n(loop2)*sz+block_n(loop3)*2*sz;
                orderSMat(:,:,loop2,loop3)=reshape(SMatrixTotal(SCOL,SROW),2*OX+1,2*OY+1);
            end
        end
    else
        orderSMat = [];
    end
        
        % get fundamental mode
        
        first_ind = 0*mlkm*nlkm+fix(mlkm*nlkm/2)+1;
        last_ind  = 3*mlkm*nlkm+fix(mlkm*nlkm/2)+1;
        
        m_ind_vec = first_ind:(mlkm*nlkm):last_ind;       
        n_ind_vec = m_ind_vec;
        
        SMat0 = SMatrixTotal(m_ind_vec,n_ind_vec);
    
    % storage variables of data
    elapsedtime = toc;
    S_MATRIX=SMat0;
    ORDER_S_MAT=orderSMat;
    FETAT=squeeze(fetat(OX+1,OY+1));
    FETAR=squeeze(fetar(OX+1,OY+1));
    FTXMN=squeeze(fTxmn(OX+1,OY+1));
    FTYMN=squeeze(fTymn(OX+1,OY+1));
    FTZMN=squeeze(fTzmn(OX+1,OY+1));
    FRXMN=squeeze(fRxmn(OX+1,OY+1));
    FRYMN=squeeze(fRymn(OX+1,OY+1));
    FRZMN=squeeze(fRzmn(OX+1,OY+1));
    
    
    % decide what to save as mat-file
    if orderSMAT_bool == true
        storageData = {FETAT, FETAR, FTXMN, FTYMN, FTZMN, FRXMN, FRYMN, FRZMN, S_MATRIX, ...
        	 ORDER_S_MAT,elapsedtime,PARAM_CELL};
        storageDataNames = {'FETAT', 'FETAR', 'FTXMN', 'FTYMN', 'FTZMN', 'FRXMN', 'FRYMN', ...
        	 'FRZMN', 'S_MATRIX','ORDER_S_MAT','elapsedtime','PARAM_CELL'};
    else
        storageData = {FETAT, FETAR, FTXMN, FTYMN, FTZMN, FRXMN, FRYMN,...
        	 FRZMN, S_MATRIX,elapsedtime,PARAM_CELL};
        storageDataNames = {'FETAT', 'FETAR', 'FTXMN', 'FTYMN', 'FTZMN', 'FRXMN',...
        	 'FRYMN', 'FRZMN', 'S_MATRIX','elapsedtime','PARAM_CELL'};
    end    
    
    % save data
    PARGRID_lastFileName = PARGRID_writeCurrentResults(storageData, storageDataNames, ...
        PARGRID_loopID, PARGRID_storageID, ...
        PARGRID_CurrentLoopVec{1}, loop, ...
        PARGRID_DataStorageFileName, PARGRID_lastFileName, PARGRID_OriginalLoopVec);
end
