function [etat,etar,Txmn,Tymn,Tzmn,Rxmn,Rymn,Rzmn,eigVecs_E,eigVecs_H,eigVals,coeff_AB,R,T,SMatrixTotal,alphamn,betamn] = ...
    fmm3d_ml_Generalized_Rauspurzel(...
    STRUC,zj,d,a,n,lambda,M,N, boundarySolvingMode, calcFieldsInsideSlabB, ... % 1-10: normal input parameters
    repeatSystem, STRUCInFront, zjInFront, STRUCBehind, zjBehind, ... % 11-15:  input for "repeated system"
    eigenModeUsage, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ...
    ReflectedMode_E, ReflectedMode_H, Incident_Vec, ... % 16-23: input for generalized input and output media
    eigVecs_E, eigVecs_H, eigVals, ... % 24-26: input for known eigensolutions of the layer system
    useTemplatesInLayer) 

% function [etat,etar,Txmn,Tymn,Tzmn,Rxmn,Rymn,Rzmn,eigVecs_E,eigVecs_H,eigVals,coeff_AB,R,T] = ...
%     fmm3d_ml_Generalized(...
%     STRUC,zj,d,a,n,lambda,M,N, boundarySolvingMode, ... % 1-9: normal input parameters
%     repeatSystem, STRUCInFront, zjInFront, STRUCBehind, zjBehind, ... % 10-14:  input for "repeated system"
%     eigenModeUsage, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ...
%     ReflectedMode_E, ReflectedMode_H, ... % 15-21: input for generalized input and output media
%     eigVecs_E, eigVecs_H, eigVals) % 22-24: input for known eigensolutions of the layer system
%
% Solves 3D-diffraction problem for multilayered structure
% using Fourier modal method with Li's formulation.
%
% INPUTS:
%
% STRUC	...
%   Cell array of refractive index distribution in different layers
% zj ...
%   Vector of depth transitions points
% d	...
%   Vector of x- and y-periods
% a	...
%   Vector of angles (phi,theta & ksi) in radians
% n	...
%   Vector of incidence and exit side refractive indices
% lambda ...
%   Wavelength
% M,N ...
%   Calculate with -M:M orders in x-direction
%   Calculate with -N:N orders in y-direction
%
% boundarySolvingMode:
%   ==0: use the original (by Kettunen) implemented S-matrix-like boundary
%        value problem approach.
%   ==1: use the classical S-Matrix approach.
%   ==2: use the classical T-Matrix approach.
%
% repeatSystem ...
%   integer value >= 1 which gives the number of repetitions the structure
%   STRUC should be repeated.
%
% STRUCInFront, zjInFront ...
%   similar to <STRUC> and <zj> (see above). The structure defined here is
%   attached in front of the "normal" structure defined by <STRUC> and
%   <zj>. This allows to define a repeating "normal" structure which has a
%   leading part defined by <STRUCInFront> and <zjInFront>. Equal
%   considerations are made for a trailing structure (see below). If you
%   don't need this option simply deliver an empty structure
%   STRUCInFront={} and zjInFront=0.
% 
% STRUCBehind, zjBehind ...
%   similar to <STRUC> and <zj> (see above). The structure defined here is
%   attached behind the "normal" structure defined by <STRUC> and
%   <zj>. This allows to define a repeating "normal" structure which has a
%   trailing part defined by <STRUCBehind> and <zjBehind>. If you
%   don't need this option simply deliver an empty structure
%   STRUCBehind={} and zjBehind=0.
%
% eigenModeUsage:
%   ==0: apply ordinary homogeneous media on the input and output side of
%        the structure. Therefore all data, like the refractive indices N
%        and field polarization vectors A are used.
%   ==1: apply the ordinary homogeneous medium only on incident side.
%        The eigenmodes of the transmitted side are defined by the
%        matrices: TransmittedMode_E, TransmittedMode_H.       
%   ==2: apply the ordinary homogeneous medium only on output side.
%        The eigenmodes of the incident side are defined by the
%        matrices: ReflectedMode_E, ReflectedMode_H and IncidentMode_E,
%        IncidentMode_H.
%   ==3: both adjacent media are "user-defined". 
%        The eigenmodes are defined by the matrices: ReflectedMode_E, 
%        ReflectedMode_H and IncidentMode_E, IncidentMode_H and
%        TransmittedMode_E, TransmittedMode_H.
% 
% IncidentMode_E, IncidentMode_H ...
%   these matrices contain the ux, uy and hx, hy components of the incident
%   field expressed in terms of the eigenmodes in the incident medium.
%
% ReflectedMode_E, ReflectedMode_H ...
%   these matrices contain eigenmodes E and H in the incident medium.
%
% TransmittedMode_E, TransmittedMode_H ...
%   these matrices contain eigenmodes E and H in the output medium.
%
% eigVecs_E, eigVecs_H, eigVals ...
%   If the same physical system was run already one time before (for
%   example only for another polarization state) then the same eigenvectors
%   and eigenvalues per layer do apply also for the new calculation run. It
%   would be very time consuming to re-calculate all eigenvalues and
%   -vectors. So by delivering these structures as input parameters they
%   are used as the eigen-solutions of the system. For normal calculation
%   simply omit these parameters.
%
%
% OUTPUTS:
%
% eta_t	Efficiency of transmitted orders
% eta_r	Efficiency of reflected orders
%
% Rxmn, Rymn, Rzmn: (Amplitude) Reflection coefficients
% Txmn, Tymn, Tzmn: (Amplitude) Transmission coefficients
%
% eigVecs_E, eigVecs_H: cell array of 2D-matrices which store the
% eigenvectors (columns correspond to different eigenvalues) of (Ex,Ey,Ez)
% or (Hx, Hy, Hz) in Fourier representation for every laye r.
%
% eigVals: cell array of vectors which store the corresponding eigenvalues
% <gamma> (normal k-vector comp) for every layer.
%
% coeff_AB: cell array of vectors which store the mode weights (for the
% meaning of A or B see the originla Turunen-paper) for every layer. For
% every layer a vector is stored which contains the values for A in the
% upper part and B in the lower part - [A,B].
%
% R,T: in principal these vectors contain only R = [Rxmn, Rymn] and T =
% [Txmn, Tymn]. They are only important in some functional mode of
% eigenModeUsage > 0, where the input and/or output media are not
% homogeneous, so that the z-components (Rzmn, Tzmn) cannot simply be
% calculated out of the divergence equation divE = 0.
%
%
% Original version : V.Kettunen 1997
% Current revision : Th.Paul 2007
%
% Use of this program must be acknowledged


% Check version, must be at least Matlab 5.0
apu=version;
if apu(1)<5,error('Your Matlab version is too old!');end

global useNumericalFFTAlgorithm;
useNumericalFFTAlgorithm = true;

% =========================================================================
if (nargin < 11)
    repeatSystem = 1;
    STRUCInFront = {};
    STRUCBehind = {};
    zjInFront = 0;
    zjBehind = 0;
end

if (nargin < 16)
    eigenModeUsage = 0;
    IncidentMode_E = [];
    IncidentMode_H = [];
    ReflectedMode_E = [];
    ReflectedMode_H = [];
    TransmittedMode_E = [];
    TransmittedMode_H = [];
    Incident_Vec = {};
end

if (nargin < 24)
    useTemplates = false;
else    
    % renorm the H-field components, so now they are the "right" valued one.
    mu0 = 4*pi * 1E-07;
    eps0 = 8.85418782E-12;
    for j = 1 : length(eigVecs_H)
        eigVecs_H{j} = eigVecs_H{j} * sqrt(mu0 / eps0);
    end    
    useTemplates = true;
end

if (nargin < 27)
    useTemplatesInLayer = cell(length(STRUC)+length(STRUCBehind)+length(STRUCInFront) , 1);
    for j = 1 : length(useTemplatesInLayer)
        useTemplatesInLayer{j} = false;
    end
else
    if (length(useTemplatesInLayer) ~= length(eigVecs_H))
        error('<useTemplatesInLayer> has to match the size of <eigVecs_H>, <eigVecs_E> and <eigVals>!');
    end
end

% =========================================================================

repeatSystem = round(repeatSystem);
if (repeatSystem < 1)
    error('repeatSystem has to be an positive integer value!');
end
if (((repeatSystem ~= 1) || (~isempty(STRUCInFront)) || (~isempty(STRUCBehind))) && (boundarySolvingMode == 0))
    error('repeatSystem unequal to one, leading or trailing structures are not allowed in boundarySolvingMode = 0!');
end

if ((eigenModeUsage > 0) && (boundarySolvingMode == 0))
    error('eigenModeUsage > 0 is not allowed in boundarySolvingMode = 0!');
end

% Unpack input vectors
dx=d(1);dy=d(2);

if (iscell(a))
    if (length(a) ~= 3)
        error(['input vector <a> must be a cell array of three elements. the first is <[phi, theta]> ', ...
               'the second is a numerical vector of polarization angles <ksi(n)> and the last one is a scalar ', ...
               'identifiying the <inputFieldAmplitude>']);
    end
    
    phi=a{1}(1); 
    theta=a{1}(2); 
    ksi=num2cell(a{2}(:)); 
    iFieldAmplitude=a{3};
else
    phi=a(1);theta=a(2);ksi={a(3)};
    if (length(a) > 3); iFieldAmplitude=a(4); else iFieldAmplitude=1; end;
end

n1=n(1);n3=n(2);
if (((sign(n1) < 0) || (sign(n3) < 0)) && (boundarySolvingMode == 0))
    error('Negative index materials are not allowed in boundarySolvingMode = 0!');
end

% Number of layers
layers=length(STRUC);
if (layers+1)~=length(zj),error('Mismatch with inputs: STRUC & h');end
J=layers;

JInFront = length(STRUCInFront); 
if ((JInFront+1) ~= length(zjInFront))
    error('Mismatch with inputs: STRUCinfront & h');
end
JBehind = length(STRUCBehind); 
if ((JBehind+1) ~= length(zjBehind))
    error('Mismatch with inputs: STRUCBehind & h');
end

% Determine number of orders
mlkm=2*M+1;
nlkm=2*N+1;

% Reshape to get correct orientation
m=reshape(-M:M,mlkm,1);
n=reshape(-N:N,1,nlkm);

ux = cell(size(ksi));
uy = cell(size(ksi));
uz = cell(size(ksi));

for j = 1 : length(ksi)
    % This defines incidence field
    ux{j}=cos(ksi{j})*cos(theta)*cos(phi)-sin(ksi{j})*sin(phi);
    uy{j}=cos(ksi{j})*cos(theta)*sin(phi)+sin(ksi{j})*cos(phi);
    uz{j}=-cos(ksi{j})*sin(theta);

    % weight the incident field to the given amplitude ...
    ux{j} = ux{j}*iFieldAmplitude;
    uy{j} = uy{j}*iFieldAmplitude;
    uz{j} = uz{j}*iFieldAmplitude;
end

% Some standard parameters
% if somebody wants to calc with complex kx0(alpha0) and/or ky0(beta0)
% parameters one must change the extraction of the real part of alpha0 and
% beta0. At this point it is only done to avoid numerical problems (for
% example having some imaginary part of 1E-15 which leads to problems
% calculating the complex square root hereafter).
k=2*pi/lambda;
alpha0 = real(n1*k*sin(theta)*cos(phi)); %% ATTENTION: here the real part of kx0(alpha0) is taken
beta0 = real(n1*k*sin(theta)*sin(phi)); %% ATTENTION: here the real part of ky0(beta0) is taken
alpham=alpha0+2*pi*m/dx;
betan=beta0+2*pi*n/dy;
alphamn=alpham*ones(size(betan));
betamn=ones(size(alpham))*betan;

rmn = sign(n1) * sqrt((n1*k)^2-(alphamn).^2-(betamn).^2);
tmn = sign(n3) * sqrt((n3*k)^2-(alphamn).^2-(betamn).^2);
rmn(imag(rmn) ~= 0) = sign(imag(rmn(imag(rmn) ~= 0))) .* rmn(imag(rmn) ~= 0);
tmn(imag(tmn) ~= 0) = sign(imag(tmn(imag(tmn) ~= 0))) .* tmn(imag(tmn) ~= 0);

alphamn=alphamn(:);
betamn=betamn(:);
zindmn=find( (alphamn==alpha0) & (betamn==beta0));

% initialize fields depending on the current modus.
if (useTemplates)
    Ez = cell(length(eigVecs_E), 1);
    Hz = cell(length(eigVecs_E), 1);
    for j = 1 : length(eigVecs_E)
        if (useTemplatesInLayer{j})
           Ez{j} = eigVecs_E{j}(2*mlkm*nlkm+1 : 3*mlkm*nlkm, :);
           eigVecs_E{j} = eigVecs_E{j}(1:2*mlkm*nlkm, :);

           Hz{j} = eigVecs_H{j}(2*mlkm*nlkm+1 : 3*mlkm*nlkm, :);
           eigVecs_H{j} = eigVecs_H{j}(1:2*mlkm*nlkm, :);    
        end
    end
else
    eigVecs_E = cell(JInFront + J + JBehind,1);
    eigVecs_H = cell(JInFront + J + JBehind,1);
    eigVals   = cell(JInFront + J + JBehind,1);
    Ez = cell(JInFront + J + JBehind,1);
    Hz = cell(JInFront + J + JBehind,1);        
end

% hx = beta0*uz - rmn(zindmn)*uy;
% hy = rmn(zindmn)*ux - alpha0*uz;
% hz = alpha0*uy - beta0*ux;

clear useTemplates;
% --------------------------------------------------------------------
%  INITIALIZE THE ADJACENT LAYER MODE MATRICES ...
% --------------------------------------------------------------------
[Incident_Vec, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, TransmittedMode_E, TransmittedMode_H] = ...
    initAdjacentMediaMatrices(eigenModeUsage, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
    TransmittedMode_E, TransmittedMode_H, Incident_Vec, mlkm, nlkm, zindmn, ux, uy, uz, rmn, tmn, alpha0, beta0, alphamn, betamn, n1, n3);

% --------------------------------------------------------------------
%  CALCULATE THE SLAB MODES AND SOLVE THE BOUNDARY VALUE PROBLEM.
% --------------------------------------------------------------------
I=eye(2*mlkm*nlkm);
switch boundarySolvingMode
    case 0
        [R,T, eigVecs_E, eigVecs_H, eigVals, Ez, Hz] = useOriginalKettunenBoundaryProblemFormulation( ...
            mlkm, nlkm, betamn, alphamn, rmn, tmn, alpha0, beta0, zindmn, M,N, STRUC, k, I, ux, uy, uz, J, zj, ...
            eigVecs_E, eigVecs_H, Ez, Hz, eigVals, useTemplatesInLayer);
    case 1
        [R,T, eigVecs_E, eigVecs_H, eigVals, Ez, Hz, SMatrixTotal] = useSMatrixBoundaryProblemFormulation( ...
            Incident_Vec, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
            mlkm, nlkm, betamn, alphamn, M,N, STRUC, k, I, J, zj, ...
            repeatSystem, STRUCInFront, zjInFront, JInFront, STRUCBehind, zjBehind, JBehind, ...
            eigVecs_E, eigVecs_H, Ez, Hz, eigVals, useTemplatesInLayer);
    case 2
        [R,T, eigVecs_E, eigVecs_H, eigVals, Ez, Hz] = useTMatrixBoundaryProblemFormulation( ...
            Incident_Vec, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
            mlkm, nlkm, betamn, alphamn, M,N, STRUC, k, I, J, zj, ...
            repeatSystem, STRUCInFront, zjInFront, JInFront, STRUCBehind, zjBehind, JBehind, ...
            eigVecs_E, eigVecs_H, Ez, Hz, eigVals, useTemplatesInLayer);
end

Rxmn = cell(size(ksi));
Rymn = cell(size(ksi));
Rzmn = cell(size(ksi));
Txmn = cell(size(ksi));
Tymn = cell(size(ksi));
Tzmn = cell(size(ksi));
etar = cell(size(ksi));
etat = cell(size(ksi));

alphamn=reshape(alphamn,[mlkm,nlkm]);
betamn=reshape(betamn,[mlkm,nlkm]);

for c = 1 : length(ksi)
    if ((eigenModeUsage == 0) || (eigenModeUsage == 1))
        Rxmn{c}=R{c}(1:mlkm*nlkm);
        Rymn{c}=R{c}(mlkm*nlkm+1:2*mlkm*nlkm);
        Rxmn{c}=reshape(Rxmn{c},[mlkm,nlkm]);
        Rymn{c}=reshape(Rymn{c},[mlkm,nlkm]);
        Rzmn{c}=(alphamn.*Rxmn{c}+betamn.*Rymn{c})./rmn;

        % Calculate efficiencies
        etar{c}=real(rmn/rmn(zindmn)).*(abs(Rxmn{c}).^2+abs(Rymn{c}).^2+abs(Rzmn{c}).^2);
    else
    %     disp(['eigenModeUsage == ', num2str(eigenModeUsage),':']);
    %     disp('  Diffraction Efficencies on the input side are not calculated!');

        Rxmn{c} = []; Rymn{c} = []; Rzmn{c} = []; etar{c} = [];
    end

    if ((eigenModeUsage == 0) || (eigenModeUsage == 2))
        Txmn{c}=T{c}(1:mlkm*nlkm);
        Tymn{c}=T{c}(mlkm*nlkm+1:2*mlkm*nlkm);
        Txmn{c}=reshape(Txmn{c},[mlkm,nlkm]);
        Tymn{c}=reshape(Tymn{c},[mlkm,nlkm]);
        Tzmn{c}=-(alphamn.*Txmn{c}+betamn.*Tymn{c})./tmn;

        % Calculate efficiencies
        etat{c}=real(tmn/rmn(zindmn)).*(abs(Txmn{c}).^2+abs(Tymn{c}).^2+abs(Tzmn{c}).^2);
    else
    %     disp(['eigenModeUsage == ', num2str(eigenModeUsage),':']);
    %     disp('  Diffraction Efficencies on the output side are not calculated!');

        Txmn{c} = []; Tymn{c} = []; Tzmn{c} = []; etat{c} = [];    
    end
end

% --------------------------------------------------------------------
%  CALCULATE THE FIELD INSIDE THE SLABS
% --------------------------------------------------------------------
if (calcFieldsInsideSlabB)
    coeff_AB = calcFieldsInsideSlab( ...
        Incident_Vec, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
        T, R, J, JInFront, JBehind, mlkm, nlkm, eigVecs_E, eigVecs_H, eigVals, k, zj, zjInFront, zjBehind, I, repeatSystem);
else
    coeff_AB = cell(size(Incident_Vec));
end

% --------------------------------------------------------------------
%  SOME ADDITIONAL CALCULATIONS
% --------------------------------------------------------------------

% calc the z-components of the fields ...
for j = 1:(J+JInFront+JBehind)
    eigVecs_E{j} = [eigVecs_E{j}; Ez{j}];
    eigVecs_H{j} = [eigVecs_H{j}; Hz{j}];
end

% renorm the H-field components, so now they are the "right" valued one.
mu0 = 4*pi * 1E-07;
eps0 = 8.85418782E-12;
for j = 1:(J+JInFront+JBehind)
    eigVecs_H{j} = eigVecs_H{j} * sqrt(eps0 / mu0);
end

if (~iscell(a))
    coeff_AB = coeff_AB{1};
    
    Rxmn = Rxmn{1};
    Rymn = Rymn{1};
    Rzmn = Rzmn{1};
    Txmn = Txmn{1};
    Tymn = Tymn{1};
    Tzmn = Tzmn{1};
    etar = etar{1};
    etat = etat{1};
    
    R = R{1};
    T = T{1};
end

% end of main function


%% *************************************************************************
%%  useSMatrixBoundaryProblemFormulation,
%%       this method solves for the layers' eigenmodes and the associated
%%       boundary value problem according to the classical S-Matrix
%%       formulation.
%% *************************************************************************
%%
function [R,T, eigVecs_E, eigVecs_H, eigVals, Ez, Hz, SMatrixTotal] = useSMatrixBoundaryProblemFormulation( ...
    Incident_Vec, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
    mlkm, nlkm, betamn, alphamn, M,N, STRUCN, k, I, JN, zjN, ...
    repeatSystem, STRUCInFront, zjInFront, JInFront, STRUCBehind, zjBehind, JBehind, ...
    eigVecs_E, eigVecs_H, Ez, Hz, eigVals, useTemplatesInLayer)
   

    SMatrixTotal = eye(4*mlkm*nlkm);
%     USubMat = [-beta0*alpha0/rmn(zindmn)*eye(mlkm*nlkm), (-beta0*beta0/rmn(zindmn) - rmn(zindmn))*eye(mlkm*nlkm); ...
%                (alpha0*alpha0/rmn(zindmn) + rmn(zindmn))*eye(mlkm*nlkm), beta0*alpha0/rmn(zindmn)*eye(mlkm*nlkm)];
%     RSubMat = [diag(betamn.*alphamn./rmn(:)),...
%                diag(betamn.^2./rmn(:)+rmn(:));...
%                diag(-rmn(:)-alphamn.^2./rmn(:)),...
%                diag(-alphamn.*betamn./rmn(:))];
%     TSubMat = [diag(-betamn.*alphamn./tmn(:)),...
%                diag(-betamn.^2./tmn(:)-tmn(:));...
%                diag(tmn(:)+alphamn.^2./tmn(:)),...
%                diag(alphamn.*betamn./tmn(:))];
          
   
    % this loop has maximal to round trips.
    % in the first one the S-Matrix of the structure which is to be
    % repeated is calculated. In the second one, the extra structure in
    % front of that thing is calculated ...
    for rep = 1:3
        
        % init values for the current trip ...
        SMatrix = eye(4*mlkm*nlkm);
        switch (rep)
            case 1
                J = JInFront;
                STRUC = STRUCInFront;
                zj = zjInFront;
                jShift = 0;            
            case 2
                J = JN;
                STRUC = STRUCN;
                zj = zjN;
                jShift = JInFront;
            case 3
                J = JBehind;
                STRUC = STRUCBehind;
                zj = zjBehind;
                jShift = JInFront + JN;
            otherwise
                error('Unknown error!');
        end
                
        
    %     -----------------------------------------------------------------
    %       old version where no repeating was concidered!
    %       tInt = [I, I; USubMat, RSubMat];
    %     -----------------------------------------------------------------

        % init the incident side "tMatrix" by a identy matrix.
        % the rest will be considered afterwards, because before adding the
        % boundary fields the repeating of the Layersystem has to be performed!
        tInt = [I, zeros(size(I)); zeros(size(I)), I];

        for j = 1:J
            
           idNr = checkForIdenticalLayer(STRUC, j, true, useTemplatesInLayer);
           if (idNr == -1)
               % Calculate needed epsilon matrices for j'th layer
               [E1,E2,E3]=internal_fcg2(STRUC{j},M,N);  

               % Construct the matrix in eigenvalue problem    
               bF = [diag(alphamn)*E1*diag(betamn),...
                     k^2*eye(mlkm*nlkm)-diag(alphamn)*E1*diag(alphamn);...
                     -k^2*eye(mlkm*nlkm)+diag(betamn)*E1*diag(betamn),...
                     -diag(betamn)*E1*diag(alphamn)];
               bG = [-diag(alphamn.*(betamn)),...
                     -k^2*E3+diag(alphamn.^2);...
                     k^2*E2-diag(betamn.^2),...
                     diag(alphamn.*(betamn))];
               WL= bF*bG;

               % Solve eigenvalue problem
               [WL,gamma]=eig(WL);

               % Scaling & sign convention for gamma
               gamma=gamma/k^2;
               gamma=sqrt(gamma);
               gamma=gamma.*sign(real(gamma)+imag(gamma));
               
               gamma=diag(gamma);
               [WL, gamma] = sortEigenStatesUpDown(WL, gamma);

               %save eigen WL gamma
               eigVecs_E{j+jShift} = WL;
               eigVecs_H{j+jShift} = bG * WL ./ repmat(k*gamma.', [size(WL,1), 1]);
               
               eigVals{j+jShift} = gamma;

               % calc z-comp. Here by the prefactor it is ensured, that Hz is always
               % in the final version. So the following renoprmalization has to be
               % performed only for Hx, Hy.
               Ez{j+jShift} = -1 / k * E1 * (diag(alphamn(:)) * eigVecs_H{j+jShift}(mlkm*nlkm+1 : 2*mlkm*nlkm, :) - diag(betamn(:)) * eigVecs_H{j+jShift}(1 : mlkm*nlkm, :));
               Hz{j+jShift} = +1 / k * (diag(alphamn(:)) * eigVecs_E{j+jShift}(mlkm*nlkm+1 : 2*mlkm*nlkm, :) - diag(betamn(:)) * eigVecs_E{j+jShift}(1 : mlkm*nlkm, :));       
           else
               eigVecs_E{j+jShift} = eigVecs_E{idNr+jShift};
               eigVecs_H{j+jShift} = eigVecs_H{idNr+jShift};
               eigVals{j+jShift} = eigVals{idNr+jShift};
               Ez{j+jShift} = Ez{idNr+jShift};
               Hz{j+jShift} = Hz{idNr+jShift};
               
               WL = eigVecs_E{j+jShift};
               gamma = eigVals{j+jShift};
           end
           
           tInt = [WL, WL; k*eigVecs_H{j+jShift}, -k*eigVecs_H{j+jShift}] \ tInt;
           XL = diag(exp(i*gamma*(zj(j+1)-zj(j))));       

           % built up the complete layer-s-matrix ...
           sLayer = [XL, zeros(size(I)); zeros(size(I)), I] * convertInterfaceTinS(tInt) * [I, zeros(size(I)); zeros(size(I)), XL];

           % append sLayer to the overall S-Matrix ...
           SMatrix = addSmatrix(SMatrix, sLayer, I);

           tInt = [WL, WL; k*eigVecs_H{j+jShift}, -k*eigVecs_H{j+jShift}];
        end

    %     -----------------------------------------------------------------
    %       old version where no repeating was concidered!
    %       tInt = [I, I; TSubMat, I] \ tInt;
    %       SMatrix = addSmatrix(SMatrix, convertInterfaceTinS(tInt), I);
    %     -----------------------------------------------------------------

        % add only the eigenvalue matrix tInt of the last layer to the current
        % S-matrix (so we end up at the right side of the last interface).
        % because we started the matrix also on the left side of the left
        % interface we have now the correct matrix which describes exactly one
        % layersystem or "unit cell"...
        SMatrix = addSmatrix(SMatrix, convertInterfaceTinS(tInt), I);
        
        switch (rep)
            case 1
                % nothing to do            
            case 2
                % repeat the unit cell to the proper size.
                SMatrix = repeatSMatrix(SMatrix, repeatSystem, I);
            case 3
                % nothing to do
            otherwise
                error('Unknown error!');
        end
        
        % append results to the final S-Matrix.
        SMatrixTotal = addSmatrix(SMatrixTotal, SMatrix, I);
    end
    
    clear SMatrix;

    
%     R = cell(size(Incident_Vec)); 
%     T = cell(size(Incident_Vec)); 
%     for c = 1 : length(Incident_Vec)
%         [cT, cR] = incorporateBoundaryMedia(Incident_Vec{c}, zeros(size(TransmittedMode_E,2),1), SMatrixTotal, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, TransmittedMode_E, TransmittedMode_H, I);
%         T{c} = cT;
%         R{c} = cR;
%     end        

    % now append also the "effect" of the transmitted side to the S-matrix.
    tInt = [TransmittedMode_E, TransmittedMode_E; TransmittedMode_H, -TransmittedMode_H] \ eye(2*size(I));
    [sLayer, midPointSLayer] = convertInterfaceTinS(tInt, size(tInt)/2);
    [SMatrixTotal, midPointS] = addSmatrix(SMatrixTotal, sLayer, I, size(SMatrixTotal)/2, midPointSLayer);

    % and finally the incident side has to be considered and has to
    % appended (but in front!) to the current S-matrix.
    tInt = [IncidentMode_E, ReflectedMode_E; IncidentMode_H, ReflectedMode_H];
    [sLayer, midPointSLayer] = convertInterfaceTinS(tInt, size(tInt)/2);
    SMatrixTotal = addSmatrix(sLayer, SMatrixTotal, I, midPointSLayer, midPointS);
   
    R = cell(size(Incident_Vec)); 
    T = cell(size(Incident_Vec)); 
    for c = 1 : length(Incident_Vec)
        B=[Incident_Vec{c}; zeros(2*mlkm*nlkm,1)];

        RTI = SMatrixTotal * B;
        T{c}=RTI(1:2*mlkm*nlkm);
        R{c}=RTI(2*mlkm*nlkm+1 : 4*mlkm*nlkm);
    end
%endfunction



%% *************************************************************************
%%  incorporateBoundaryMedia,
%%       this method incorporates the boundary media.
%%       this function is specially designed for the case, that the 
%%       eigenvector matrices for the boundary media are no square-
%%       matrices. 
%%       This function causes problems in the normal case, where we deal
%%       with quadaratic matrices and ordinary homog. input and output
%%       media. Problems in this connection means: that some of the block
%%       matrices which are incorporated in the calculation have no full
%%       rank. Especially the blocks of the S-Matrix itself are not
%%       "zwangsl�ufig" full-rank matrices. It semms to be (from numerical
%%       experience) that only s12 and s21 are fullrank, but not s11, s22.
%%       In the usual code, where S-matrices are connected to each other
%%       this fact causes no problem (why: good question, but it must be
%%       the appearance of the identity matrix inside the code (see the
%%       associated function, which causes everything to go right!) But in
%%       this function we have problems: because s11 and s22 are used to
%%       calc the M and Q matrices (see below) which makes them to be no
%%       full-rank matrices. This for itself makes problems during the
%%       inversion which causes numerical errors in the overall
%%       calculation.
%%       FAZIT: This function is obsolete and the "old" function to
%%       concatinate s-matrices (<addSMatrix>) together with the conversion
%%       function for T in S matrices <> is adapated to work also with
%%       non-quadratic matrices. In the usual case of quadratic matrices
%%       this code is stable, of course. In the case of non-quadratic
%%       matrices it needs to be investigated further.
%% *************************************************************************
%%
function [F3, B0] = incorporateBoundaryMedia(F0, B3, SMatrix, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, TransmittedMode_E, TransmittedMode_H, I)

    T2 = [TransmittedMode_E, eye(size(TransmittedMode_E)); TransmittedMode_H, eye(size(TransmittedMode_H))];
    T1 = [IncidentMode_E, ReflectedMode_E; IncidentMode_H, ReflectedMode_H];
    
    si = size(T1)/2;    
    
    t1_11 = T1(1:si(1), 1:si(2));
    t1_12 = T1(1:si(1), si(2)+1 : 2*si(2));
    t1_21 = T1(si(1)+1 : 2*si(1), 1:si(2));
    t1_22 = T1(si(1)+1 : 2*si(1), si(2)+1 : 2*si(2));
    
    si = size(T2)/2;    
    
    t2_11 = T2(1:si(1), 1:si(2));
    t2_12 = T2(1:si(1), si(2)+1 : 2*si(2));
    t2_21 = T2(si(1)+1 : 2*si(1), 1:si(2));
    t2_22 = T2(si(1)+1 : 2*si(1), si(2)+1 : 2*si(2));
       
    si = size(SMatrix,1)/2;    
    
    s11 = SMatrix(1:si, 1:si);
    s12 = SMatrix(1:si, si+1 : 2*si);
    s21 = SMatrix(si+1 : 2*si, 1:si);
    s22 = SMatrix(si+1 : 2*si, si+1 : 2*si);
    
    clear T1 T2 SMatrix;
    clear IncidentMode_E IncidentMode_H ReflectedMode_E ReflectedMode_H TransmittedMode_E TransmittedMode_H;
    
    M = s22 * t2_21;
    N = s12 * t2_21 - t2_11;
    P = s21 * t1_12 - t1_22;
    Q = s11 * t1_12;
    
    F3Mat = (P\s22 - Q\s12) * t2_21 + Q\t2_11;
    
    F3 = (F3Mat \ ((Q\s11 - P\s21) * t1_11 + P\t1_21)) * F0 + ...
         (F3Mat \ ((Q\s12 - P\s22) * t2_22 - Q\t2_12)) * B3;
     
    B0Mat = (M\s21 - N\s11) * t1_12 - M\t1_22;
    
    B0 = (B0Mat \ ((N\s11 - M\s21) * t1_11 + M\t1_21)) * F0 + ...
         (B0Mat \ ((N\s12 - M\s22) * t2_22 - N\t2_12)) * B3;

% endfunction



%% *************************************************************************
%%  useTMatrixBoundaryProblemFormulation,
%%       this method solves for the layers' eigenmodes and the associated
%%       boundary value problem according to the classical T-Matrix
%%       formulation.
%% *************************************************************************
%%
function [R,T, eigVecs_E, eigVecs_H, eigVals, Ez, Hz] = useTMatrixBoundaryProblemFormulation( ...
    Incident_Vec, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
    mlkm, nlkm, betamn, alphamn, M,N, STRUCN, k, I, JN, zjN, ...
    repeatSystem, STRUCInFront, zjInFront, JInFront, STRUCBehind, zjBehind, JBehind, ...
    eigVecs_E, eigVecs_H, Ez, Hz, eigVals, useTemplatesInLayer)


    TMatrixTotal = eye(4*mlkm*nlkm);
%     USubMat = [-beta0*alpha0/rmn(zindmn)*eye(mlkm*nlkm), (-beta0*beta0/rmn(zindmn) - rmn(zindmn))*eye(mlkm*nlkm); ...
%                (alpha0*alpha0/rmn(zindmn) + rmn(zindmn))*eye(mlkm*nlkm), beta0*alpha0/rmn(zindmn)*eye(mlkm*nlkm)];
%     RSubMat = [diag(betamn.*alphamn./rmn(:)),...
%                diag(betamn.^2./rmn(:)+rmn(:));...
%                diag(-rmn(:)-alphamn.^2./rmn(:)),...
%                diag(-alphamn.*betamn./rmn(:))];
%     TSubMat = [diag(-betamn.*alphamn./tmn(:)),...
%                diag(-betamn.^2./tmn(:)-tmn(:));...
%                diag(tmn(:)+alphamn.^2./tmn(:)),...
%                diag(alphamn.*betamn./tmn(:))];
              
    % this loop has maximal to round trips.
    % in the first one the S-Matrix of the structure which is to be
    % repeated is calculated. In the second one, the extra structure in
    % front of that thing is calculated ...
    for rep = 1:3
        
        % init values for the current trip ...
        TMatrix = eye(4*mlkm*nlkm);
        switch (rep)
            case 1
                J = JInFront;
                STRUC = STRUCInFront;
                zj = zjInFront;
                jShift = 0;            
            case 2
                J = JN;
                STRUC = STRUCN;
                zj = zjN;
                jShift = JInFront;
            case 3
                J = JBehind;
                STRUC = STRUCBehind;
                zj = zjBehind;
                jShift = JInFront + JN;
            otherwise
                error('Unknown error!');
        end
        
        % tInt = [I, I; USubMat, RSubMat];    

        % init the incident side "tMatrix" by a identy matrix.
        % the rest will be considered afterwards, because before adding the
        % boundary fields the repeating of the Layersystem has to be performed!
        tInt = [I, zeros(size(I)); zeros(size(I)), I];

        for j = 1:J
            
           idNr = checkForIdenticalLayer(STRUC, j, true, useTemplatesInLayer);
           if (idNr == -1) 
               % Calculate needed epsilon matrices for j'th layer
               [E1,E2,E3]=internal_fcg2(STRUC{j},M,N);  

               % Construct the matrix in eigenvalue problem    
               bF = [diag(alphamn)*E1*diag(betamn),...
                     k^2*eye(mlkm*nlkm)-diag(alphamn)*E1*diag(alphamn);...
                     -k^2*eye(mlkm*nlkm)+diag(betamn)*E1*diag(betamn),...
                     -diag(betamn)*E1*diag(alphamn)];
               bG = [-diag(alphamn.*(betamn)),...
                     -k^2*E3+diag(alphamn.^2);...
                     k^2*E2-diag(betamn.^2),...
                     diag(alphamn.*(betamn))];
               WL= bF*bG;

               % Solve eigenvalue problem
               [WL,gamma]=eig(WL);

               % Scaling & sign convention for gamma
               gamma=gamma/k^2;
               gamma=sqrt(gamma);
               gamma=gamma.*sign(real(gamma)+imag(gamma));       

               %save eigen WL gamma
               eigVecs_E{j+jShift} = WL;
               eigVecs_H{j+jShift} = bG * WL ./ repmat(k*diag(gamma).', [size(WL,1), 1]);

               gamma=diag(gamma);
               eigVals{j+jShift} = gamma;

               % calc z-comp. Here by the prefactor it is ensured, that Hz is always
               % in the final version. So the following renoprmalization has to be
               % performed only for Hx, Hy.
               Ez{j+jShift} = -1 / k * E1 * (diag(alphamn(:)) * eigVecs_H{j+jShift}(mlkm*nlkm+1 : 2*mlkm*nlkm, :) - diag(betamn(:)) * eigVecs_H{j+jShift}(1 : mlkm*nlkm, :));
               Hz{j+jShift} = +1 / k * (diag(alphamn(:)) * eigVecs_E{j+jShift}(mlkm*nlkm+1 : 2*mlkm*nlkm, :) - diag(betamn(:)) * eigVecs_E{j+jShift}(1 : mlkm*nlkm, :));
           else
               eigVecs_E{j+jShift} = eigVecs_E{idNr+jShift};
               eigVecs_H{j+jShift} = eigVecs_H{idNr+jShift};
               eigVals{j+jShift} = eigVals{idNr+jShift};
               Ez{j+jShift} = Ez{idNr+jShift};
               Hz{j+jShift} = Hz{idNr+jShift};
               
               WL = eigVecs_E{j+jShift};
               gamma = eigVals{j+jShift};               
           end
           
           tInt = [WL, WL; k*eigVecs_H{j+jShift}, -k*eigVecs_H{j+jShift}] \ tInt;
           XL = diag(exp(i*gamma*(zj(j+1)-zj(j))));
           XLn = diag(exp(-i*gamma*(zj(j+1)-zj(j))));

           TMatrix = [XL, zeros(size(I)); zeros(size(I)), XLn] * tInt * TMatrix;

           tInt = [WL, WL; k*eigVecs_H{j+jShift}, -k*eigVecs_H{j+jShift}];
        end

        % add only the eigenvalue matrix tInt of the last layer to the current
        % S-matrix (so we end up at the right side of the last interface).
        % because we started the matrix also on the left side of the left
        % interface we have now the correct matrix which describes exactly one
        % layersystem or "unit cell"...    
        TMatrix = tInt * TMatrix;
        
        switch (rep)
            case 1
                % nothing to do            
            case 2
                % repeat the unit cell to the proper size.
                TMatrix = TMatrix^(repeatSystem-1);
            case 3
                % nothing to do                            
            otherwise
                error('Unknown error!');
        end        
        
        TMatrixTotal = TMatrix * TMatrixTotal;
    end
    
    clear TMatrix;
    
    % now append also the "effect" of the transmitted side to the T-matrix.
    tInt = [TransmittedMode_E, eye(size(TransmittedMode_E)); TransmittedMode_H, eye(size(TransmittedMode_H))] \ eye(2*size(I));
    TMatrixTotal = tInt * TMatrixTotal;
    
    % and finally the incident side has to be considered and has to
    % appended (but in front!) to the current T-matrix (Attention:
    % physically "in front" means matrix multiplkication from right).
    tInt = [IncidentMode_E, ReflectedMode_E; IncidentMode_H, ReflectedMode_H];
    TMatrixTotal = TMatrixTotal * tInt;    
    
    % now at the end the T-Matrix is converted into an S-matrix, so that we
    % can use the incident field as the input vector for the calculation.
    % This conversion is only an rearrangement of the underlying linear
    % equation system and does not change anything in the principles of the
    % previous T-matrix calculations and their numerical properties.
    sMat = convertInterfaceTinS(TMatrixTotal, size(TMatrixTotal)/2);
    
    R = cell(size(Incident_Vec)); 
    T = cell(size(Incident_Vec)); 
    for c = 1 : length(Incident_Vec)
        B=[Incident_Vec{c}; zeros(2*mlkm*nlkm,1)];
        RTI = sMat * B;
        T{c}=RTI(1:2*mlkm*nlkm);
        R{c}=RTI(2*mlkm*nlkm+1 : 4*mlkm*nlkm);
    end
%endfunction


%% *************************************************************************
%%  useOriginalKettunenBoundaryProblemFormulation,
%%       this method solves for the layers' eigenmodes and the associated
%%       boundary value problem according to the original formulation by
%%       Ketunnen.
%% *************************************************************************
%%
function [R,T, eigVecs_E, eigVecs_H, eigVals, Ez, Hz] = useOriginalKettunenBoundaryProblemFormulation( ...
    mlkm, nlkm, betamn, alphamn, rmn, tmn, alpha0, beta0, zindmn, M,N, STRUC, k, I, ux, uy, uz, J, zj, ...
    eigVecs_E, eigVecs_H, Ez, Hz, eigVals, useTemplatesInLayer)

    A=I;
    FG=[eye(2*mlkm*nlkm);...
          [diag(-betamn.*alphamn./tmn(:)),...
             diag(-betamn.^2./tmn(:)-tmn(:));...
             diag(tmn(:)+alphamn.^2./tmn(:)),...
             diag(alphamn.*betamn./tmn(:))]];

    % This is the main loop
    for j=J:-1:1,
        
       idNr = checkForIdenticalLayer(STRUC, j, false, useTemplatesInLayer);       
       if (idNr == -1)
           % Calculate needed epsilon matrices for j'th layer
           [E1,E2,E3]=internal_fcg2(STRUC{j},M,N);  

           % Construct the matrix in eigenvalue problem    
           bF = [diag(alphamn)*E1*diag(betamn),...
                 k^2*eye(mlkm*nlkm)-diag(alphamn)*E1*diag(alphamn);...
                 -k^2*eye(mlkm*nlkm)+diag(betamn)*E1*diag(betamn),...
                 -diag(betamn)*E1*diag(alphamn)];
           bG = [-diag(alphamn.*(betamn)),...
                 -k^2*E3+diag(alphamn.^2);...
                 k^2*E2-diag(betamn.^2),...
                 diag(alphamn.*(betamn))];
           WL= bF*bG;

           % Solve eigenvalue problem
           [WL,gamma]=eig(WL);

           % Scaling & sign convention for gamma
           gamma=gamma/k^2;
           gamma=sqrt(gamma);
           gamma=gamma.*sign(real(gamma)+imag(gamma));   

           %save eigen WL gamma
           eigVecs_E{j} = WL;
           eigVecs_H{j} = bG * WL ./ repmat(k*diag(gamma).', [size(WL,1), 1]);
           VL = k*eigVecs_H{j};

           gamma=diag(gamma);
           eigVals{j} = gamma;
           % Implement boundary conditions between layers (numerically stabile formulation)

           % calc z-comp. Here by the prefactor it is ensured, that Hz is always
           % in the final version. So the following renoprmalization has to be
           % performed only for Hx, Hy.
           Ez{j} = -1 / k * E1 * (diag(alphamn(:)) * eigVecs_H{j}(mlkm*nlkm+1 : 2*mlkm*nlkm, :) - diag(betamn(:)) * eigVecs_H{j}(1 : mlkm*nlkm, :));
           Hz{j} = +1 / k * (diag(alphamn(:)) * eigVecs_E{j}(mlkm*nlkm+1 : 2*mlkm*nlkm, :) - diag(betamn(:)) * eigVecs_E{j}(1 : mlkm*nlkm, :));
       else
               eigVecs_E{j} = eigVecs_E{idNr};
               eigVecs_H{j} = eigVecs_H{idNr};
               eigVals{j} = eigVals{idNr};
               Ez{j} = Ez{idNr};
               Hz{j} = Hz{idNr};
               
               WL = eigVecs_E{j};
               VL = k*eigVecs_H{j};
               gamma = eigVals{j};            
       end
       
       % this is the diagonal matrix: exp{i gamma h}
       XL=diag(exp(i*gamma*(zj(j+1)-zj(j))));

       % this is the product of the layer eigenmode matrix [E E; H -H] and
       % current right border field distribution. In the initial case (j=J)
       % this corresponds to the "transmitted field" ...
       albl=[WL WL;VL -VL]\FG;

       ial=inv(albl(1:2*mlkm*nlkm,1:2*mlkm*nlkm));
       bl=albl(2*mlkm*nlkm+1:4*mlkm*nlkm,1:2*mlkm*nlkm);
       FG=[WL*(I+XL*bl*ial*XL);VL*(I-XL*bl*ial*XL)];
       A=A*ial*XL;
    end

    % Implement boundary conditions for input and output layers
    R_M=[eye(2*mlkm*nlkm);[...
             diag(betamn.*alphamn./rmn(:)),...
             diag(betamn.^2./rmn(:)+rmn(:));...
             diag(-rmn(:)-alphamn.^2./rmn(:)),...
             diag(-alphamn.*betamn./rmn(:))]];
        
    R = cell(size(ux));
    T = cell(size(ux));
    for c = 1 : length(ux)
        B=zeros(4*mlkm*nlkm,1);
        B(N*mlkm+M+1)=ux{c};
        B(mlkm*nlkm+N*mlkm+M+1)=uy{c};
        B(2*mlkm*nlkm+N*mlkm+M+1)=(beta0*uz{c}-rmn(zindmn)*uy{c});
        B(3*mlkm*nlkm+N*mlkm+M+1)=(rmn(zindmn)*ux{c}-alpha0*uz{c});

        % Solve boundary problem
        RTI=[-R_M FG]\B;

        % Coefficients for reflected and transmitted orders
        R{c}=RTI(1:2*mlkm*nlkm);
        T{c}=A*RTI(2*mlkm*nlkm+1:4*mlkm*nlkm);
    end   
%endfunction


%% *************************************************************************
%%  calcFieldsInsideSlab,
%%       this method solves for the layers' eigenmodes' coefficients, which
%%       with the field distribution inside the layer structure can be
%%       calculated.
%% *************************************************************************
%%
function coeff_AB = calcFieldsInsideSlab( ...
    Incident_Vec, TransmittedMode_E, TransmittedMode_H, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
    Tl, Rl, JN, JInFront, JBehind, mlkm, nlkm, eigVecs_E, eigVecs_H, eigVals, k, zjN, zjInFront, zjBehind, I, repeatSystem)

    %%% In this section the field distribution inside the layers is
    %%% calculated. So let's try to give the idea of this calculation: In
    %%% principle it shouldn't be any problem to calculate the A and B
    %%% coefficients inside the layers, when the reflection and
    %%% transmission coefficients are known. But during "propagation"
    %%% through the slabs one needs to use the propagation kernels 
    %%% exp(-i * gamma * z) and exp(+i * gamma * z), respectively. The
    %%% first one has the problem that it grows for evanescent waves and
    %%% the second one goes to zero (--> numerical representation problem).
    %%% In principle the property of the second one (going to zero) is no
    %%% serious problem but during "propagation" by T-matrix we need to
    %%% invert this factor (so we need the first one). Therefore in
    %%% numerical representation we have the problem that zero*infinity
    %%% must be equal to one. The best way to avoid the usage of any kernel
    %%% of the type exp(-i * gamma * z) is to use a S-matrix algorithm for
    %%% this problem. So the underlying algorithm does the following: 
    %%%     1. By simple Interface-T-Matrix the coefficient vector 
    %%%        [A_1, B_1*exp(i*gamma*h)] of the first layer (coming from the
    %%%        incident side) and the vector [A_L*exp(i*gamma*h), B_L] of
    %%%        the last layer L (coming from transmitted side) can be
    %%%        calculated. So we have direct access to the parameters A_1
    %%%        (foward going field in the first layer) and B_L (backward
    %%%        going field in the layer L).
    %%%     2. If we now represent the remaining structure (layers from 2
    %%%        to L-1) by an S-matrix we can see that A1 and B_L determine
    %%%        the input vector for this S-matrix problem. More in detail
    %%%        the equations are as follows: 
    %%%        [A_L, B_1] = S * [A_1*exp(i*gamma*h), B_L*exp(i*gamma*h)].
    %%%        By this procedure we have now the complete coefficient set
    %%%        of the first (1) and the last (L) layer.
    %%%     3. Now by returning to point 1. we can go to the second and the
    %%%        second last (L-1) layer and use the whole procedure
    %%%        recursively until we meet in the middle of the layer system.

    coeff_AB = cell(size(Incident_Vec));
    for polState = 1 : length(Incident_Vec)
        
        thicknessTotal = [diff(zjInFront), diff(zjN), diff(zjBehind)];

        coeff_AB{polState} = cell(JInFront + JN*repeatSystem + JBehind,1);

        % init values for the current left/right layer.
        leftLayer = 1;
        rightLayer = JInFront + JN*repeatSystem + JBehind;

        % FRight is a vector of length 4*M*N which contains the Fourier
        % coefficients (Ex, Ey, Hx, Hy) of the field on the right interface
        % (--> init value: transmitted field).
        FRight = [TransmittedMode_E, eye(size(TransmittedMode_E)); TransmittedMode_H, eye(size(TransmittedMode_H))] * [Tl{polState}; zeros(length(Tl{polState}), 1)];

        % FLeft is a vector of length 4*M*N which contains the Fourier
        % coefficients (Ex, Ey, Hx, Hy) of the field on the left interface
        % (--> init value: incident + reflected field).    
        FLeft = [IncidentMode_E, ReflectedMode_E; IncidentMode_H, ReflectedMode_H] * [Incident_Vec{polState}; Rl{polState}];

        % loop over the structure, until the left and right layer counter met
        % in the middle of the structure.
        while (leftLayer <= rightLayer)

            % init coefficient vector for the current layers.
            coeff_AB{polState}{leftLayer} = zeros(4*mlkm*nlkm,1);
            coeff_AB{polState}{rightLayer} = zeros(4*mlkm*nlkm,1);

            %read in the stored eigenvectors ...
            [WL, VL] = getFieldAndStructureDataForSpecificLayer(rightLayer, JN, JInFront, JBehind, repeatSystem, ...
                eigVecs_E, eigVecs_H, eigVals, thicknessTotal);
            VL = k*VL;

            % calc the coefficient vector [A_ri*exp(i*gamma*h), B_ri]. B_ri
            % can be extracted directly and is saved to the coefficient
            % vector ...
            AexpB = [WL WL;VL -VL] \ FRight;
            coeff_AB{polState}{rightLayer}(2*mlkm*nlkm+1 : 4*mlkm*nlkm) = AexpB(2*mlkm*nlkm+1 : 4*mlkm*nlkm);

            [WL, VL] = getFieldAndStructureDataForSpecificLayer(leftLayer, JN, JInFront, JBehind, repeatSystem, ...
                eigVecs_E, eigVecs_H, eigVals, thicknessTotal);
            VL = k*VL;

            % calc the coefficient vector [A_le, B_le*exp(i*gamma*h)]. A_le
            % can be extracted disrectly and is saved to the coefficient
            % vector ...        
            ABexp = [WL WL;VL -VL] \ FLeft;
            coeff_AB{polState}{leftLayer}(1 : 2*mlkm*nlkm) = ABexp(1 : 2*mlkm*nlkm);

            % now we need to built up the S-Matrix for the sandwiched
            % structure. So in principle we do this layer by layer. We
            % calculate the layer S-matrix which includes the propagation from
            % the (left+1)- side to the (right-1)+ side. After this the effect
            % of the (left+0)+ side and the (right-0)- side are appended in the
            % next if-statement.

            SMat = createSMatrixBetweenGivenLayersFromStoredData(eigVecs_E, eigVecs_H, eigVals, ...
                leftLayer+1, rightLayer-1, JN, JInFront, JBehind, thicknessTotal, ...
                mlkm, nlkm, repeatSystem, I, k);

            if (leftLayer < rightLayer)

                % calculate the eigenmode matrix of the left layer and append
                % it to the S-matrix, see below ...
                [WL, VL] = getFieldAndStructureDataForSpecificLayer(leftLayer, JN, JInFront, JBehind, repeatSystem, ...
                    eigVecs_E, eigVecs_H, eigVals, thicknessTotal);
                VL = k*VL;

                tIntLeft = [WL WL;VL -VL];

                % calculate the inverse eigenmode matrix of the right layer and append
                % it to the S-matrix, see below ...            
                [WL, VL] = getFieldAndStructureDataForSpecificLayer(rightLayer, JN, JInFront, JBehind, repeatSystem, ...
                    eigVecs_E, eigVecs_H, eigVals, thicknessTotal);
                VL = k*VL;

                tIntRight = [WL WL;VL -VL] \ eye(2*size(I));

                SMat = addSmatrix(convertInterfaceTinS(tIntLeft), SMat, I);
                SMat = addSmatrix(SMat, convertInterfaceTinS(tIntRight), I);
            end

            % now the coefficients that we are looking for are extracted
            % multiplying the S-matrix on the according input vector which
            % gives the equation:
            % [A_ri, B_le] = S * [A_le*exp(i*gamma*h), B_ri*exp(i*gamma*h)]

            if (leftLayer < rightLayer)

                % this vector will contain the input values:
                % [A_le*exp(i*gamma*h), B_ri*exp(i*gamma*h)].
                IN = zeros(4*mlkm*nlkm,1);

                [WL, VL, thickness, gamma] = getFieldAndStructureDataForSpecificLayer(leftLayer, JN, JInFront, JBehind, repeatSystem, ...
                    eigVecs_E, eigVecs_H, eigVals, thicknessTotal);

                XL = diag(exp(i*gamma*thickness));
                IN(1 : 2*mlkm*nlkm) = XL * coeff_AB{polState}{leftLayer}(1 : 2*mlkm*nlkm);

                [WL, VL, thickness, gamma] = getFieldAndStructureDataForSpecificLayer(rightLayer, JN, JInFront, JBehind, repeatSystem, ...
                    eigVecs_E, eigVecs_H, eigVals, thicknessTotal);

                XL = diag(exp(i*gamma*thickness));
                IN(2*mlkm*nlkm+1 : 4*mlkm*nlkm) = XL * coeff_AB{polState}{rightLayer}(2*mlkm*nlkm+1 : 4*mlkm*nlkm);

                % solve the boundary value problem and store the received
                % coefficients ...
                ArBl = SMat * IN;
                coeff_AB{polState}{leftLayer}(2*mlkm*nlkm+1 : 4*mlkm*nlkm) = ArBl(2*mlkm*nlkm+1 : 4*mlkm*nlkm);        
                coeff_AB{polState}{rightLayer}(1 : 2*mlkm*nlkm) = ArBl(1 : 2*mlkm*nlkm);

                % now the vectors (inhomogenities) for the left and right field
                % distribution must be updated to reflect the distributions of
                % one layer onwardly.

                % update the field to the right side distribution inside the
                % (current) left layer ...
                [WL, VL] = getFieldAndStructureDataForSpecificLayer(leftLayer, JN, JInFront, JBehind, repeatSystem, ...
                    eigVecs_E, eigVecs_H, eigVals, thicknessTotal);
                VL = k*VL;
                FLeft(1 : 2*mlkm*nlkm) = IN(1 : 2*mlkm*nlkm);
                FLeft(2*mlkm*nlkm+1 : 4*mlkm*nlkm) = coeff_AB{polState}{leftLayer}(2*mlkm*nlkm+1 : 4*mlkm*nlkm);
                FLeft = [WL WL;VL -VL] * FLeft;

                % update the field to the left side distribution inside the
                % (current) right layer ...            
                [WL, VL] = getFieldAndStructureDataForSpecificLayer(rightLayer, JN, JInFront, JBehind, repeatSystem, ...
                    eigVecs_E, eigVecs_H, eigVals, thicknessTotal);
                VL = k*VL;
                FRight(2*mlkm*nlkm+1 : 4*mlkm*nlkm) = IN(2*mlkm*nlkm+1 : 4*mlkm*nlkm);
                FRight(1 : 2*mlkm*nlkm) = coeff_AB{polState}{rightLayer}(1 : 2*mlkm*nlkm);
                FRight = [WL WL;VL -VL] * FRight;
            end

            % go to the next step ...
            leftLayer = leftLayer + 1;
            rightLayer = rightLayer - 1;
        end   
    end %polState

%endfunction



%% *************************************************************************
%%  createSMatrixBetweenGivenLayersFromStoredData,
%%       this creates a S-matrix between the given layers 
%%       The data is extracted out of the data structures
%%       eigVecs_E, eigVecs_H, eigVals.
%% *************************************************************************
%%
function SMatrixTotal = createSMatrixBetweenGivenLayersFromStoredData(eigVecs_E, eigVecs_H, eigVals, ...
    leftLayer, rightLayer, JN, JInFront, JBehind, ...
    thicknessTotal, mlkm, nlkm, repeatSystem, I, k)

    SMatrixTotal = eye(4*mlkm*nlkm);
    repSystemLocal = fix(double((min(rightLayer-JInFront, repeatSystem*JN) + JInFront+1) - max(JInFront+1, leftLayer)) / JN);
    repSystemLocal = max(repSystemLocal, 1);
      
    % this loop has maximal to round trips.
    % in the first one the S-Matrix of the structure which is to be
    % repeated is calculated. In the second one, the extra structure in
    % front of that thing is calculated ...
    for rep = 1:3
        
        % init values for the current trip ...
        SMatrix = eye(4*mlkm*nlkm);
        switch (rep)
            case 1
                % structure in front ...
                lIndex = leftLayer;
                rIndex = min(JInFront, rightLayer);
            case 2                
                 lIndex = max(JInFront+1, leftLayer);
                 rIndex = min(lIndex + JN - 1, rightLayer);
            case 3
                 lIndex = max(JInFront+1, leftLayer) + repSystemLocal*JN + 0;
                 rIndex = rightLayer;
            otherwise
                error('Unknown error!');
        end
                
        % init the incident side "tMatrix" by a identy matrix.
        % the rest will be considered afterwards, because before adding the
        % boundary fields the repeating of the Layersystem has to be performed!
        tInt = [I, zeros(size(I)); zeros(size(I)), I];

        for j = lIndex : +1 : rIndex

           [eigE, eigH, thickness, eigval] = ...
               getFieldAndStructureDataForSpecificLayer(j, JN, JInFront, JBehind, ...
               repeatSystem, eigVecs_E, eigVecs_H, eigVals, thicknessTotal);
            
           tInt = [eigE, eigE; k*eigH, -k*eigH] \ tInt;
           XL = diag(exp(i*eigval*thickness));       

           % built up the complete layer-s-matrix ...
           sLayer = [XL, zeros(size(I)); zeros(size(I)), I] * convertInterfaceTinS(tInt) * [I, zeros(size(I)); zeros(size(I)), XL];

           % append sLayer to the overall S-Matrix ...
           SMatrix = addSmatrix(SMatrix, sLayer, I);

           tInt = [eigE, eigE; k*eigH, -k*eigH];
        end

        % add only the eigenvalue matrix tInt of the last layer to the current
        % S-matrix (so we end up at the right side of the last interface).
        % because we started the matrix also on the left side of the left
        % interface we have now the correct matrix which describes exactly one
        % layersystem or "unit cell"...
        SMatrix = addSmatrix(SMatrix, convertInterfaceTinS(tInt), I);
        
        switch (rep)
            case 1
                % nothing to do            
            case 2
                % repeat the unit cell to the proper size.
                SMatrix = repeatSMatrix(SMatrix, repSystemLocal, I);
            case 3
                % nothing to do
            otherwise
                error('Unknown error!');
        end
        
        % append results to the final S-Matrix.
        SMatrixTotal = addSmatrix(SMatrixTotal, SMatrix, I);
    end
    
%endfunction



%% ************************************************************************
%%  repeatSMatrix,
%%       repeat a S-Matrix N times.
%%       The idea here is now, do not repeat the S-Matrix by a simple loop
%%       but, decompose the number of repeatings into a sum over powers
%%       regarding to a special base (i.e. 10). Then if we want to repeat
%%       the S-Matrix i.e. 263 times we decompose it into 1*10^2 + 2*10^1 +
%%       3*10^0. Then the strategy is to biult up an SMatrix for the single
%%       structure (and repeat it 3 times), a second S-matrix for a
%%       10-superlayer structure (and repeat it 6 times) and a
%%       100-superlayer-SMatrix (and repeat it 2 times). The advantage this
%%       technique is the reduction of S-matrix concatenation operations
%%       which take all lot of calculation time due to the involved matrix
%%       conversions. 
%% ************************************************************************
%% 
function SFinal = repeatSMatrix(SMatrix, N, I)

    % choose an appropriate numerical base ... 
    if (N <= 500)
        base = 5;
    else
        base = 10;
    end

    % if the repeatings are too less, then perform a simple loop to repeat
    % the S-Matrices ...
    if (N <= 2*base)
        SFinal = repeatSMatrixInternal(SMatrix, N, I);
        return;
    end

    % decompose the number of repeatings into the proper numerical base ...
    % the value returned by this function is an array, whereby the first
    % element is the multiple of base^0, the second one that of base^1 and
    % so on. I.e. 263 --> [3,6,2].
    powerArray = powerDecomposition(N, base);

    % intitialize the final matrix to unity ...
    SFinal = eye(size(SMatrix));
    
    % SOfPower stores the S-matrices for the given number of
    % superlayers (i.e. for base=10: 1, 10, 100, ... and for base=5: 1, 5,
    % 25, ...) ...
    SOfPower = cell(length(powerArray), 1);
    
    % The lowest order always corresponds to initial S-matrix ...
    SOfPower{1} = SMatrix;
        
    % loop over the power decomposition ...
    for n = 1 : length(powerArray)-1

        % in the first step calc: current multiple times current power ...
        dummy = repeatSMatrixInternal(SOfPower{n}, powerArray(n), I);
        
        % add it to SFinal ...
        SFinal = addSmatrix(SFinal, dummy, I);
        
        % here: the next {n+1} superlayer S-matrix is generated out of the
        % ... previous ones ...
        SOfPower{n+1} = dummy;
        clear dummy;
        for c = powerArray(n)+1 : base
            SOfPower{n+1} = addSmatrix(SOfPower{n+1}, SOfPower{n}, I);
        end
    end
    
    % don't forget to add the last one ...
    SFinal = addSmatrix(SFinal, repeatSMatrixInternal(SOfPower{end}, powerArray(end), I), I);  
    
%endfunction


%% ************************************************************************
%%  logBase,
%%       log regarding to given base ...
%% ************************************************************************
%% 
function out = logBase(number, base)
    out = log10(number) ./ log10(base);
%endfunction



%% ************************************************************************
%%  powerDecomposition,
%%       decompose a given numerical value into its representation
%%       regarding to a given numerical base ...
%%       i.e. 163 --> 1*10^2 + 6*10^1 + 3*10^0 or
%%            163 --> 1*5^3 + 1*5^2 + 2*5^1 + 3*5^0
%% ************************************************************************
%% 
function [powerArray, basePower] = powerDecomposition(number, base)
    basePower = fix(logBase(number, base));
    powerArray = fix(number / base^basePower);
    number = number - (powerArray * base^basePower);
    if (number > 0)
        [lowerPowerArray, previousBasePower] = powerDecomposition(number, base);
        lowerPowerArray = [lowerPowerArray; zeros(basePower-previousBasePower-1, 1)];    
        powerArray = [lowerPowerArray; powerArray];
    else
        powerArray = [zeros(basePower, 1); powerArray];
    end
%endfunction



%% ************************************************************************
%%  repeatSMatrixInternal,
%%       simply repeat a S-Matrix N times by a loop ...
%% ************************************************************************
%% 
function SNew = repeatSMatrixInternal(SMatrix, N, I)

    if (N <= 0)
        SNew = eye(size(SMatrix));
    else
        SNew = SMatrix;
        for n = 2 : N
            SNew = addSmatrix(SNew, SMatrix, I);
        end
    end
%endfunction


%% ************************************************************************
%%  getFieldAndStructureDataForSpecificLayer,
%%       extracts the field data for a specific layer out of the structures
%%       eigVecs_E, eigVecs_H, eigVals, zjTotal.
%% ************************************************************************
%% 
function [eigE, eigH, thickness, eigval] = getFieldAndStructureDataForSpecificLayer(layer, JN, JInFront, JBehind, repeatSystem, eigVecs_E, eigVecs_H, eigVals, thicknessTotal)


    if (layer <= JInFront)
        % do nothing;
    elseif (layer > (JInFront + JN*repeatSystem))
        layer = layer - JN * (repeatSystem-1);
    else
        layer = layer - JInFront;
        layer = rem(layer-1, JN)+1;
        layer = layer + JInFront;
    end
    
    eigE = eigVecs_E{layer};
    eigH = eigVecs_H{layer};
    thickness = thicknessTotal(layer);
    eigval = eigVals{layer};
    
%endfunction


%% ************************************************************************
%%  initAdjacentMediaMatrices,
%%       initialize the "mode matrices" of the adjacent media
%% ************************************************************************
%% 
function [Incident_Vec, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, TransmittedMode_E, TransmittedMode_H] = ...
    initAdjacentMediaMatrices(eigenModeUsage, IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H, ...
    TransmittedMode_E, TransmittedMode_H, Incident_Vec, mlkm, nlkm, zindmn, ux, uy, uz, rmn, tmn, alpha0, beta0, alphamn, betamn, n1, n3)

    switch (eigenModeUsage)    

        % mode 0: ordinary FMM. Here input and output side are homogeneous, so
        % that the according plane wave modes are used.
        case 0
%             disp('eigenModeUsage == 0:');
%             disp('===========================');
%             disp('  The eigen modes of all adjacent media are calculated automatically! ');
%             disp('  Consequently all inputs like: TransmittedMode_E, TransmittedMode_H');
%             disp('  IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H are ignored!');
%             disp(' ');

            Incident_Vec = cell(size(ux));
            for j = 1 : length(Incident_Vec)
                Incident_Vec{j} = zeros(2*mlkm*nlkm, 1);
                Incident_Vec{j}(zindmn) = ux{j};
                Incident_Vec{j}(mlkm*nlkm + zindmn) = uy{j};
            end
            
            IncidentMode_E = eye(2*mlkm*nlkm);
            % IncidentMode_E = zeros(2*mlkm*nlkm);
            % IncidentMode_E(zindmn, zindmn) = ux;
            % IncidentMode_E(mlkm*nlkm + zindmn, mlkm*nlkm + zindmn) = uy;

            IncidentMode_H = [-diag(betamn.*alphamn./rmn(:)),...
                   diag(-betamn.^2./rmn(:)-rmn(:));...
                   diag(rmn(:)+alphamn.^2./rmn(:)),...
                   diag(alphamn.*betamn./rmn(:))];
            
            % IncidentMode_H = zeros(2*mlkm*nlkm);
            % IncidentMode_H(zindmn, zindmn) = sign(n1) * (beta0*uz - rmn(zindmn)*uy);
            % IncidentMode_H(mlkm*nlkm + zindmn, mlkm*nlkm + zindmn) = sign(n1) * (rmn(zindmn)*ux - alpha0*uz);

            ReflectedMode_E = eye(2*mlkm*nlkm);
            ReflectedMode_H = sign(n1) * [diag(betamn.*alphamn./rmn(:)),...
                   diag(betamn.^2./rmn(:)+rmn(:));...
                   diag(-rmn(:)-alphamn.^2./rmn(:)),...
                   diag(-alphamn.*betamn./rmn(:))];

            TransmittedMode_E = eye(2*mlkm*nlkm);
            TransmittedMode_H = sign(n3) * [diag(-betamn.*alphamn./tmn(:)),...
                   diag(-betamn.^2./tmn(:)-tmn(:));...
                   diag(tmn(:)+alphamn.^2./tmn(:)),...
                   diag(alphamn.*betamn./tmn(:))];

        % mode 1: Here only the input side is homogeneous, so
        % that the according plane wave modes are used. The modes for the
        % transmitted side directly taken from the delivered data structure ...
        case 1
            
%             disp('eigenModeUsage == 1:');
%             disp('===========================');
%             disp('  The eigen modes of the medium in the output region are taken ');
%             disp('  from user data (TransmittedMode_E, TransmittedMode_H)! ');
%             disp('  Consequently all inputs like: IncidentMode_E, IncidentMode_H, ');
%             disp('  ReflectedMode_E, ReflectedMode_H and N(2) are ignored!'); 
%             disp(' ');

            Incident_Vec = cell(size(ux));
            for j = 1 : length(Incident_Vec)
                Incident_Vec{j} = zeros(2*mlkm*nlkm, 1);
                Incident_Vec{j}(zindmn) = ux{j};
                Incident_Vec{j}(mlkm*nlkm + zindmn) = uy{j};
            end
            
            IncidentMode_E = eye(2*mlkm*nlkm);
            % IncidentMode_E = zeros(2*mlkm*nlkm);
            % IncidentMode_E(zindmn, zindmn) = ux;
            % IncidentMode_E(mlkm*nlkm + zindmn, mlkm*nlkm + zindmn) = uy;
           
            IncidentMode_H = [-diag(betamn.*alphamn./rmn(:)),...
                   diag(-betamn.^2./rmn(:)-rmn(:));...
                   diag(rmn(:)+alphamn.^2./rmn(:)),...
                   diag(alphamn.*betamn./rmn(:))];
            
            % IncidentMode_H = zeros(2*mlkm*nlkm);
            % IncidentMode_H(zindmn, zindmn) = sign(n1) * (beta0*uz - rmn(zindmn)*uy);
            % IncidentMode_H(mlkm*nlkm + zindmn, mlkm*nlkm + zindmn) = sign(n1) * (rmn(zindmn)*ux - alpha0*uz);

            ReflectedMode_E = eye(2*mlkm*nlkm);
            ReflectedMode_H = sign(n1) * [diag(betamn.*alphamn./rmn(:)),...
                   diag(betamn.^2./rmn(:)+rmn(:));...
                   diag(-rmn(:)-alphamn.^2./rmn(:)),...
                   diag(-alphamn.*betamn./rmn(:))];

        % mode 2: Here only the output side is homogeneous, so
        % that the according plane wave modes are used. The modes for the
        % reflected side (incident + refl) are directly taken from the
        % delivered data structure ... 
        case 2
%             disp('eigenModeUsage == 2:');
%             disp('===========================');
%             disp('  The eigen modes of the medium in the input region are taken ');
%             disp('  from user data (IncidentMode_E, IncidentMode_H, ReflectedMode_E, ReflectedMode_H)! ');
%             disp('  Consequently all inputs like: TransmittedMode_E, TransmittedMode_H and ');
%             disp('  N(2), ux, uy, uz are ignored!'); 
%             disp(' ');

            %%% Incident_Vec = {ones([size(IncidentMode_E, 2), 1])};
            if ~iscell(Incident_Vec)
                Incident_Vec = {Incident_Vec};
            end
            
            for iv = 1 : length(Incident_Vec)
                if (length(Incident_Vec{iv}) ~= size(IncidentMode_E, 2))
                    error('<Incident_Vec> and <IncidentMode_E> sizes do not match!');
                end
            end
                
            
            TransmittedMode_E = eye(2*mlkm*nlkm);
            TransmittedMode_H = sign(n3) * [diag(-betamn.*alphamn./tmn(:)),...
                   diag(-betamn.^2./tmn(:)-tmn(:));...
                   diag(tmn(:)+alphamn.^2./tmn(:)),...
                   diag(alphamn.*betamn./tmn(:))];   

        % mode 3: incident and transmitted side are abitrary and the modes are
        % taken out of the delivered data structure.
        case 3
            disp('eigenModeUsage == 3:');
            disp('===========================');
            disp('  The eigen modes of the medium in the input AND output region are taken ');
            disp('  from user data! Consequently all inputs like: N(1), N(2), ux, uy, uz are ignored!'); 
            disp(' ');
            
            %%% Incident_Vec = {ones([size(IncidentMode_E, 2), 1])};
            if ~iscell(Incident_Vec)
                Incident_Vec = {Incident_Vec};
            end
            
            for iv = 1 : length(Incident_Vec)
                if (length(Incident_Vec{iv}) ~= size(IncidentMode_E, 2))
                    error('<Incident_Vec> and <IncidentMode_E> sizes do not match!');
                end
            end
            
            
            % nothing is to be done.
        otherwise
            error('Unknown adjacent media type usage!');
    end

%endfunction



%% ********************************************************************
%%  INTERNAL FUNCTIONS.
%% ********************************************************************
%%

% Additional functions for construction of epsilon matrices.
function [E1,E2,E3]=internal_fcg2(STRUC,M,N)
% function [E1,E2,E3]=internal_fcg2(STRUC,M,N)
%
% Calculates three epsilon matrices required by Li's reformulation.
%

    global useNumericalFFTAlgorithm;

    mlkm=2*M+1;
    nlkm=2*N+1;
    E=internal_fcg2s4(STRUC.^2,2*M,2*N);
    [mmax,nmax]=size(STRUC);

    zmeps=2*M+1;
    zneps=2*N+1;

    eps=1./(STRUC.^2);
    
    RXY = zeros((2*M+1)^2, nmax);
    for xyind=1:nmax,
       feps_rivi=internal_fcg2s4(eps(:,xyind),2*M,0);
       apu1=flipud(feps_rivi(zmeps-2*M:zmeps));
       apu2=feps_rivi(zmeps:zmeps+2*M);
       apu=toeplitz(apu2,apu1);
       apu=apu\eye(size(apu));   
       RXY(:,xyind)=apu(:);
    end
    
    RYX = zeros((2*N+1)^2, mmax);
    for xyind=1:mmax,
       feps_rivi=internal_fcg2s4(eps(xyind,:),0,2*N);
       apu1=fliplr(feps_rivi(zneps-2*N:zneps));
       apu2=feps_rivi(zneps:zneps+2*N);
       apu=toeplitz(apu2,apu1);
       apu=apu\eye(size(apu));
       RYX(:,xyind)=apu(:);
    end
    
    if (useNumericalFFTAlgorithm)
        % if the numerical FFT should be used, ... then the remaining
        % fourier transform is now performed on the RXY and RYX data
        % structures ... The functions <internal_fcg2s2> and
        % <internal_fcg2s3> will later deliver "only" the right element of
        % this fields. They do not perform any tranformation at all in this
        % operation mode. See the implementation of these functions for
        % detailes.
        dim = 2;
        
        si = size(RXY);
        RXY = FFTC(RXY, ones(size(dim)), dim) / sqrt(prod(si(dim)));
        newSize = size(RXY);
        newSize(dim) = 2*(2*N)+1;
        RXY = EmbedExtractCentral2D(RXY, newSize);
        
        si = size(RYX);
        RYX = FFTC(RYX, ones(size(dim)), dim) / sqrt(prod(si(dim)));
        newSize = size(RYX);
        newSize(dim) = 2*(2*M)+1;
        RYX = EmbedExtractCentral2D(RYX, newSize);
    end

    for n1=-N:N,
       for n2=-N:N,
          apu1=flipud(E(zmeps-2*M:zmeps,zneps+(n1-n2)));
          apu2=E(zmeps:zmeps+2*M,zneps+(n1-n2));
          E1((n1+N)*mlkm+1:(n1+N+1)*mlkm,(n2+N)*mlkm+1:(n2+N+1)*mlkm)=toeplitz(apu2,apu1);      
          E2((n1+N)*mlkm+1:(n1+N+1)*mlkm,(n2+N)*mlkm+1:(n2+N+1)*mlkm)=internal_fcg2s2(RXY,n1,n2,M,N);   
       end
    end
    for m1=-M:M,
       for m2=-M:M,
          E3(m1+M+1:2*M+1:mlkm*nlkm,m2+M+1:2*M+1:mlkm*nlkm)=internal_fcg2s3(RYX,m1,m2,M,N);
       end
    end
    E1=E1\eye(size(E1));
% endfunction


function Y=internal_sinc(X)
% function Y=internal_sinc(X)
%
% Calculates sinc-function as used in the main program
%
    Y=(X~=0).*(X)+(X==0);
    Y=sin(X)./Y;
    Y=(X~=0).*Y+(X==0);
    %Y=H;
% endfunction


function E2=internal_fcg2s2(RXY,n1,n2,M,N)
% sub-function for internal_fcg2
    
    global useNumericalFFTAlgorithm;
    
    if (~useNumericalFFTAlgorithm) 
        E2=zeros(2*M+1);
        [mmax,nmax]=size(RXY);
        for yind=1:nmax, 
           E2=E2+reshape(RXY(:,yind),2*M+1,2*M+1)*(1/nmax*exp(-i*2*pi*(n1-n2)*(yind-1)/nmax)*internal_sinc(pi*(n1-n2)/nmax)*exp(-i*pi*(n1-n2)/nmax));
        end
    else
        E2 = reshape(RXY(:,(n1-n2) + 2*N+1),2*M+1,2*M+1);
    end
        
% endfunction


function E3=internal_fcg2s3(RYX,m1,m2,M,N)
% sub-function for internal_fcg2
    E3=zeros(2*N+1);
    global useNumericalFFTAlgorithm;
    
    if (~useNumericalFFTAlgorithm)
        E3=zeros(2*N+1);
        [nmax,mmax]=size(RYX);
        %zmeps=2*N+1;
        for xind=1:mmax,
           E3=E3+(1/mmax*exp(-i*2*pi*(m1-m2)*(xind-1)/mmax)*internal_sinc(pi*(m1-m2)/mmax)*exp(-i*pi*(m1-m2)/mmax))*reshape(RYX(:,xind),2*N+1,2*N+1);
        end
    else
        E3 = reshape(RYX(:,(m1-m2) + 2*M+1),2*N+1,2*N+1);
    end
% endfunction


function T=internal_fcg2s4(C,M,N)
% sub-function for internal_fcg2

    global useNumericalFFTAlgorithm;

    if (~useNumericalFFTAlgorithm)
        C=internal_fcg2s5(C);
        dx=C(1,1);
        dy=C(1,2);
        T=zeros(2*M+1,2*N+1);
        for u=-M:M,
           for v=-N:N,      
              T(u+M+1,v+N+1)=T(u+M+1,v+N+1)+sum( (C(:,5).*(C(:,1)-C(:,3)).*(C(:,2)-C(:,4))/(dx*dy)).*internal_sinc(pi*u*(C(:,1)-C(:,3))/dx).*internal_sinc(pi*v*(C(:,2)-C(:,4))/dy).*exp(-i*pi*u*(C(:,1)+C(:,3))/dx).*exp(-i*pi*v*(C(:,2)+C(:,4))/dy) );
           end
        end
    else
        si = size(C);
        dim = [1,2];
        
        T = FFTC(C, ones(size(dim)), dim) / sqrt(prod(si(dim)));

        newSize = size(T);
        newSize(dim) = 2*[M,N]+1;
        T = EmbedExtractCentral2D(T, newSize);
    end
% endfunction


function C=internal_fcg2s5(NS,Tdir)
% sub-function for internal_fcg2
    
    dx=1;dy=1;
    if (nargin<2),
       Tdir=0;
    end
    flipped=false;
    if Tdir==1,NS=rot90(NS);end
    [nx,ny]=size(NS);
    if (nx==1)&&(Tdir==0),
       NS=NS(:);
       flipped=true;
       [nx,ny]=size(NS);
    end
    px=dx/nx;
    py=dy/ny;
    nmin=min(min(NS));
    C(1,:)=[dx dy 0 0 nmin];
    NS=NS-nmin;
    counter=2;
    while ~(isempty(find(NS~=0))),
       n0=min(NS(NS~=0));
       apu=find(NS==n0);
       sind=apu(1);
       for t=1:length(apu)-1,      
          apu2=( fix(apu(t)/nx)+(~(fix(apu(t)/nx)==(apu(t)/nx))))~=( fix(apu(t+1)/nx)+(~(fix(apu(t+1)/nx)==(apu(t+1)/nx))));
          if apu2||((apu(t)+1)~=apu(t+1)),  
             lind=apu(t);         
             by=fix(lind/nx)+(~(fix(lind/nx)==(lind/nx)));
             bx=lind-(by-1)*nx;
             ey=fix(sind/nx)+(~(fix(sind/nx)==(sind/nx)))-1;
             ex=sind-ey*nx-1;
             C(counter,:)=[bx*px by*py ex*px ey*py n0];
             counter=counter+1;
             sind=apu(t+1);
          end
       end
       lind=apu(length(apu));
       by=fix(lind/nx)+(~(fix(lind/nx)==(lind/nx)));
       bx=lind-(by-1)*nx;
       ey=fix(sind/nx)+(~(fix(sind/nx)==(sind/nx)))-1;
       ex=sind-ey*nx-1;
       C(counter,:)=[bx*px by*py ex*px ey*py n0];
       NS(NS==n0)=0;
    end
    if (flipped)
       apu=C(:,1:2:3);
       C(:,1:2:3)=C(:,2:2:4);
       C(:,2:2:4)=apu;
    end
% endfunction
    

% this function is for conversion of a T-matrix into a S-Matrix.
% Here the relation described in the Lifeng Li paper,
% "Formulation and comparison of two recursive matrix algorithms for
% modeling layered diffraction gratings", JOSAA (1996), is used

% NEW: This function contains to branches of code: First: for the normal
% case that the T matrix to convert is quadratic and the it can be
% separated in 4 block matrices which have the same size  and which are
% quadratic, too. The second branch is the general case, where the T-matrix
% is no longer partitioned in four equal block matrices. Now the
% partitioning is defined by the midPointT argument, which gives the
% coordinate (row, col) where the four block matrices of T intersect. 
% Therefore also the outcoming S-matrix is no longer partioned into euqal
% blocks. The block size of S is connected to that of T and the the correct
% partition parameter <midPointS> is given back.
function [S, midPointS] = convertInterfaceTinS(T, midPointT)

    if (nargin < 2)
        % if only one argument is given ...
        % it is expected that the T matrix is a quadratic one which has its
        % midpoint at size(T)/2.
        midPointT = size(T)/2;
    end

    if ((size(T, 1) == size(T, 2)) && isequal(midPointT, size(T)/2))
        si = midPointT(1);

        t11 = T(1:si, 1:si);
        t12 = T(1:si, si+1 : 2*si);
        t21 = T(si+1 : 2*si, 1:si);
        t22inv = inv(T(si+1 : 2*si, si+1 : 2*si));

        S = [t11 - t12*t22inv*t21, t12*t22inv; -t22inv*t21, t22inv];    
        midPointS = midPointT;
    else
        si1 = midPointT(1);
        si2 = midPointT(2);

        t11 = T(1:si1, 1:si2);
        t12 = T(1:si1, si2+1 : end);
        t21 = T(si1+1 : end, 1:si2);
        t22inv = pinv(T(si1+1 : end, si2+1 : end));

        S = [t11 - t12*t22inv*t21, t12*t22inv; -t22inv*t21, t22inv];
        
        % achtung: S has not the same size compared to T, so the
        % partitioning is another one. Nethertheless, the midpoint is the
        % same as T, becuase the s11-block has the same size compared to
        % the t11-block.
        midPointS = midPointT;
    end
% endfunction    


% this function appends a new S-Matrix (sLayer) to an former one (SOld).
% Here the recursive relation described in the Lifeng Li paper,
% "Formulation and comparison of two recursive matrix algorithms for
% modeling layered diffraction gratings", JOSAA (1996), is used
function [SNew, midPointSNew] = addSmatrix(SOld, sLayer, I, midPointSOld, midPointSLayer)

    if (nargin < 4)
        % if only one argument is given ...
        % it is expected that the S matrices to connect are quadratic which
        % have their midpoints at size/2.
        midPointSOld = size(SOld)/2;
        midPointSLayer = size(sLayer)/2;
    end

    si1 = midPointSOld(1);
    si2 = midPointSOld(2);
    
    Tuu = SOld(1:si1, 1:si2);
    Rud = SOld(1:si1, si2+1 : end);
    Rdu = SOld(si1+1 : end, 1:si2);
    Tdd = SOld(si1+1 : end, si2+1 : end);
    
    si1 = midPointSLayer(1);
    si2 = midPointSLayer(2);

    tuu = sLayer(1:si1, 1:si2);
    rud = sLayer(1:si1, si2+1 : end);
    rdu = sLayer(si1+1 : end, 1:si2);
    tdd = sLayer(si1+1 : end, si2+1 : end);
    
    TuuNew = tuu * inv(I - Rud*rdu) * Tuu;
    RudNew = rud + (tuu * Rud * inv(I - rdu*Rud) * tdd);
    RduNew = Rdu + (Tdd * rdu * inv(I - Rud*rdu) * Tuu);
    TddNew = Tdd * inv(I - rdu*Rud) * tdd;
    
    SNew = [TuuNew, RudNew; RduNew, TddNew];
    midPointSNew = size(TuuNew);
%endfunction


% this function looks for identical epsilon distributions in the layer
% structure ...
function idNr = checkForIdenticalLayer(S, layerNr, lookDownward, useTemplatesInLayer)

    if (useTemplatesInLayer{layerNr})
        % if TRUE all layer data (eigenvectors, eigenvalues) are already
        % given, so that the idNr is always set to the layerNr to search
        % for ...
        idNr = layerNr;
    else                
        if (lookDownward)
            % look all previous elements in the structure and compare to
            % the current one ...

            c = layerNr-1;   
            idNr = -1;
            while ((idNr == -1) && (c > 0))
                if (isequal(S{c}, S{layerNr}))
                    idNr = c;
                end

                c = c - 1;
            end
        else
            % look all following elements in the structure and compare to
            % the current one ...
            
            c = layerNr+1;   
            idNr = -1;
            while ((idNr == -1) && (c <= length(S)))
                if (isequal(S{c}, S{layerNr}))
                    idNr = c;
                end

                c = c + 1;
            end        
        end
    end

%endfunction


function [M, gamma] = sortEigenStatesUpDown(M, gamma)
   
    % use sort to sort "ascending" and then flipdim to get a "descending"
    % order. This manual flipping is done due to combatibility reasons with
    % ols matlab versions which do not support sort(..., 'descend').
    [dummy, index] = sort(real(gamma) + imag(gamma));
    dim = (size(index,1) == 1) + 1;
    index = flipdim(index, dim);

    gamma = gamma(index);    
    M = M(:, index);
    
%endfunction