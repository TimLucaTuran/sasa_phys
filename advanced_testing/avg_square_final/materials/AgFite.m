function efkt = AgFite(minlambda,maxlambda,anzahl)
% Anzahl... Anzahl an Stuezstellen zwischen minlambda un maxlambda

% alternativer Funktionsaufruf
% [lambda,efkt,lambdaori,efktori,efktdrude]=AgFite(minlambda,maxlambda,anzahl);

% lambda in m
lambda=linspace(minlambda,maxlambda,anzahl);
c0=2.99792458e8;

eV=[0.125 0.13 0.14 0.15 .16 .17 .18 .19 .20 .22 .24 .26 .28 .30 .32 .34 ...
    .36 .38 .40 .42 .44 .46 .48 .50 .52 .54 .56 0.64 0.77 0.89 1.02 1.14 ...
    1.26 1.39 1.51 1.64 1.76 1.88 2.01 2.13 2.26 2.38 2.50 2.63 2.75 2.88 ...
    3.00 3.12 3.25 3.37 3.50 3.62 3.74 3.87 3.99 4.12 4.24 4.36 4.49 4.61 ...
    4.74 4.86 4.98 5.11 5.23 5.36 5.48 5.60 5.73 5.85 5.98 6.10 6.22 6.35 6.47 6.60];

U=eV.*1.60217733e-19;
nuspace=U./6.6260755e-34;

n= [13.11 12.21 10.69 9.441 8.376 7.461 6.67 5.96 5.355 4.425 3.732 3.202 ...
    2.786 2.446 2.16 1.915 1.71 1.536 1.387 1.265 1.168 1.083 1.007 .939 ...
    .878 .823 .774 0.24 0.15 0.13 0.09 0.04 0.04 0.04 0.04 0.03 0.04 0.05 ...
    0.06 0.05 0.06 0.05 0.05 0.05 0.04 0.04 0.05 0.05 0.05 0.07 0.10 0.14 ...
    0.17 0.81 1.13 1.34 1.39 1.41 1.41 1.38 1.35 1.33 1.31 1.30 1.28 1.28 ...
    1.26 1.25 1.22 1.20 1.18 1.15 1.14 1.12 1.10 1.07];

k=[53.7 52.2 49.4 47.1 44.8 42.5 40.4 38.6 37 34 31.3 29 26.9 25.1 23.5 ...
    22.1 20.9 19.8 18.8 17.9 17.1 16.4 15.7 15.1 14.5 14.0 13.5 14.08 ...
    11.85 10.10 8.828 7.795 6.992 6.312 5.727 5.242 4.838 4.483 4.152 ...
    3.858 3.586 3.324 3.093 2.869 2.657 2.462 2.275 2.070 1.864 1.657 ...
    1.419 1.142 0.829 0.392 0.616 0.964 1.161 1.264 1.331 1.372 1.387 ...
    1.393 1.389 1.378 1.367 1.357 1.344 1.342 1.336 1.325 1.312 1.296 ...
    1.277 1.255 1.232 1.212];

lambdaspace=c0./nuspace;
efktexp = (n+1i*k).^2;
efkt = zeros(anzahl,1);
for lauf1=1:anzahl
    v1    = find(lambdaspace-lambda(lauf1)>0);
    v2    = find(lambdaspace-lambda(lauf1)<0);
    uvlam = lambdaspace(v1(length(v1)));
    lvlam = lambdaspace(v2(1));
    uvexp = efktexp(v1(length(v1)));
    lvexp = efktexp(v2(1));
    efkt(lauf1) = (uvexp-lvexp)/(uvlam-lvlam)*(lambda(lauf1)-lvlam)+lvexp;
end;

% ======= extra output that is not necessary for FMM simulations ======== %
% lambdaori=lambdaspace;
% efktori=efktexp;
% nudrude=c0./lambda;
% U=nudrude.*6.6260755e-34;
% omega=U./1.60217733e-19;
% efktdrude=8.926-11.585.^2./(omega.^2+1i.*0.203.*omega);
%efktdrude=1-2175e12^2./(omega.^2+4.35e12^2);

