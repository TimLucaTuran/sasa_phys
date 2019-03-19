function efkt = aufitek(kmin,kmax,anzahl)


% function [lambda,efkt,lambdaori,efktori]=aufitek(kmin,kmax,anzahl)


%clear all
%minlambda=300e-9;
%maxlambda=900e-9;
%anzahl=200;

knum=linspace(kmin,kmax,anzahl);
lambda=1./knum*1e-2;
c0=2.99792458e8;
eV=[0.1771 0.2066 0.248 0.31 0.4133 0.64 0.77 0.89 1.02 1.14 1.26 1.39 1.51 1.64 1.76 1.88 2.01 2.13 2.26 2.38 2.50 2.63 2.755 2.883 3.024 3.179 3.351 3.543 3.757 4 4.276];
U=eV.*1.60217733e-19;
nuspace=U./6.6260755e-34;
n=[6.26 4.70 3.27 2.04 1.17 0.92 0.56 0.43 0.35 0.27 0.22 0.17 0.16 0.14 0.13 0.14 0.21 0.29 0.43 0.62 1.04 1.31 1.4 1.607 1.641 1.671 1.706 1.751 1.813 1.83 1.750];
k=[48.15 41.7 35.20 27.9 21.0 13.78 11.21 9.519 8.145 7.150 6.350 5.663 5.083 4.542 4.103 3.697 3.272 2.863 2.455 2.081 1.833 1.849 1.88 1.9341 1.9575 1.9401 1.8833 1.8471 1.8704 1.9159 1.9044];

lambdaspace=c0./nuspace;
efktexp=(n+1i*k).^2;
for lauf1=1:anzahl
   v1=find(lambdaspace-lambda(lauf1)>0);
   v2=find(lambdaspace-lambda(lauf1)<0);
   if (isempty(v1))
       v2=1:(length(eV)-1);
       v1=length(eV);
   end;
   uvlam=lambdaspace(v2(1));
   lvlam=lambdaspace(v1(length(v1)));
   uvexp=efktexp(v2(1));
   lvexp=efktexp(v1(length(v1)));
   efkt(lauf1)=(uvexp-lvexp)/(uvlam-lvlam)*(lambda(lauf1)-lvlam)+lvexp;
end;
% lambdaori=lambdaspace;
% efktori=efktexp;




