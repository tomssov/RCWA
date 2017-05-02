%function eps = Faryad(nmz, nmlambda)

% aSiGe139=  

aSiHGC('eVEg', 1.95, 'nmlambda', 400, 'optimum',1)
temp = xlsread('epsaSiC195.xlsx');
temp(1,1)
temp(1,2)
temp(1,3)

%aSiGe158=xlsread('epsaSiGe158.xlsx');
% aSiC195=xlsread('epsaSiC195.xlsx');
% aSi18=xlsread('epsaSi18.xlsx');
% nazO=xlsread('nAZO.xlsx');



eps = zeros(length(nmz),1);

