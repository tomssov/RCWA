function val = Silver(varargin)
%lam0 = wavelength (in nm)
%val = refractive index for wavelength lam0 (experimental data)

% Setup data structure and default values
loc = struct(...
    'nmlambda', 400 ...      % Incident light loc nm
);

loc = ApplyVarargin(loc, varargin);


%data = csvread('silvernk1200.csv');
%data = xlsread('epsAg.xlsx');
data = load('silver.mat');


nmlambda=data.data(:,1);
ag_n=data.data(:,2);
ag_k=data.data(:,3);
val=interp1(nmlambda,ag_n+1i.*ag_k,loc.nmlambda);  % linear interpolation into the table


