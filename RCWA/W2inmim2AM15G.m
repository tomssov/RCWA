function val = W2inmim2AM15G(varargin)
% Input: wavelength (in nm)
%val = refractive index for wavelength lam0 (experimental data)

% Setup data structure and default values
loc = struct(...
    'nmlambda', 400 ...      % Incident light loc nm
);

loc = ApplyVarargin(loc, varargin);

data = csvread('AM15G.csv');

nmlambda = data(:,1);
flux = data(:,2);
val = interp1(nmlambda,flux,loc.nmlambda);  % linear interpolation into the table
