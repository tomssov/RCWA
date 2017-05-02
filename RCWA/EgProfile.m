function [ Eg ] = EgProfile(nmz, varargin)
% EgProfile produces a vector of bandgap values at each nmz location, i.e.
% Eg(nmz)
%
% print((array) nmz, (structure) loc)

loc = struct(...
    'Eg0', 1.2, ...     % Baseline Bandgap
    'A', 0, ...         % Perturbation Amplitude
    'kappa', 3, ...     % Periods
    'phi',0.75, ...     % Phase
    'alpha',0, ...      % Shaping Parameter
    'nmLi', 100, ...    % i-layer thickness
    'nmLn', 10, ...     % n-layer thickness
    'nmLp', 10 ...      % p-layer thickness
    );

loc = ApplyVarargin(loc,varargin);

Eg = zeros(1,length(nmz));

 
% Specify the bandgap profile here Eg(nmz) for the p-, i-, and n-layers
player = sign(nmz - (loc.nmLp));
player = (1/2)*(-player + abs(player));
ilayer =  sign(nmz - (loc.nmLp+loc.nmLi)) + player;
ilayer = (1/2)*(-ilayer + abs(ilayer));
nlayer = sign(nmz - (loc.nmLp+loc.nmLi+loc.nmLn)) + ilayer + player;
nlayer = (1/2)*(-nlayer + abs(nlayer));

Eg(nlayer == 1) = 1.95;
Eg(ilayer == 1) = loc.Eg0 + loc.A .* ((0.5*(sin(2* pi * (nmz(ilayer == 1) - loc.nmLn) * loc.kappa/loc.nmLi + 2*pi*loc.phi )+1)).^loc.alpha);
Eg(player == 1) = 1.95;


end

