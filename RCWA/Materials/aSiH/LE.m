function output = LE(eVlight, eVEg, optimum)
% GC Cody Band Edge Function
%  Takes the final input if light is specified multiple ways
%
%  Input String    Default   Description
%  'nmlambda'      400       Incident light in nm
%  'mlambda'       400*10^-9 Incident light in m
%  'eVlight'       1.1       Incident light in eV
%  'eVEg'          1.1       a-Si:H Bandgap
%  'optimum'       0         Is material optimal 0 or 1

% Calculations
E0val = E0(eVEg, optimum);
Gammaval = Gamma( eVEg, optimum);
Aval = A( eVEg, optimum);

output = Aval.*E0val.*Gammaval.* eVlight./(( eVlight.^2-E0val.^2).^2+Gammaval.^2.* eVlight.^2);

end

