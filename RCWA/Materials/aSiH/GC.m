function output = GC(eVlight, eVEg, optimum)
% GC Cody Band Edge Function
%  Takes the final input if light is specified multiple ways
%
%  Input String    Default   Description
%  'eVlight'       1.1       Incident light in eV
%  'eVEg'          1.1       a-Si:H Bandgap
%  'optimum'       0         Is material optimal 0 or 1

% Calculations
output = ( eVlight- eVEg).^2./(( eVlight- eVEg).^2+Ep(eVEg, optimum).^2);

end

