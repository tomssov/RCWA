function geps = gepsz(nmx, nmz, loc)
% gx Grating relief
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the grating function z = g(x) here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmgval = gx(nmx, loc);

kron(ones(length(nmx),1),nmz);

idx = sign(nmz - loc.nmLm -nmgval);

geps(idx==1) = loc.epsd;
geps(idx==0) = (loc.epsd+loc.epsm)/2;
geps(idx==-1) = loc.epsm;


