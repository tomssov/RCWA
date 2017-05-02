function g = gx(nmx, loc)
% gx Grating relief
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the grating function z = g(x) here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (loc.type == 0) % Square grating    
    idx = sign(abs(nmx) - 0.5 * loc.zeta * loc.nmLx);
    g(idx==1) = 0;
    g(idx==0) = loc.nmLg;
    g(idx==-1) = loc.nmLg;

elseif (loc.type == 1) % Sinusoidal   
    idx = sign(abs(nmx) - 0.5 * loc.zeta * loc.nmLx);
    g(idx==1) = 0;
    g(idx==0) = 0;
    g(idx==-1) = loc.nmLg*sin(pi*(nmx(idx==-1) + 0.5 * loc.zeta * loc.nmLx)/(loc.zeta*loc.nmLx));
elseif (loc.type == 2) % Pyramid
    idx = sign(abs(nmx) - 0.5 * loc.zeta * loc.nmLx);
    g(idx==1) = 0;
    g(idx==0) = 0;
    m = loc.nmLg/(loc.nmLx*loc.zeta*0.5);
    g(idx==-1) = loc.nmLg - m*abs(nmx(idx==-1));    
else % Spherical - note that the back of the spheres are filled
    idx = sign(abs(nmx) - 0.5 * loc.zeta * loc.nmLx);
    g(idx==1) = 0;
    g(idx==0) = 0;
    c = 0.5*loc.zeta*loc.nmLx;        % Grating profile halfwidth
    a = 0.5*(c^2/loc.nmLg-loc.nmda);  % Distance sphere is sunk into mirror
    r = a + loc.nmLg;                    % Radius of circle
    g(idx==-1) = sqrt(r^2-nmx(idx==-1).^2) - a; % Eqn. of surface
end