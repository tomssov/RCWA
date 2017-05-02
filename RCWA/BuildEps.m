function [eps,loc]=BuildEps(mat_cat, nmx, nmz, nmlambda, Eg, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc = DefaultLoc;
loc = ApplyVarargin(loc,varargin);


if(loc.dmat==1) % Dielectric material in the grating
    glasseps = (glass('nmlambda', nmlambda)+(10^-3)*1i).^2;
    loc.epsd = glasseps;     % permittivity of the glass
elseif(loc.dmat==2)
    glasseps = AZO(nmlambda);
    loc.epsd = glasseps;
end

if(loc.wmat==1) % Dielectric material in the window
    glasseps = (glass('nmlambda', nmlambda)+(10^-3)*1i).^2;
    loc.epsw = glasseps;     % permittivity of the glass
elseif(loc.dmat==2)
    glasseps = AZO(nmlambda);
    loc.epsw = glasseps;
end

if(loc.gmat == 1) % Metallic material in the grating
    silvereps = (Silver('nmlambda', nmlambda).^2);
    loc.epsm = silvereps;    % permittivity of the circle
end

% Pass the material function, argument of position here, to the solver.
% Solved once nmz is populated in BuildEps
if(loc.material == 1 && loc.Nz ~=0)
    loc.epsJ = aSiHGC(nmlambda, Eg(mat_cat==2)).';
elseif(loc.material ==2)
    loc.epsJ = Faryad(nmz, nmlambda);
end


% Build elist
% Begin building permittivity (eps) vectors here
Nz = length(nmz);
Nx = length(nmx);
eps=zeros(Nx, Nz);


eps(:, mat_cat == 4) = loc.nsa;

eps(:, mat_cat == 3) = loc.epsW;

if length(loc.epsJ)==1
    eps(:, mat_cat == 2) = loc.epsJ;
else
    eps(:, mat_cat == 2) = kron(ones(length(nmx),1),loc.epsJ);
end

eps(:, mat_cat == 0) = loc.epsm;

zgrating= nmz(mat_cat==1);
epsgrating = zeros(length(nmx), length(zgrating));
for i=1:sum(mat_cat==1)
    epsgrating(:,i) = gepsz(nmx, zgrating(i), loc);
end
eps(:, mat_cat == 1) = epsgrating;