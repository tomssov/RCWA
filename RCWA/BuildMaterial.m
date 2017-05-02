function [mat_cat, nmz, nmdz] = BuildMaterial(varargin)
% BUILDMATERIAL Categorises the material regions

%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc = DefaultLoc;
loc = ApplyVarargin(loc,varargin);

nmLz = loc.nmLp+loc.nmLi+loc.nmLn;

% Only allow gratings that do not contain semiconductor material
if(loc.nmda < loc.nmLg)
    loc.nmda = loc.nmLg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Create a list of positions      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(loc.Nm >= 1 && loc.nmLm > 0)
    dm = loc.nmLm/loc.Nm;
    nmz_mirror = linspace(0.5*dm, loc.nmLm - 0.5*dm , loc.Nm);
 else
    dm=0;
    nmz_mirror = [];
end

if(loc.Ng >= 1 && loc.nmda > 0)
    dg = loc.nmda/loc.Ng;
    nmz_grating = linspace(loc.nmLm  + 0.5*dg, loc.nmLm + loc.nmda - 0.5*dg, loc.Ng);
 else
    dg=0;
    nmz_grating = [];
end

if(loc.Nz >= 1 && nmLz)
    dz=nmLz/loc.Nz;
    nmz_junction = linspace((loc.nmda+loc.nmLm)+0.5*dz, (loc.nmda+loc.nmLm)+nmLz -0.5*dz, loc.Nz);
else
    dz = 0;
    nmz_junction = [];
end

if(loc.Nw >= 1 && loc.nmLw > 0)
    dw = loc.nmLw/(loc.Nw);
    nmz_window = linspace((loc.nmda+loc.nmLm)+nmLz+0.5*dw, (loc.nmda+loc.nmLm)+nmLz+loc.nmLw - 0.5*dw, loc.Nw);
 
else
    dw = 0;
    nmz_window = [];
end

if(loc.Nair >= 1 && loc.nmLair > 0)
    dair = loc.nmLair/(loc.Nair);
    nmz_air = linspace((loc.nmda+loc.nmLm)+nmLz+loc.nmLw+0.5*dair, (loc.nmda+loc.nmLm)+nmLz+loc.nmLw+loc.nmLair - 0.5*dair, loc.Nair);
else
    dair = 0;
    nmz_air = [];
end

nmz = [ ... 
    nmz_mirror,....
    nmz_grating, ...
    nmz_junction, ...
    nmz_window, ...
    nmz_air];

nmdz = [ ...
    0*nmz_mirror+dm, ...
    0*nmz_grating+dg, ...
    0*nmz_junction+dz, ...
    0*nmz_window + dw, ...
    0*nmz_air + dair ...
    ];

nmz = fliplr(nmz);
nmdz = fliplr(nmdz);



% Categorise locations
mirror = sign(nmz - (loc.nmLm));
mirror = (1/2)*(-mirror + abs(mirror));

grating = sign(nmz - (loc.nmLm+loc.nmda)) + mirror;
grating = (1/2)*(-grating + abs(grating));

junction =  sign(nmz - (nmLz+(loc.nmda+loc.nmLm))) + mirror + grating;
junction = (1/2)*(-junction + abs(junction));

window = sign(nmz - (nmLz+(loc.nmda+loc.nmLm)+loc.nmLw)) + mirror + junction + grating;
window = (1/2)*(-window + abs(window));

air = sign(nmz - (nmLz+(loc.nmda+loc.nmLm)+loc.nmLw+loc.nmLair)) + mirror + junction + grating+window;
air = (1/2)*(-air + abs(air));

mat_cat = 0*mirror + 1*grating + 2*junction + 3 * window + 4* air;

end

