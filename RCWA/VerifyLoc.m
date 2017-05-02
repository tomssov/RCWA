function [ loc ] = VerifyLoc( loc )
%VERIFYLOC Summary of this function goes here
%   Detailed explanation goes here


if(loc.Nair == 0)
    loc.nmLair = 0;
end

if(loc.Nw == 0)
    loc.nmLw = 0;
end
if(loc.Nz == 0)
    loc.Lp =0;
    loc.Li =0;
    loc.Ln =0;
end
if(loc.Ng == 0)
    loc.nmda = 0;
    loc.nmLg = 0;
end

if(loc.Nm == 0)
    loc.nmLm = 0;
end

if loc.nmda < loc.nmLg
    loc.nmda = loc.nmLg;
end
    
end

