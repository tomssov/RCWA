function out = Ep(Eg,opt)
% Bandgap parameter Ep
% Eg:   a-Si:H bandgap of Eg (eV/V)
% opt: Is material opt (0 for no, 1 for yes)

idx = sign(Eg-1.803);
if(opt==1)
    out(idx==-1) = 1.134 + 1.001.*(Eg(idx==-1)-1.803);
else
    out(idx==-1) = 1.134 -0.3157.*(Eg(idx==-1)-1.803);
end
out(idx==1) = 1.134 + 9.731.*(Eg(idx==1)-1.803);
end

