function out = Eu(Eg,opt)
% Urbach Tail parameter Eu
% Eg:   a-Si:H bandgap of Eg (eV/V)
% opt: Is material opt (0 for no, 1 for yes)

idx = sign(Eg-1.803);
if(opt==1)
    out(idx==-1) = 49.03 - 4.866.*(Eg(idx==-1)-1.803);
else
    out(idx==-1) = 49.03 - 28.31.*(Eg(idx==-1)-1.803);
end
out(idx==1) = 49.03 + 90.63.*(Eg(idx==1)-1.803);
end

