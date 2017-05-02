function out = A(Eg,opt)
% Lorents Oscillator Width Gamma
% Eg:   a-Si:H bandgap of Eg (eV/V)
% opt: Is material opt (0 for no, 1 for yes)

idx = sign(Eg-1.803);
if(opt==1)
    out(idx==-1) = 74.94 + 1.505.*(Eg(idx==-1)-1.803);
else
    out(idx==-1) = 74.94 -60.7.*(Eg(idx==-1)-1.803);
end
out(idx==1) = 74.94 + 387.3.*(Eg(idx==1)-1.803);
end

