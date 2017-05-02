function out = Gamma(Eg,opt)
% Lorents Oscillator Width Gamma
% Eg:   a-Si:H bandgap of Eg (eV/V)
% opt: Is material opt (0 for no, 1 for yes)

idx = sign(Eg-1.803);
if(opt==1)
    out(idx==-1) = 2.122 - 0.9931.*(Eg(idx==-1)-1.803);
else
    out(idx==-1) = 2.122 - 2.197.*(Eg(idx==-1)-1.803);
end
out(idx==1) = 2.122 + 4.737.*(Eg(idx==1)-1.803);
end

