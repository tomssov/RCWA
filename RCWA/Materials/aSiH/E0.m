function out = E0(Eg,opt)
% Lorentz Oscillator Energy E0
% Eg:   a-Si:H bandgap of Eg (eV/V)
% opt: Is material opt (0 for no, 1 for yes)

idx = sign(Eg-1.803);
if(opt==1)
    out(idx==-1) = 3.832 + 0.2914.*(Eg(idx==-1)-1.803);
else
    out(idx==-1) = 3.832 + 0.0155.*(Eg(idx==-1)-1.803);
end
out(idx==1) = 3.832 - 0.9354.*(Eg(idx==1)-1.803);
end

