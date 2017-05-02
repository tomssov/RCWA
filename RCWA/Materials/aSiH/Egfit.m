function out = Egfit(Eg)
% Bandgap fit parameter Egfit
% Eg:   a-Si:H bandgap of Eg (eV/V)

out = 1.727 + 0.8153*(Eg - 1.803);

end

