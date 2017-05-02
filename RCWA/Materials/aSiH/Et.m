function out = Et(Eg)
% Urbach Tail Parameter Et
% Eg:   a-Si:H bandgap of Eg (eV/V)

out = 1.85 + 0.8601*(Eg - 1.803);

end

