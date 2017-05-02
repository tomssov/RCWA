function val = AZO(nmlambda)
%lam0 = wavelength (in nm)
%val = refractive index for wavelength lam0 (experimental data)

data = load('AZO.mat');

nmlambdas=data.data(:,1);
AZO_n=data.data(:,2);
AZO_k=data.data(:,3);
val=interp1(nmlambdas,AZO_n+1i.*AZO_k,nmlambda);  % linear interpolation into the table