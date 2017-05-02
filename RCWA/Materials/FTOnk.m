function val = FTOnk(lam0)
%lam0 = wavelength (in nm)
%val = refractive index for wavelength lam0 (experimental data)

% if ((lam0<400)||(lam0>1100))
%     error('No data for this value of lamnda')
% end

data = load('FTOnk1200.mat');

lambda=data.data(:,1);
FTO_n=data.data(:,2);
FTO_k=data.data(:,3);
val=interp1(lambda,FTO_n+1i.*FTO_k,lam0);  % linear interpolation into the table
