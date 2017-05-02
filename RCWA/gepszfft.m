function fgeps = gepszfft(nmz,loc)
% gx Grating relief
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the grating function z = g(x) here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc.NFFT = max(loc.Nt, loc.NFFT); % NFFT is the number of FFT terms taken. Must be at least Nt

nmx = linspace(-loc.nmLx/2, loc.nmLx/2-loc.nmLx/(4*loc.NFFT+1), 4*loc.NFFT+1);

fgeps = zeros(length(nmz), 4*loc.Nt+1);

for i = 1:length(nmz)
    geps = gepsz(nmx, nmz(i), loc);
    fgeps_temp = fft(geps)/(4*loc.NFFT+1);
    fgeps(i,:) = [fgeps_temp(end-2*loc.Nt+1:end),fgeps_temp(1:2*loc.Nt+1)];    
end

fgeps(:,2:2:end) = -fgeps(:,2:2:end);