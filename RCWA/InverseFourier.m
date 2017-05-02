function [ output ] = InverseFourier(x, fcomp, Lx, method, varargin)
% Calculate the inverse fourier of a dataset at points x
% x positions lie in range -Lx/2:Lx/2
% Nt terms of fourier sum taken

fcomp_temp =  squeeze(fcomp);
terms = length(fcomp_temp);
NT = (terms-1)/2;

if (length(varargin)==1)
    temp=floor(varargin{1});
    NT0 = min(NT,temp);
else
    NT0 = NT;
end

%method=0;
%hold off;
%for method = 0:1
if(method ==1)
    
    fcomp_temp(2:2:end) = -(fcomp_temp(2:2:end));
           
    FFTorder = (2*NT+1)*[fcomp_temp(NT+1:NT+1+NT0);fcomp_temp(NT-NT0+1:NT)];
    
    invff = ifft(FFTorder); % Correct Scaling for zeroth order
    output = interpft(invff, length(x));
    output = circshift(output,-0*length(x)/4);
  
    
else
    output=0;
    for j = -NT0:NT0
       k =  2*pi * (j) /Lx; % add incident angle
       output = output + fcomp(j+NT+1).*exp(1i*k*x);
    end
   % output = output/sqrt(NT+1);
end

%plot(x,abs(output))
%hold
%end
%legend('sum','FFT');
