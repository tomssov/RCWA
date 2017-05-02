function output = aSiHGC(varargin)
% GC Cody Band Edge Function
%  Input       Description
%  nmlambda    Incident light nm
%  eVEg        a-Si:H Bandgap profile
%  optimum     Is material optimal? 0 no or 1 yes (use 0 for thesis)
%
% Output
% - a matrix of permittivities with rows at set bandgap, columns at set
% wavelength, e.g. output(:,1) gives the permittivities at the first input
% wavelength.

%  Permittivity at infinity
epsinf = 1;
optimum = 0;

if length(varargin)==2
   nmlambda = varargin{1};
   eVEg = varargin{2};
elseif length(varargin)==3
   nmlambda = varargin{1};
   eVEg = varargin{2};
   optimum = varargin{3};
elseif length(varargin)==4
   nmlambda = varargin{1};
   eVEg = varargin{2};
   optimum = varargin{3};
   epsinf = varargin{4};
else
    error('aSiHGC takes either 2, 3, or 4 arguments');
end


if(size(eVEg,1)>size(eVEg,2))
    eVEg=eVEg.';
end

eVEg_unique = unique(eVEg); % Only calculate for unique Eg then rebuild

output = zeros(length(eVEg), length(nmlambda));

for ll = 1:length(nmlambda)
eVlight = eV_from_nm(nmlambda(ll));


eVEu = Eu(eVEg_unique,optimum);
eVEt = Et(eVEg_unique);

E1 = eVEt.*LE(eVEt, eVEg_unique, optimum).*GC(eVEt, eVEg_unique, optimum);


    IU = E1.*(exp((eVlight - eVEt)./eVEu).*(ei((eVEt - eVlight)./eVEu) -...
        ei(-eVlight./eVEu)) - exp(-(eVlight + eVEt)./eVEu).*(ei((eVEt + eVlight)./eVEu) -...
        ei(eVlight./eVEu)))./(pi.*eVlight);
    
    
    eps1 = epsinf + IU + ICL(eVlight, eVEg_unique, optimum);
    
    eps2 = (eVlight <= eVEt).*(E1./eVlight).*exp((eVlight - eVEt)./eVEu);
    eps2 = eps2 + (eVlight > eVEt).*GC(eVlight, eVEg_unique, optimum).*LE(eVlight, eVEg_unique, optimum);
    
    
    for i=1:length(eVEg_unique)
        output(eVEg == eVEg_unique(i), ll) = eps1(i) + 1i.*eps2(i);
    end
end

end
