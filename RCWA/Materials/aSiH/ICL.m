function output = ICL(eVlight, eVEg, optimum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GC Cody Band Edge Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Takes the final input if light is specified 
% multiple ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Possible inputs:
%  Input String    Default   Description
%  'nmlambda'      400       Incident light in nm
%  'mlambda'       400.*10.^-9 Incident light in m
%  'eVlight'       1.1       Incident light in eV
%  'eVEg'          1.1       eVA-Si:H Bandgap
%  'optimum'       0         Is material optimal 0 or 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
% ICL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations begin here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eVA = A( eVEg, optimum);
eVE0 = E0( eVEg, optimum);
eVGamma = Gamma( eVEg, optimum);
eVEp = Ep( eVEg, optimum);
%eVEgfit = Egfit(in.eVEg);
eVEt = Et( eVEg);
eVEE =  eVlight;


LDE = (eVEE.^2 - eVE0.^2).^2 + eVGamma.^2.*eVEE.^2;
LDT = (eVEt.^2 - eVE0.^2).^2 + eVGamma.^2.*eVEt.^2;
I0C = (pi./2 - atan((eVEt -  eVEg)./eVEp))./eVEp;
zet = (eVE0.^2 - eVGamma.^2./2).^(1/2);
Chi = (4*eVE0.^2 - eVGamma.^2).^(1/2);
FF2 = eVEp.^2 +  eVEg.^2;
KK2 = 2*FF2 + 2*zet.^2 - 4* eVEg.^2;
YY4 = eVE0.^4 + FF2.*(KK2 - FF2) - 4* eVEg.^2.*KK2;
c0C = eVEE.*GC(eVEE,  eVEg,  optimum)./(2*LDE);
d0C = -eVEE.*(eVEE +  eVEg).^2./(2*LDE.*((eVEE +  eVEg).^2 + eVEp.^2));


b0C = YY4.*FF2.*(LDE.*((c0C - d0C)./eVEE + 2* eVEg.*KK2.*(c0C + d0C)./YY4) - 1)./((KK2 - FF2).*FF2.*YY4 + eVE0.^4.*YY4 + 4* eVEg.^2.*FF2.*KK2.^2);
b1C = (2* eVEg.*KK2.*b0C - LDE.*(c0C + d0C))./YY4;
a3C = -(b1C + c0C + d0C);
a2C = -(b0C + 2* eVEg.*b1C + eVEE.*(c0C - d0C));
a1C = -(eVEE.^2 - 2*zet.^2).*(c0C + d0C) - 2* eVEg.*b0C + (KK2 - FF2).*b1C;
a0C = 1 - eVEE.*(eVEE.^2 - 2*zet.^2).*(c0C - d0C) + (KK2 - FF2).*b0C + 2* eVEg.*KK2.*b1C;
I1T = (pi - 2*atan(2*(eVEt.^2 - zet.^2)./(Chi.*eVGamma)))./(2*Chi.*eVGamma);
I0AT = (pi - atan((2*eVEt + Chi)./eVGamma) + atan((-2*eVEt + Chi)./eVGamma))./(2*eVGamma);
I0BT = log((eVEt.^2 + eVE0.^2 + Chi.*eVEt)./(eVEt.^2 + eVE0.^2 - Chi.*eVEt))./(4*Chi);

ITL = (2*eVA.*eVE0.*eVGamma./pi).*(a3C .*(zet.^2.*I1T - log(LDT).^(1./4)) + a2C.*(I0AT + I0BT) + a1C.*I1T + a0C.*(I0AT - I0BT)./eVE0.^2 - c0C.*log(abs(eVEE - eVEt)) - d0C.*log(eVEE + eVEt));

output = ITL + (2*eVA.*eVE0.*eVGamma./pi).*(b1C.*( eVEg.*I0C - log((eVEt -  eVEg).^2 + eVEp.^2).^(1./2)) + b0C.*I0C);

end

