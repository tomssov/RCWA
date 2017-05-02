function [Z,Tn,R,E,H] ...
    = RCWAold(inmk0, epsf, nmdz, varargin)
% [Z, Tn, R, E, H] = RCWA(inmk0, epsf, nmdz, (opt) varargin)
% Z is the Z matrix from EMSW
% Tn the tranfter matrix for each mode through each slice
% R the reflected wave
% E is the resulting electric field
% H is the resulting magnetic field

loc=DefaultLoc;
loc = ApplyVarargin(loc,varargin);
Ns = length(nmdz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompose light wave vectors (for incident, reflected and transmitted) - in air with refractive index nsa, grating
% with period nmLx
radtheta = rad_from_deg(loc.degtheta);
inmk0x = loc.nsa*inmk0*sin(radtheta)+(-loc.Nt:loc.Nt)*2*pi/loc.nmLx;
inmk0z = sqrt(loc.nsa^2*inmk0^2-inmk0x.^2);

% Create Y matrices depending on polarization (ESW - eqns. 2.128 and 2.129)
if (loc.pol == 1) % nonzero right columns of Y matrices
    Ye_inc = - diag(inmk0z/(loc.nsa*inmk0));
    Ye_ref = - Ye_inc;
    
    Yh_inc = - diag(ones(2*loc.Nt+1,1));
    Yh_ref = Yh_inc;
    
else % nonzero left columns of Y matrices
    Ye_inc = loc.nsa*diag(ones(2*loc.Nt+1,1));
    Ye_ref = Ye_inc;
    
    Yh_inc = - diag(inmk0z/(loc.nsa*inmk0));
    Yh_ref = - Yh_inc;
end

Y_inc = [Ye_inc; Yh_inc]; % Y matrix for incident and transmitted field
Y_ref = [Ye_ref; Yh_ref];

Z  = cell(Ns+1,1);
Z{Ns+1} = Y_inc; % Eqn. 2.421

% Initialise matrices
Vn = cell(Ns+1,1);
inmGn = cell(Ns+1,1);
X_Upper = cell(Ns+1,1);
inmG_Upper = cell(Ns+1,1);

NullMat=zeros(2*loc.Nt+1);

inmKX=diag(inmk0x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     The implementation of the stable algorithm begins
% (The notation is taken from Electromagnetic Surface Waves -
%  T. Mackay, J. Polo, and A. Lakhtakia, p75-77)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E=(Ex,Ey,Ez)=(Ex,0,Ez)+(0,Ey,0) and H=(Hx,Hy,Hz) = (Hx,0,Hz) +(0,Hy,0)
% Since E and H are orthogonal and by linearity, you can solve Maxwell
% for either E1=(Ex,0,Ez) and H1 = (0,Hy,0)
% or E2 = (0,Ey,0) and H2 = (Hx,0,Hz).
%(Ex,0,Ez)  and (0,Hy,0) corresponds to P-polaraization
%(0,Ey,0) and  (Hx,0,Hz) corresponds to S-polaraization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for n = Ns:-1:1
    % Select the optical permittivity, and the inverse, for the slice
    epsz = epsf(:,n);
    % epsz_recip = epsmat_recip(:,n);
    
       EPSZ=toeplitz(epsz(2*loc.Nt+1:-1:1),epsz(2*loc.Nt+1:4*loc.Nt+1));             % Eqn. 2.114
        %   EPSZ_recip=toeplitz(epsz_recip(2*loc.Nt+1:-1:1),epsz_recip(2*loc.Nt+1:4*loc.Nt+1));
        
        if (loc.pol == 1)
            inmP14 = inmk0*eye(2*loc.Nt + 1) - (1/inmk0)*inmKX*inv(EPSZ)*inmKX;
            %inmP41 = inmk0*inv_epsf; % Why is this used?
            inmP41 = inmk0*EPSZ;
            
            inmP = [NullMat, inmP14; inmP41, NullMat];
        else
            inmP23 = - inmk0*eye(2*loc.Nt + 1);
            inmP32 = - inmk0*EPSZ + (1./inmk0)*inmKX*inmKX;
            
            inmP = ([NullMat, inmP23; inmP32, NullMat]);
        end
        
        % decompose P into form given in 2.138, order by the eigen values and
        
        [V, inmG] = eig(inmP);
   
    
    inmGD=diag(inmG);
    [V, inmG] = sortmat(V,inmGD);
    
    
    % store V and G in Vn and Nn for current layer
    Vn{n+1} = V;
    inmGn{n+1} = inmG;
    
    % Split G into upper and lower diagonal matrices
    inmG_Upper{n + 1} = inmG(1:2*loc.Nt + 1);
    inmG_Lower = inmG(2*loc.Nt + 2:4*loc.Nt + 2);
    
    inmexpG_U = diag(exp( 1i.*nmdz(n).*inmG_Upper{n + 1}));
    inmexpG_L = diag(exp(-1i.*nmdz(n).*inmG_Lower));
    
    % Define the X matrix as in Eqn. 2.144
    X = V\Z{n+1}; % X = inv(V).Z{n+1}
    X_Upper{n + 1} = X(1:2*loc.Nt + 1, :);
    
    % Define the U matrix as in Eqn. 2.146
    U = inmexpG_L*(X(2*loc.Nt + 2 : (4*loc.Nt + 2), :)*(X_Upper{n + 1}\inmexpG_U));
    
    % Calculate the next Z using Eqn. 2.145
    Z{n} = V*[eye(2*loc.Nt + 1); U];
end

% We can now calculate Tn{i}

% A vector describing the excited modes of the incident field.
A=zeros(2*loc.Nt+1,1);  % Y and Z are for the correct pol only, so problem dims are halved
A(loc.Nt+1) = 1;

% Eqn. 2.149 - [T0; R] = [Z(0), - Y_ref] * Y_inc * A

T0R = [Z{1}, - Y_ref]\(Y_inc*A); % Why is this the inverse...??

T0 = T0R(1 : 2*loc.Nt + 1, :);
R = T0R(2*loc.Nt + 2 : (4*loc.Nt + 2), :);

% Build Transmission Vectors
Tn=cell(Ns+1,1);
Tn{1} = T0;


for i = 1: Ns
    inmG_Upper = inmGn{i+1}(1: 2*loc.Nt + 1);
    Tn{i+1} = X_Upper{i+1}\(diag(exp(1i*nmdz(i)*inmG_Upper))*Tn{i}); % Eqn. 2.143
end


% Calculate the electric fields
% Create the f vector which contains [e_x; e_y; h_x; h_y]
f = zeros(Ns,4*loc.Nt+2);
e_x = zeros(Ns,2*loc.Nt+1);
e_y = zeros(Ns,2*loc.Nt+1);
e_z = zeros(Ns,2*loc.Nt+1);
h_x = zeros(Ns,2*loc.Nt+1);
h_y = zeros(Ns,2*loc.Nt+1);
h_z = zeros(Ns,2*loc.Nt+1);

E = zeros(3, Ns, 2*loc.Nt+1);
H = zeros(3, Ns, 2*loc.Nt+1);

for i = 1: Ns
    
    f(i,:) = 0.5*(Z{i}*Tn{i}+Z{i+1}*Tn{i+1}); % Eqn. 2.143
    
    
    if(loc.pol == 1) % p-polarisation
        e_x(i,:) = f(i,1:2*loc.Nt+1);
        h_y(i,:) = f(i,2*loc.Nt+2:4*loc.Nt+2);
        eps = toeplitz(epsf(2*loc.Nt+1:-1:1,i),epsf(2*loc.Nt+1:4*loc.Nt+1,i));
        %inv_eps = toeplitz(inv_epsmat(2*loc.Nt+1:-1:1,i),inv_epsmat(2*loc.Nt+1:4*loc.Nt+1,i));
        %e_z(i,:) = - 1./(im_from_inm(inmk0)*const.c) *  const.eps0 * inv(eps) * im_from_inm(inmKX)*(h_y(i,:)).';
        e_z(i,:) = -1/inmk0 * inv(eps) * (inmKX)*(h_y(i,:)).';
        
        
        h_x(i,:) = 0*e_x(i);
        e_y(i,:) = 0*h_y(i);
        h_z(i,:) = 0*e_z(i);
    elseif(loc.pol == 0) % s-polarisation
        e_y(i,:) = f(i,1:2*loc.Nt+1);
        h_x(i,:) = f(i,2*loc.Nt+2:4*loc.Nt+2);
        % h_z(i,:) = (im_from_inm(inmk0)/(const.c*const.mu0)) * inmKX*(e_y(i,:)).';
        h_z(i,:) = im_from_inm(inmk0) * inmKX*(e_y(i,:)).';
        
        
        e_x(i,:) = 0*h_x(i);
        h_y(i,:) = 0*e_y(i);
        e_z(i,:) = 0*h_z(i);
    end
end

E(1,:,:) = e_x;
E(2,:,:) = e_y;
E(3,:,:) = e_z;

H(1,:,:) = h_x;
H(2,:,:) = h_y;
H(3,:,:) = h_z;
