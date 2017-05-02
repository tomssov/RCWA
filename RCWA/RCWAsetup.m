function [Gx,mAicm2JOpt] = RCWAsetup(varargin)
%Run the RCWA algorithm for a range of wavelengths to build the absorption spectrum

addpath(genpath('./'))
%loc = MishaLoc;
loc = DefaultLoc;
loc = ApplyVarargin(loc,varargin);
loc = VerifyLoc(loc);

% Create wavelength and wavenumber vectors
nmlambda = linspace(loc.nmlambda0, loc.nmlambda1, loc.nlambda);
inmk0 = 2.*pi./nmlambda;
radtheta = rad_from_deg(linspace(loc.degtheta0, loc.degtheta1, loc.ntheta));

% Create Material Vectors
nmx = linspace(-0.25*loc.nmLx, 0.75*loc.nmLx, loc.Nx);
nmx = linspace(-0.5*loc.nmLx, 0.5*loc.nmLx, loc.Nx);
[mat_cat, nmz, nmdz] = BuildMaterial(loc);
Nz = length(nmz);
nmgtop = loc.nmda + loc.nmLm;

nmLz = loc.nmLp+loc.nmLi+loc.nmLn;

% Calculate Eg at junction locations
Eg = zeros(length(nmz),1);
Eg(mat_cat == 2) = fliplr(EgProfile(nmz(mat_cat == 2)-nmgtop, loc));

% Create blank matrices for output
r0 = zeros(loc.nlambda, loc.ntheta, 1);
r1 = zeros(loc.nlambda, loc.ntheta, 1);
t0 = zeros(loc.nlambda, loc.ntheta, 1);
t1 = zeros(loc.nlambda, loc.ntheta, 1);

Efl = zeros(loc.nlambda, loc.ntheta, 3, Nz, 2*loc.Nt+1);

% Matrices for permittivity and fourier of permittivity
epsl = zeros(loc.nlambda, loc.Nx, Nz);
epsf = zeros(loc.nlambda, Nz, 4*loc.Nt+1);

if 1 ==1%loc.recalc ~=0
if loc.parforArg 
    disp('parallel computation starting ...')
    for ll = 1 : loc.nlambda
        
        % Build eps and epsf matrices
        [epsl(ll,:,:),loc] = BuildEps(mat_cat, nmx, nmz, nmlambda(ll), Eg, loc);
        epsf(ll,:,:) = BuildEpsF(nmz, squeeze(epsl(ll,:,:)), loc);
        
        for tt = 1 : loc.ntheta
            
            % Perform the RCWA calculation at this wavelength
            [~, Tn, R, Efl(ll,tt,:,:,:)] ...
                =RCWA(...
                inmk0(ll), radtheta(tt), squeeze(epsf(ll,:,:)), nmdz,...
                loc ...
                );
            
            
            % If plotting output is on, calculate the reflection and transmittion
            % vectors
            if loc.plotting == 1
                [r0(ll,tt),r1(ll,tt),t0(ll,tt),t1(ll,tt)]...
                    =specular(...
                    R,...
                    Tn{end},...
                    inmk0(ll), radtheta(tt), ...
                    loc ...
                    );
            end
        end
    end
else
    disp('serial computation starting ...')
    for ll = 1 : loc.nlambda
        
        % Build eps and epsf matrices
        [epsl(ll,:,:),loc] = BuildEps(mat_cat, nmx, nmz, nmlambda(ll), Eg, loc);
        epsf(ll,:,:) = BuildEpsF(nmz, squeeze(epsl(ll,:,:)), loc);
        
        for tt = 1 : loc.ntheta
            
            % Perform the RCWA calculation at this wavelength
            [~, Tn, R, Efl(ll,tt,:,:,:)] ...
                =RCWA(...
                inmk0(ll), radtheta(tt), squeeze(epsf(ll,:,:)), nmdz,...
                loc ...
                );
            
            
            % If plotting output is on, calculate the reflection and transmittion
            % vectors
            if loc.plotting == 1
                [r0(ll,tt),r1(ll,tt),t0(ll,tt),t1(ll,tt)]...
                    =specular(...
                    R,...
                    Tn{end},...
                    inmk0(ll), radtheta(tt), ...
                    loc ...
                    );
            end
        end
    end
end
end

nmdlambda = (loc.nmlambda1-loc.nmlambda0)/loc.nlambda;
nmdlambda=max(nmdlambda,1);

nmdJ = (loc.nmLp+loc.nmLi+loc.nmLn)/loc.Nz;
nmdx = loc.nmLx/(loc.Nx-1);
nmdg = (loc.nmda)/loc.Ng;
nmdm = loc.nmLm/loc.Nm;

W2inmim2Sl = W2inmim2AM15G('nmlambda', nmlambda);
Enormconst = sqrt(2*const.eta0); % Scale E by this and H by 1/this such that incident power density is 1 W/m^2


% Electric Field Variables
El = zeros(loc.nlambda, loc.ntheta, 3, loc.Nx, Nz);
normE2l = zeros(loc.nlambda, loc.ntheta,loc.Nx, Nz);

% Energy Absorption Locations
Ql = zeros(loc.nlambda,loc.Nx, Nz);
Q_temp = zeros(1,loc.Nx, Nz);
Q = zeros(loc.Nx, Nz);
Gl = zeros(loc.nlambda, loc.Nx, length(nmz(mat_cat==2)));
G = zeros(loc.Nx, length(nmz(mat_cat==2)));

% Where the energy is absorbed
die_absl = zeros(loc.nlambda,loc.ntheta);
metal_absl = zeros(loc.nlambda,loc.ntheta);


if loc.parforArg >0
   parfor (ll = 1:loc.nlambda, loc.parforArg)
        Q_temp = zeros(1,loc.Nx, Nz);
        for tt = 1:loc.ntheta
            
            Ef = squeeze(Efl(ll,tt,:,:,:));

            Et = zeros(loc.ntheta, 3, loc.Nx, Nz);
            normE2t= zeros(loc.ntheta, loc.Nx, Nz);

            die_abst = zeros(loc.ntheta,1);
            metal_abst = zeros(loc.ntheta,1);
    
     
            E_temp = zeros(3, loc.Nx, Nz);
            eps_temp = zeros(loc.Nx, Nz);
            
            if loc.pol == 0
                for zz = 1:Nz
                    E_temp(2, :, zz) = Enormconst * InverseFourier(nmx', Ef(2,zz,:), loc.nmLx, 0, loc.Nt);
                end
            else
                for zz = 1:Nz
                    E_temp(1, :, zz) = Enormconst * InverseFourier(nmx', Ef(1,zz,:), loc.nmLx, 0, loc.Nt);
                    E_temp(3, :, zz) = Enormconst * InverseFourier(nmx', Ef(3,zz,:), loc.nmLx, 0, loc.Nt);
                end
            end
            
            if loc.epscalc == 1
                for zz = 1:Nz
                    eps_temp(:,zz) = InverseFourier(nmx', epsf(ll,zz,:), loc.nmLx, 0, loc.Nt); %error here
                end
                epsl(ll,:,:) = eps_temp;
            end

            Et(tt,:,:,:) = E_temp;
            normE2t(tt,:,:) = squeeze(abs(E_temp(1,:,:)).^2 + abs(E_temp(2,:,:)).^2 + abs(E_temp(3,:,:)).^2);
            
            % Add the generation at this wavelength if photon energy > Eg
            Q_temp_theta = (pi/(const.eta0*m_from_nm(nmlambda(ll)))) * abs(imag(epsl(ll,:,:))).* reshape(normE2t(tt,:,:),1,loc.Nx,Nz);
            
            % Calculate where the light is absorbed
            die_abst(tt) = sum(sum(Q_temp_theta(1,:,mat_cat==2)))* m_from_nm(nmdJ) /loc.Nx;
            die_abst(tt) = die_abst(tt) + sum(sum(Q_temp_theta(1,ceil(length(nmx)/2):length(nmx), mat_cat==1)))* m_from_nm(nmdg)/loc.Nx;
            
            metal_abst(tt) = sum(sum(Q_temp_theta(1,:,mat_cat==0)))* m_from_nm(nmdm) /loc.Nx;
            metal_abst(tt) = metal_abst(tt) + sum(sum(Q_temp_theta(1,1:floor(length(nmx)/2), mat_cat==1)))* m_from_nm(nmdg) /loc.Nx;
            
            Q_temp = Q_temp + (Q_temp_theta);
        end
        
        El(ll,:,:,:,:) = Et;
        normE2l(ll,:,:,:) = normE2t;
        die_absl(ll,:) = die_abst;
        metal_absl(ll,:) = metal_abst;
        
        Q_temp = (Q_temp/loc.ntheta);
        
        Ql(ll,:,:) = Ql(ll,:,:) + Q_temp;
        Q = Q + squeeze(Ql(ll,:,:));
        
        above_bandgap= sign(eV_from_nm(nmlambda(ll)) - Eg(mat_cat == 2));
        
        gamma = const.h*const.c/m_from_nm(nmlambda(ll));
        Gl(ll,:,:) = (W2inmim2Sl(ll) * Q_temp(1, :, mat_cat==2)/gamma);
        Gl(ll,:,:) = bsxfun(@times, squeeze(Gl(ll,:,:)), (above_bandgap==1)');
        G = G + nmdlambda*squeeze(Gl(ll,:,:));
    end 
else
    for ll = 1:loc.nlambda
           Q_temp = zeros(1,loc.Nx, Nz);
        for tt = 1:loc.ntheta
            
            Ef = squeeze(Efl(ll,tt,:,:,:));

            Et = zeros(loc.ntheta, 3, loc.Nx, Nz);
            normE2t= zeros(loc.ntheta, loc.Nx, Nz);

            die_abst = zeros(loc.ntheta,1);
            metal_abst = zeros(loc.ntheta,1);
    
     
            E_temp = zeros(3, loc.Nx, Nz);
            eps_temp = zeros(loc.Nx, Nz);
            
            if loc.pol == 0
                for zz = 1:Nz
                    E_temp(2, :, zz) = Enormconst * InverseFourier(nmx', Ef(2,zz,:), loc.nmLx, 0, loc.Nt);
                end
            else
                for zz = 1:Nz
                    E_temp(1, :, zz) = Enormconst * InverseFourier(nmx', Ef(1,zz,:), loc.nmLx, 0, loc.Nt);
                    E_temp(3, :, zz) = Enormconst * InverseFourier(nmx', Ef(3,zz,:), loc.nmLx, 0, loc.Nt);
                end
            end
            
            if loc.epscalc == 1
                for zz = 1:Nz
                    eps_temp(:,zz) = InverseFourier(nmx', epsf(ll,zz,:), loc.nmLx, 0, loc.Nt); %error here
                end
                epsl(ll,:,:) = eps_temp;
            end

            Et(tt,:,:,:) = E_temp;
            normE2t(tt,:,:) = squeeze(abs(E_temp(1,:,:)).^2 + abs(E_temp(2,:,:)).^2 + abs(E_temp(3,:,:)).^2);
            
            % Add the generation at this wavelength if photon energy > Eg
            Q_temp_theta = (pi/(const.eta0*m_from_nm(nmlambda(ll)))) * abs(imag(epsl(ll,:,:))).* reshape(normE2t(tt,:,:),1,loc.Nx,Nz);
            
            % Calculate where the light is absorbed
            die_abst(tt) = sum(sum(Q_temp_theta(1,:,mat_cat==2)))* m_from_nm(nmdJ) /loc.Nx;
            die_abst(tt) = die_abst(tt) + sum(sum(Q_temp_theta(1,ceil(length(nmx)/2):length(nmx), mat_cat==1)))* m_from_nm(nmdg)/loc.Nx;
            
            metal_abst(tt) = sum(sum(Q_temp_theta(1,:,mat_cat==0)))* m_from_nm(nmdm) /loc.Nx;
            metal_abst(tt) = metal_abst(tt) + sum(sum(Q_temp_theta(1,1:floor(length(nmx)/2), mat_cat==1)))* m_from_nm(nmdg) /loc.Nx;
            
            Q_temp = Q_temp + (Q_temp_theta);
        end
        
        El(ll,:,:,:,:) = Et;
        normE2l(ll,:,:,:) = normE2t;
        die_absl(ll,:) = die_abst;
        metal_absl(ll,:) = metal_abst;
        
        Q_temp = (Q_temp/loc.ntheta);
        
        Ql(ll,:,:) = Ql(ll,:,:) + Q_temp;
        Q = Q + squeeze(Ql(ll,:,:));
        
        above_bandgap= sign(eV_from_nm(nmlambda(ll)) - Eg(mat_cat == 2));
        
        gamma = const.h*const.c/m_from_nm(nmlambda(ll));
        Gl(ll,:,:) = (W2inmim2Sl(ll) * Q_temp(1, :, mat_cat==2)/gamma);
        Gl(ll,:,:) = bsxfun(@times, squeeze(Gl(ll,:,:)), (above_bandgap==1)');
        G = G + nmdlambda*squeeze(Gl(ll,:,:));
    end
end

if loc.plotting ==1
    % If we have specified a spectrum, plot the absorption spectrums
    if(loc.nlambda ~=1)
        figure(1)
        tot_abs = abs(1-t0 - r0 - t1 - r1);

        hold off;
        for tt=1:loc.ntheta
            plot(nmlambda,tot_abs(:,tt),nmlambda,die_absl(:,tt),nmlambda,metal_absl(:,tt),nmlambda,metal_absl(:,tt)+die_absl(:,tt));
        hold on;
        end
        hold off;
        
        xlabel('\lambda (nm)','FontSize',16);
        ylabel('Absorption','FontSize',16);
        if(loc.pol == 0)
            title('Specular absorption for s-polarization');
        elseif(loc.pol == 1)
            title('Specular absorption for p-polarization');
        end
        legend ('A=1-R-T','Dielectric','Metal','Die'' + Metal','location','southwest');
        ylim([0,1]);
    end
    
    
    figure(2);
    colormap(jet);
 
    %contourf(nmx,nmz - nmgtop, squeeze(abs(El(1,1,1,:,:))).','LevelStep',0.2,'linestyle','none');
    %contourf(nmx,nmz - nmgtop, abs(squeeze(eps_temp))'.*squeeze(abs(El(1,1,3,:,:))).','LevelStep',1,'linestyle','none');
    %contourf(nmx,nmz - nmgtop, real(squeeze(eps_temp))','LevelStep',0.1,'linestyle','none');
    contourf(nmx, nmz(mat_cat==2) - nmgtop, log10(abs(G))','LevelStep', 0.01, 'linestyle', 'none');
    %contourf(nmx, nmz - nmgtop, real(squeeze(eps_temp(:,:)))',100,'linestyle','none');
    %contourf(nmx, nmz - nmgtop, sqrt(squeeze(normE2l(1,1,:,:)))','LevelStep',0.1,'linestyle','none');
    %contourf(nmx, nmz - nmgtop, sqrt(squeeze(Q(1,1,:,:)))','LevelStep',0.2,'linestyle','none');
    %axis([-inf,inf,-inf,inf,0,inf])
    %caxis([0,inf])
    caxis([0,30]);
    axis equal;
    xlim([nmx(1),nmx(end)]);
   
    colorbar;
    %DrawOverlay(loc);
    
    if(loc.ntheta > 1)
        figure(3);
        contourf(nmlambda, linspace(loc.degtheta0,loc.degtheta1,loc.ntheta), die_absl', 200, 'linestyle','none');
    end
    
end

% Calculate optical short circuit current density
mAicm2JOpt = const.q * sum(sum(G)) * m_from_nm(nmdJ) * nmdx/(10*loc.nmLx);

if(loc.pol == 0)
    fprintf('JOpt_s: %f \n',mAicm2JOpt)
else
    fprintf('JOpt_p: %f \n',mAicm2JOpt)
end

%figure(1);
%plot(nmx, sqrt(squeeze(normE2l(1,1,:,loc.Nw+loc.Nair+loc.Nz+3)))');

Gx = [nmLz - (nmz(mat_cat==2) - nmgtop) ; sum(G,1) * nmdx/(loc.nmLx)]';