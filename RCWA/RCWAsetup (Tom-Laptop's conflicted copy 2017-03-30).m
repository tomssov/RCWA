%function [Gx,mAicm2JOpt] = RCWAsetup(varargin)
%Run the RCWA algorithm for a range of wavelengths to build the absorption spectrum

addpath(genpath('./'))
loc = DefaultLoc;
%loc = ApplyVarargin(loc,varargin);
loc = VerifyLoc(loc);

% Create wavelength and wavenumber vectors
nmlambda = linspace(loc.nmlambda0, loc.nmlambda1, loc.nlambda);
inmk0 = 2.*pi./nmlambda;

% Create Material Vectors
nmx = linspace(-0.25*loc.nmLx, 0.75*loc.nmLx, loc.Nx);
[mat_cat, nmz, nmdz] = BuildMaterial(loc);
Nz = length(nmz);
nmgtop = loc.nmda + loc.nmLm;

nmLz = loc.nmLp+loc.nmLi+loc.nmLn;

% Calculate Eg at junction locations
Eg = zeros(length(nmz),1);
Eg(mat_cat == 2) = fliplr(EgProfile(nmz(mat_cat == 2)-nmgtop, loc));

% Create blank matrices for output
r0 = zeros(loc.nlambda,1);
r1 = zeros(loc.nlambda,1);
t0 = zeros(loc.nlambda,1);
t1 = zeros(loc.nlambda,1);

Efl = zeros(loc.nlambda, 3, Nz, 2*loc.Nt+1);

% Matrices for permittivity and fourier of permittivity
epsl = zeros(loc.nlambda,loc.Nx, Nz);
epsf = zeros(loc.nlambda, 4*loc.Nt+1, Nz);

if loc.parforArg > 0
    disp('parallel computation starting ...')
    parfor (ll = 1 : loc.nlambda, loc.parforArg)
        %for ll = 1 : loc.nlambda
        
        % Build eps and epsf matrices
        epsl(ll,:,:) = BuildEps(mat_cat,nmx,nmz, nmlambda(ll), Eg, loc);
        epsf(ll,:,:) = BuildEpsF(nmz, epsl(ll,:,:), loc);
        
        % Perform the RCWA calculation at this wavelength
        [~, Tn, R, Efl(ll,:,:,:),~] ...
            =RCWA(...
            inmk0(ll), epsf(ll,:,:), nmdz,...
            loc ...
            );
        
        % If plotting output is on, calculate the reflection and transmittion
        % vectors
        if loc.plotting == 1
            [r0(ll),r1(ll),t0(ll),t1(ll)]...
                =specular(...
                R,...
                Tn{end},...
                loc, ...
                'inmk0', inmk0(ll), ...
                'nmdz', nmdz);
        end
    end
else
    disp('serial computation starting ...')
    for ll = 1 : loc.nlambda
        
        % Build eps and epsf matrices
        [epsl(ll,:,:),loc] = BuildEps(mat_cat, nmx, nmz, nmlambda(ll), Eg, loc);
        epsf(ll,:,:) = BuildEpsF(nmz, squeeze(epsl(ll,:,:)), loc);
        
        % Perform the RCWA calculation at this wavelength
         [~, Tn, R, Efl(ll,:,:,:)] ...
             =RCWA(...
             inmk0(ll), squeeze(epsf(ll,:,:)), nmdz,...
             loc ...
             );
        
       
        % If plotting output is on, calculate the reflection and transmittion
        % vectors
        if loc.plotting == 1
            [r0(ll),r1(ll),t0(ll),t1(ll)]...
                =specular(...
                R,...
                Tn{end},...
                loc, ...
                'inmk0', inmk0(ll), ...
                'nmdz', nmdz);
        end
    end
end

nmdlambda = (loc.nmlambda1-loc.nmlambda0)/loc.nlambda;
nmdlambda=max(nmdlambda,1);

nmdz = (loc.nmLp+loc.nmLi+loc.nmLn)/loc.Nz;
nmdx = loc.nmLx/(loc.Nx-1);
nmdg = (loc.nmda)/loc.Ng;
nmdm = loc.nmLm/loc.Nm;

W2inmim2Sl = W2inmim2AM15G('nmlambda', nmlambda);
Enormconst = sqrt(2*const.eta0); % Scale E by this and H by 1/this such that incident power density is 1 W/m^2


% Electric Field Variables
El = zeros(loc.nlambda, 3, loc.Nx, Nz);
normE2l = zeros(loc.nlambda,loc.Nx, Nz);

% Energy Absorption Locations
Ql = zeros(loc.nlambda,loc.Nx, Nz);
Q = zeros(loc.Nx, Nz);
Gl = zeros(loc.nlambda, loc.Nx, length(nmz(mat_cat==2)));
G = zeros(loc.Nx, length(nmz(mat_cat==2)));

% Where the energy is absorbed
die_absl = zeros(loc.nlambda,1);
metal_absl = zeros(loc.nlambda,1);

if loc.parforArg >0
    parfor (ll = 1:loc.nlambda, loc.parforArg)
        %for ll = 1:loc.nlambda
        
        Ef = squeeze(Efl(ll,:,:,:));
        
        W2inmim2S = W2inmim2Sl(ll);
        
        El_temp = squeeze(El(ll,:,:,:));
        eps_temp = epsl(ll,:,:);
        
        if loc.pol == 0
            for zz = 1:Nz
                El_temp(2, :, zz) = Enormconst * InverseFourier(nmx', Ef(2,zz,:), loc.nmLx, 1, loc.Nt/2);
            end
        else
            for zz = 1:Nz
                El_temp(1, :, zz) = Enormconst * InverseFourier(nmx', Ef(1,zz,:), loc.nmLx, 1, loc.Nt);
                El_temp(3, :, zz) = Enormconst * InverseFourier(nmx', Ef(3,zz,:), loc.nmLx, 1, loc.Nt);
            end
        end
        
        if loc.epscalc == 1
            for zz = 1:Nz
                eps_temp(:,zz) = InverseFourier(nmx', epsf(:,zz,ll), loc.nmLx, 1, loc.Nt);
            end
        end
        
        El(ll,:,:,:) = El_temp;
        normE2l(ll,:,:) = squeeze(abs(El_temp(1,:,:)).^2 + abs(El_temp(2,:,:)).^2 + abs(El_temp(3,:,:)).^2);
        
        % Add the generation at this wavelength if photon energy > Eg
        Q_temp = (pi/(const.eta0*m_from_nm(nmlambda(ll)))) * abs(imag(epsl(ll,:,:))).* normE2l(ll,:,:);
        Ql(ll,:,:) = Q_temp;
        Q = Q + nmdlambda * squeeze(Ql(ll,:,:));
        
        above_bandgap= sign(eV_from_nm(nmlambda(ll)) - Eg(mat_cat == 2));
        gamma = const.h*const.c/m_from_nm(nmlambda(ll));
        
        Gl(ll,:,:) = (W2inmim2S * Q_temp(:, mat_cat==2)/gamma);
        Gl(ll,:,:) = bsxfun(@times, Gl(ll,:,:), (above_bandgap==1)');
        G = G + nmdlambda * squeeze(Gl(ll,:,:));
        
        % Calculate where the light is absorbed
        die_absl(ll) = sum(sum(Q_temp(:,mat_cat==2)))* m_from_nm(nmdz) /loc.Nx;
        die_absl(ll) = die_absl(ll) + sum(sum(Q_temp(ceil(length(nmx)/2):length(nmx), mat_cat==1)))* m_from_nm(nmdg)/loc.Nx;
        
        metal_absl(ll) = sum(sum(Q_temp(:,mat_cat==0)))* m_from_nm(nmdm) /loc.Nx;
        metal_absl(ll) = metal_absl(ll) + sum(sum(Q_temp(1:floor(length(nmx)/2), mat_cat==1)))* m_from_nm(nmdg) /loc.Nx;
        
    end
else
    for ll = 1:loc.nlambda
        
        Ef = squeeze(Efl(ll,:,:,:));
        
        W2inmim2S = W2inmim2Sl(ll);
        
        El_temp = squeeze(El(ll,:,:,:));
        eps_temp = squeeze(epsl(ll,:,:));
        
        if loc.pol == 0
            for zz = 1:Nz
                El_temp(2, :, zz) = Enormconst * InverseFourier(nmx', Ef(2,zz,:), loc.nmLx, 1, loc.Nt/2);
            end
        else
            for zz = 1:Nz
                El_temp(1, :, zz) = Enormconst * InverseFourier(nmx', Ef(1,zz,:), loc.nmLx, 1, loc.Nt/2);
                El_temp(3, :, zz) = Enormconst * 0*InverseFourier(nmx', Ef(3,zz,:), loc.nmLx, 1, loc.Nt/2);
            end
        end
        
        if loc.epscalc == 1
            for zz = 1:Nz
                eps_temp(:,zz) = InverseFourier(nmx', squeeze(epsf(ll,:,zz)).', loc.nmLx, 1, loc.NFFT);
            end
        end
        
        El(ll,:,:,:) = El_temp;
        normE2l(ll,:,:) = squeeze(abs(El_temp(1,:,:)).^2 + abs(El_temp(2,:,:)).^2 + abs(El_temp(3,:,:)).^2);
        
        % Add the generation at this wavelength if photon energy > Eg
        Q_temp = squeeze((pi/(const.eta0*m_from_nm(nmlambda(ll)))) * abs(imag(epsl(ll,:,:))).* normE2l(ll,:,:));
        Ql(ll,:,:) = Q_temp;
        Q = Q + nmdlambda * squeeze(Ql(ll,:,:));
       
        above_bandgap= sign(eV_from_nm(nmlambda(ll)) - Eg(mat_cat == 2));
        gamma = const.h*const.c/m_from_nm(nmlambda(ll));
        
        Gl_temp = (W2inmim2S * Q_temp(:, mat_cat==2)/gamma);
        Gl(ll,:,:) = bsxfun(@times, Gl_temp, (above_bandgap==1)');
        G = G + nmdlambda * squeeze(Gl(ll,:,:));
        
        % Calculate where the light is absorbed
        die_absl(ll) = sum(sum(Q_temp(:,mat_cat==2)))* m_from_nm(nmdz) /loc.Nx;
        die_absl(ll) = die_absl(ll) + sum(sum(Q_temp(ceil(length(nmx)/2):length(nmx), mat_cat==1)))* m_from_nm(nmdg)/loc.Nx;
        
        metal_absl(ll) = sum(sum(Q_temp(:,mat_cat==0)))* m_from_nm(nmdm) /loc.Nx;
        metal_absl(ll) = metal_absl(ll) + sum(sum(Q_temp(1:floor(length(nmx)/2), mat_cat==1)))* m_from_nm(nmdg) /loc.Nx;
    end
end

if loc.plotting ==1
    % If we have specified a spectrum, plot the absorption spectrums
    if(loc.nlambda ~=1)
        figure(1)
        tot_abs = real(1-t0 - r0 - t1 - r1);
        plot(nmlambda,tot_abs,nmlambda,die_absl,nmlambda,metal_absl,nmlambda,metal_absl+die_absl);
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
    %contourf(nmx,nmz - nmgtop, squeeze(abs(El(1,3,:,:))).',100,'linestyle','none');
    %contourf(nmx, nmz(mat_cat==2) - nmgtop, log10(abs(G))',200,'linestyle','none');
    contourf(nmx, nmz - nmgtop, real(squeeze(epsl(1,:,:)))',100,'linestyle','none');
    %contourf(nmx, nmz - nmgtop, sqrt(squeeze(normE2l(1,:,:)))',100,'linestyle','none');
    %axis([-inf,inf,-inf,inf,0,inf])
    %caxis([0,inf])
    
    axis equal;
    xlim([nmx(1),nmx(end)]);
   % caxis([0,100]);
    colorbar;
    DrawOverlay(loc);
end

% Calculate optical short circuit current density
mAicm2JOpt = const.q * sum(sum(G)) * m_from_nm(nmdz) * nmdx/(10*loc.nmLx);

if(loc.pol == 0)
    fprintf('JOpt_s: %f \n',mAicm2JOpt)
else
    fprintf('JOpt_p: %f \n',mAicm2JOpt)
end

Gx = [nmLz - (nmz(mat_cat==2) - nmgtop) ; sum(G,1) * nmdx/(loc.nmLx)]';