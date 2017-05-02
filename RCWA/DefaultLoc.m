function [loc] = DefaultLoc()
%DEFAULTLOC Summary of this function goes here
%   Detailed explanation goes here

loc = struct(...
    ...% Computation variablesd
    ...% 'recalc', 0, ...   % If 0, only plots are reconstructed
    'parforArg', inf, ...     % Numer of parallel agents
    'plotting', 0, ...      % Produce plots (0: no, 1: yes)
    'epscalc', 0, ...       % Eps exact (0) of fourier (1) 
    'pol', 1, ...           % Polarization State (0: s-polarization, 1: p-polaization)
    'degtheta0', 0, ...     % Minimum angle of Incident Light
    'degtheta1', 0, ...     % Maximum angle of Incident Light
    'ntheta', 1, ...        % Number of angles
    'Nt', 2, ...            % Number of modes
    'nmlambda0', 300, ...   % Minimum wavelength in nm
    'nmlambda1', 1240/1.6, ...   % Maximum wavelength in nm
    'nlambda', 10, ...     % Number of wavelengths
    'Nx', 50, ...          % Number x-direction sample points   
    ...% Grating Specification
    'nmLx', 400, ...        % Grating period
    'zeta', 0.5, ...        % Graing Duty Cycle
    'type', 0, ...          % Grating Shape (0 for rect., 1 for sine)
    'nmda', 0, ...          % Grating region thickness
    'nmLg', 0, ...          % Grating relief - peak to trough
    'Ng', 1, ...            % Number of grating slices
    'nmLm', 100, ...        % Backing mirror thickness (from grating trough to back)
    'Nm', 1, ...            % Number of grating slices
    'gmat', 1, ...          % Grating material ( 1 for silver, const otherwise)
    'epsm', 1+3i, ...          % Grating metal permittivity
    'dmat', 2, ...          % Grating material (2 for AZO, 1 for glass, const otherwise)
    'epsd', 11.559 + 0.204i, ... % Grating dielectric permittivity
    'NFFT', 900, ...        % Number of terms for FFT of grating
    'ftype', 1, ...         % Grating Fourier method (0 for FFT, 1 for explicit)
    ...% Junction Specification
    'nmLn', 15, ...         % n-layer thickness
    'nmLi', 200, ...        % i-layer thickness
    'nmLp', 15, ...         % p-layer thickness
    'Nz', 230/5, ...          % Number of slices in junction
    'epsJ', 11.559 + 0.204i,...   % Junction permittivity
    'material', 1,...       % Material Specification (0 constant at epsJ,
    ...% 1 aSiH with bandgap as in eVEgz.m, 2 Faryad SPIE testcase)
    ...% Junction Material Specification (used if materlial = 1)
    'Eg0', 1.6, ...         % Baseline Bandgap
    'A', 0, ...             % Perturbation Amplitude
    'kappa', 0, ...         % Periods
    'phi', 0, ...           % Phase
    'alpha', 0, ...         % Shaping Parameter
    ...% Window Specification
    'nmLw', 75, ...        % Thickness of Window
    'epsW', 4+1e-6i, ...    % Permittivity of Window
    'Nw', 1, ...            % Number of window slices
    'wmat', 2, ...          % Window material (2 AZO, 1 for glass, const otherwise)
    ...% Air Specification
    'nmLair', 1000, ...     % Air in front of solar cell
    'Nair', 1, ...          % Number of z-slices in air
    'nsa', 1+1e-6i ...      % Refractive index of the containing medium (approx 1 for air)
    );

end