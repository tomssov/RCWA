function [loc] = MishaLoc()
%DEFAULTLOC Summary of this function goes here
%   Detailed explanation goes here

loc = struct(...
    ...% Computation variables
   ...% 'recalc',0, ...         % If 0, only plots are reconstructed
    'parforArg',0,...       % Numer of parallel agents
    'plotting', 1, ...      % Produce plots (0: no, 1: yes)
    'epscalc', 0, ...       % Eps exact (0) of fourier (1) 
    'pol', 1, ...           % Polarization State (0: s-polarization, 1: p-polaization)
    'degtheta0', 0, ...     % Minimum angle of Incident Light
    'degtheta1', 0, ...     % Maximum angle of Incident Light
    'ntheta', 1, ...        % Number of angles
    'Nt', 14, ...           % Number of modes
    'nmlambda0', 562, ...   % Minimum wavelength in nm
    'nmlambda1', 562, ...   % Maximum wavelength in nm
    'nlambda', 1, ...     % Number of wavelengths
    'Nx', 300, ...          % Number x-direction sample points   
    ...% Grating Specification
    'nmLx', 500, ...        % Grating period
    'zeta', 0.5, ...        % Graing Duty Cycle
    'type', 0, ...          % Grating Shape (0 for rect., 1 for sine)
    'nmda', 0, ...         % Grating region thickness
    'nmLg', 40, ...         % Grating relief - peak to trough
    'Ng', 20, ...           % Number of grating slices
    'nmLm', 150, ...        % Backing mirror thickness (from grating trough to back)
    'Nm', 20, ...           % Number of grating slices
    'gmat', 1, ...          % Grating material ( 1 for silver, const otherwise)
    'epsm', 0, ...          % Grating metal permittivity
    'dmat', 0, ...          % Grating material ( 1 for glass, const otherwise)
    'epsd', 11.559 + 0.204i, ... % Grating dielectric permittivity
    'NFFT', 900, ...        % Number of terms for FFT of grating
    'ftype', 1, ...         % Grating Fourier method (0 for FFT, 1 for explicit)
    ...% Junction Specification
    'nmLn', 0, ...          % n-layer thickness
    'nmLi', 275, ...        % i-layer thickness
    'nmLp', 0, ...          % p-layer thickness
    'Nz', 100, ...          % Number of slices in junction
    'epsJ', 11.559 + 0.204i,...   % Junction permittivity
    'material', 0,...       % Material Specification (0 constant at epsJ,
    ... % 1 aSiH with bandgap as in eVEgz.m, 2 Faryad SPIE testcase)
    ...% Junction Material Specification (used if materlial = 1)
    'Eg0', 0, ...           % Baseline Bandgap
    'A', 0, ...             % Perturbation Amplitude
    'kappa', 0, ...         % Periods
    'phi', 0, ...           % Phase
    'alpha', 0, ...         % Shaping Parameter
    ...% Window Specification
    'nmLw', 75, ...         % Thickness of Window
    'epsW', 4+1e-6i, ...    % Permittivity of Window
    'Nw', 20, ...           % Number of window slices
    'wmat', 0, ...          % Window material (1 for glass, const otherwise)
    ...% Air Specification
    'nmLair', 1000, ...     % Air in front of solar cell
    'Nair', 100, ...        % Number of z-slices in air
    'nsa', 1+1e-6i ...      % Refractive index of the containing medium (approx 1 for air)
    );

end

