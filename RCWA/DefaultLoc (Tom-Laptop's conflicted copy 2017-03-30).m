function [loc] = DefaultLoc()
%DEFAULTLOC Summary of this function goes here
%   Detailed explanation goes here

loc = struct(...
    ...% Simulation variables
    'parforArg', 0, ...     % Set the maximum numer of parallel worker. Set to 0 for serial code.
    'plotting', 1, ...      % Produce plots (0: no, 1: yes)
    'epscalc', 1, ...       %
    'pol', 1, ...           % Polarization State (0: s-polarization, 1: p-polaization)
    'degtheta', 0, ...      % Angle of Incident Light
    'Nt', 2, ...            % Number of modes
    'nmlambda0', 700, ...   % Minimum wavelength in nm
    'nmlambda1', 700, ...   % Maximum wavelength in nm
    'nlambda', 1, ...     % Number of wavelengths
    'nsa', 1+1e-6i, ...     % Refractive index of the containing medium (approx 1 for air)
    'Nx', 200, ...           % Number x-direction sample points
    'nmLair', 0*1000, ...      % Air in front of solar cell
    'Nair', 0, ...          %
    ...% Grating Specification
    'nmLx', 800, ...        % Grating period
    'zeta', 0.5, ...        % Graing Duty Cycle
    'type', 2, ...          % Grating Shape
    'nmLg', 150, ...         % Grating relief - peak to trough
    'nmLm', 1*150, ...        % Backing mirror thickness (from grating trough to back)
    'Nm', 10, ...            % Number of grating slices
    'nmda', 100, ...         % Grating region thickness
    'Ng', 200, ...           % Number of grating slices
    'gmat', 1, ...          % Grating material ( 1 for silver, const otherwise)
    'epsm', 0, ...          % Grating metal permittivity
    'dmat', 1, ...          % Grating material ( 1 for glass, const otherwise)
    'epsd', 11.559 + 0.204i, ... % Grating dielectric permittivity
    'NFFT', 900, ...        % Number of terms for FFT of grating
    'ftype', 0, ...         % Grating Fourier method (0 for FFT, 1 for explicit)
    ...% Junction Specification
    'epsJ', 11.559 + 0.204i,...   % Junction permittivity
    'Nz', 20, ...          % Number of slices in junction
    'material', 1,...       % Material Specification (0 constant at epsJ,
    ... % 1 aSiH with bandgap as in eVEgz.m, 2 Faryad SPIE testcase)
    'nmLn', 0*15, ...         % n-layer thickness
    'nmLi', 200, ...        % i-layer thickness
    'nmLp', 0*15, ...         % p-layer thickness
    'Eg0', 1.5, ...         % Baseline Bandgap
    'A', 0.0, ...           % Perturbation Amplitude
    'kappa', 2, ...         % Periods
    'phi', 0.75, ...        % Phase
    'alpha', 5, ...         % Shaping Parameter
    ...% Window Specification
    'nmLw', 0*70, ...         % Thickness of Window
    'epsw', 4+1e-6i, ...    % Permittivity of Window
    'Nw', 0, ...            % Number of window slices
    'wmat', 1 ...           % Window material (1 for glass, const otherwise)
    );

end

