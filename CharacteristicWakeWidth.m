function dw = CharacteristicWakeWidth(x_dist,ti0,ct,A)
%CHARACTERISTICWAKEWIDTH calculates the quantity sigma/D-epsilon, with sigma being the characteristic wake width at a
%certain distance downstream of a rotor with diameter D
%
% Inputs:
%   x_dist (vector): 
%       Downstream distance from rotor in units of the rotor diameter 
%
%   ti0 (vector or scalar):
%       Ambient turbulence intensity. Value must be between 0 an 1 (1=100%)
%
%   ct (two dimensional matrix, first dimension must be length(x_dist)):
%       Thrust coefficients for wake-generating turbines at the relevant wind speeds (second dimension)
%
%   A (scalar):
%       Wake expansion calibration parameter
%
% Outputs:
%   dw (matrix of same size as ct): 
%       Wake width (sigma/D-epsilon) at the specified downstream positions at the relevant wind speeds
%
% Examples:
%   ::
%
%       x_dist = [6; 12]; % two upstream turbines
%       ti0 = [0.07, 0.1 0.12];
%       ct = [0.8, 0.7, 0.4; 0.8 0.7, 0.4]; % Two identical thrust curves
%       dw = CharacteristicWakeWidth(x_dist, ti0, ct, 0.04);
%
% References:
%   - https://github.com/OrstedRD/TurbOPark/blob/main/TurbOPark%20description.pdf
%
% See also:
%   - TurbOPark
%
%   Author: Søren Trads Steen, Jesper Grønnegaard, Nicolai Gayle Nygaard, Sidse Damgaard Hansen
%   Checked by: Cecilia Mortensen Kobæk
%
%   Copyright (c) 2021 by Ørsted


arguments
    x_dist (:,1) double
    ti0 (1,:) double {mustBeGreaterThanOrEqual(ti0,0), mustBeLessThanOrEqual(ti0,1)}
    ct (:,:) double
    A (1,1) double {mustBeNonnegative}
end

%% Check lengths
assert(size(ct,2) == length(ti0) || length(ti0) == 1, 'ti0 must have the same length as second dimension of ct, or ti0 must be a scalar')
assert(size(ct,1) == length(x_dist), 'First dimension of ct must have the same length as x_dist')

%% Define parameters
c1 = 1.5; 
c2 = 0.8; % S. T. Frandsen, “Risø-R-1188(EN) Turbulence and turbulence generated structural loading in wind turbine clusters” Risø, Roskilde, Denmark, 2007.

ti0 = repmat(ti0, size(ct,1), 1); % To enable multiple ti0 values. Would not be necessary if ti0 was always a scalar
alpha = ti0 * c1;
beta = c2 * ti0./sqrt(ct); 

% Formula for characteristic wake width: sigma/rotor_diameter = epsilon + dw
dw = A*ti0./beta .* ( sqrt((alpha+beta.*x_dist).^2+1) - sqrt(1+alpha.^2)...
    - log( ((sqrt((alpha+beta.*x_dist).^2+1)+1).*alpha) ./ ((sqrt(1+alpha.^2)+1).*(alpha+beta.*x_dist)) ) );
end