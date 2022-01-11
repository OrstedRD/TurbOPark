function [power, u] = TurbOPark(u0, direction, u_corr, x_utm, y_utm, hub_height, power_curve,...
    power_curve_index, ti0, options)
%TurbOPark Calculates wind farm wakes given wind conditions and turbine positions and characteristics.
% The function uses the TurbOPark turbulence intensity augmented Jensen model developed by Søren Trads Steen, Jesper
% Grønnegaard Pedersen, and Nicolai Gayle Nygaard.
%
% Inputs:
%   u0 (vector of length n_u0): 
%       Free wind speeds to calculate waked wind speeds and turbine power for
%
%   direction (scalar): 
%       Direction of free wind
%
%   u_corr (vector of length n_WTG):
%       Correction factor relating the free wind speed at the reference height and location to the free wind speed at
%       the hub height and location of each turbine 
%
%   x_utm (vector of length n_WTG):
%       x-coordinates of the turbines in UTM coordinates
%
%   y_utm (vector of length n_WTG):
%       y-coordinates of the turbines in UTM coordinates
%
%   hub_height (vector of length n_WTG):
%       Hub height of turbines
%
%   power_curve (array of structs):
%       The unique "power curves" (turbine types) present in the wind farm
%       Each struct must have the fields
%       - rotor_diameter (scalar): rotor diameter of the given turbine type
%       - interpolant_ct (griddedInterpolant): interpolant describing the
%           relation between wind speed and thrust of the turbine type
%       - interpolant_power (griddedInterpolant): interpolant describing the
%           relation between wind speed and power of the turbine type
%
%   power_curve_index (vector of integers of length n_WTG):
%       Index for each turbine pointing to the relevant entry in the power_curve array defined above
%
%   ti0 (vector of length n_u0 or scalar):
%       Ambient turbulence intensity. Value must be between 0 an 1 (1=100%)
%
%   'A' [opt](scalar):
%       Wake expansion calibration parameter
%
%   'sigma_max_rel' [opt](scalar):
%       We do not include the wake impact from turbine i on turbine j if the radial distance between the two turbines
%       is larger than 0.5*sigma_max_rel*sigma+R. Here sigma is the characteristic wake width and R is the rotor-radius
%       of turbine j.
%       The default value is 4 which is what Ørsted's validation of the TurbOPark model is based on
%
% Outputs:
%   power (matrix of size [n_WTG, n_u0]): 
%       Power output for each turbine accounting for wake effects
%
%   u (matrix of size [n_WTG, n_u0]): 
%       Wind speed at each turbine position accounting for wake effects
%
% Example:
%
%       u0 = [4 8 12 16];
%       direction = 270;
%       u_corr = [1.01 1.02 1.03];
%       x_utm = [500 0 1000];
%       y_utm = [0 100 200];
%       hub_height = [150 150 120];
%       pc_A.interpolant_power = griddedInterpolant(0:26, [linspace(0,3000,26) 0], 'linear','nearest'); 
%       pc_A.interpolant_ct = griddedInterpolant(0:26, [linspace(1,0,26) 0], 'linear','nearest');
%       pc_A.rotor_diameter = 100;
%       pc_B = pc_A; pc_B.rotor_diameter = 130;
%       power_curve = [pc_A pc_B]; 
%       power_curve_index = [1 1 2];
%       ti0 = [0.5 0.4 0.3 0.2];
%       [power, u] = TurbOPark(...
%           u0,direction,u_corr,x_utm,y_utm,hub_height,power_curve,power_curve_index,ti0);
%
% References: 
%   - https://github.com/OrstedRD/TurbOPark/blob/main/TurbOPark%20description.pdf
%
% See also:
%   - CHARACTERISTICWAKEWIDTH
%   - CREATEGAUSSOVERLAPLOOKUPTABLE
%
%   Author: Søren Trads Steen, Jesper Grønnegaard, Nicolai Gayle Nygaard, Sidse Damgaard Hansen
%   Checked by: Cecilia Mortensen Kobæk
%
%   Copyright (c) 2021 by Ørsted

arguments
    u0 (1,:) double {mustBeNonnegative}
    direction (1,1) double {mustBeGreaterThanOrEqual(direction,0), mustBeLessThan(direction,360)}
    u_corr (1,:) double 
    x_utm (1,:) double 
    y_utm (1,:) double
    hub_height (1,:) 
    power_curve (1,:) struct
    power_curve_index (1,:) {mustBeInteger} 
    ti0 (1,:) double {mustBeGreaterThanOrEqual(ti0,0), mustBeLessThanOrEqual(ti0,1)}
    options.A (1,1) double {mustBeNonnegative} = 0.04
    options.sigma_max_rel (1,1) double {mustBeNonnegative} = 4
end

% Test lengths
n_u0 = length(u0);
n_WTG = length(x_utm);

assert(length(y_utm) == n_WTG, 'x and y must have the same length')
assert(length(u_corr) == n_WTG, 'ws_corr must have the same length as the number of turbines')
assert(length(hub_height) == n_WTG, 'hub_height must have the same length as the number of turbines')
assert(length(power_curve_index) == n_WTG, 'power_curve_index must have the same length as the number of turbines')
assert(max(power_curve_index) <= length(power_curve), 'You cannot assign a turbine to a power curve that does not exist')
assert(length(ti0) == n_u0 | length(ti0) == 1, 'ti0 must be same length as u0 or be a scalar')

% Check the fields of the power_curve
for i = 1:length(power_curve)
    assert( all( isfield(power_curve(i), {'rotor_diameter', 'interpolant_power', 'interpolant_ct'}) ), 'power_curve must have the fields: rotor_diameter, interpolant_power, interpolant_ct' )
    assert( isscalar( power_curve(i).rotor_diameter ) && power_curve(i).rotor_diameter > 0, 'rotor_diameter must be positive scalar')
    assert( isa( power_curve(i).interpolant_ct, 'griddedInterpolant' ), 'interpolant_ct must be a griddedInterpolant' )
    assert( isa( power_curve(i).interpolant_power, 'griddedInterpolant' ), 'interpolant_power must be a griddedInterpolant' )
end

%% Load results of numerical integration of gaussian deficit over rotor disk, and set up interpolator
val = load('gauss_lookup_table.mat'); % Pre-calculated values for overlap between rotor disk and gaussian wake 
% For more details look at CREATEGAUSSOVERLAPLOOKUPTABLE (used for
% generating 'gauss_lookup_table.mat') and in the reference note.

% Set up interpolant. Note that dist and radius_down is in units of sigma 
overlap_gauss = griddedInterpolant({val.overlap_lookup_table.dist, ...
    val.overlap_lookup_table.radius_down}, ...
    val.overlap_lookup_table.overlap_gauss,'linear','nearest');

%% Define turbine positions in rotated coordinate system with the x-axis aligned with the wind direction
% Note that the wind direction is measured from the North (the positive y-axis in the original coordinate system)
% increasing in the clockwise direction.

% Coordinate transformation matrix
m = [-sind(direction) -cosd(direction);
    cosd(direction) -sind(direction)];
% Note that this transforms the coordinate system - not the vectors
% (see http://mathworld.wolfram.com/RotationMatrix.html)

% Transform the turbine coordinates
x_y = m*[x_utm ; y_utm];
x = x_y(1,:);
y = x_y(2,:);

%% Sort input from most upwind to most downwind
% Sort turbine positions from most upwind to most downwind
[~, x_index] = sort(x,'ascend');
x = x(:,x_index);
y = y(:,x_index);

% Sort power curve indices, hub heights, and correction factor in same order as x and y
power_curve_index = power_curve_index(x_index);
hub_height = hub_height(x_index);
u_corr = u_corr(x_index);

% Get the rotor diameter of each turbine and sort in the same order as x and y (and make sure it's a row vector)
rotor_diameter = [power_curve.rotor_diameter];
rotor_diameter = rotor_diameter(power_curve_index);
rotor_diameter = rotor_diameter(:)';

%% Initiate variables 

% Predicted wind speed at each turbine position. Initially set equal to the
% freestream reference wind speed
u = u0.*ones(n_WTG,n_u0); 

% Predicted Ct for each turbine
ct = zeros(n_WTG,n_u0);

% Predicted power for each turbine
power = zeros(n_WTG,n_u0);

% Total wake deficit affecting each turbine
delta_total = zeros(n_WTG,n_u0);

%% Loop through all turbines
% Calculate the wake effect on each turbine from all other upstream
% turbines
for i = 1:n_WTG   
    % Distances along x between the i'th turbine and all other turbines
    % x_dist(j) are non-negative for upwind turbines (j < i)
    x_dist = x(i)-x;
    
    % Normalize x_dist(j) with rotor diameter of turbine j
    x_dist = x_dist./rotor_diameter;
    
    % Radial distance between the rotor center of turbine i and the centerlines of wakes from all turbines
    r_dist = sqrt((y(i)-y).^2+(hub_height(i)-hub_height).^2);    
    % Radial distance between the rotor center of turbine i and the centerlines of wakes from all image turbines
    r_dist_image = sqrt((y(i)-y).^2+(hub_height(i)-(-hub_height)).^2);
    % (Image turbines with negative hub height are used to account for ground effects)
    
    % Characteristic width of wakes from all turbines at the position of the i'th turbine (ct for upstream turbines has
    % been calculated in previous iterations, and the width is NaN for downstream WTGs where ct is still set to 0) 
    dw = CharacteristicWakeWidth(x_dist, ti0, ct, options.A); % sigma/rotor_diameter = epsilon + dw
    epsilon = 0.25*sqrt(min(0.5*(1+sqrt(1-ct))./sqrt(1-ct),3));
    sigma = rotor_diameter.*(epsilon + dw)';
    
    % Peak wake deficits
    C = 1-sqrt(1-ct'./(8*(sigma./rotor_diameter).^2));
    % (sigma and thereby C is NaN for downstream WTGs where ct is still set to 0)
    
    % Find upstream WTGs with wakes overlapping the rotor of turbine i and calculate the deficit-contribution from each.
    % (The Gaussian wake is in principle infinitely wide, but in practice we only consider its impact out to a certain
    % distance (defined by sigma_max_rel). Numerical experiments show that sigma_max_rel >= 4 gives adequately 
    % converged AEP results.)
    effective_width = options.sigma_max_rel*sigma;
    is_overlapping = effective_width/2 + rotor_diameter(i)/2 > r_dist;
    
    % We loop over upstream WTGs that for any wind speed (u0) is found to wake WTG i
    wtg_overlapping = x_dist > 0 & any(is_overlapping,1);
    wtg_overlapping_idx = find(wtg_overlapping);
    
    if any(wtg_overlapping)
    
        delta_real = nan(length(u0), n_WTG); % Wake deficit from an actual turbine
        delta_image = nan(length(u0), n_WTG); % Wake deficit from an image turbine
        for j = wtg_overlapping_idx
            
            % Integrate the wake deficit from turbine j over the rotor of turbine i
            sigma_j = sigma(:,j);
            r_dist_j = r_dist(j);
            r_dist_image_j = r_dist_image(j);            
            C_j = C(:,j).*is_overlapping(:,j);
            % (We only include the wake from turbine j for wind speeds where it overlaps with turbine i according to the
            % definition above.) 
            
            delta_real(:,j) = C_j.*overlap_gauss(r_dist_j./sigma_j, rotor_diameter(i)/2./sigma_j); % r_dist_j and rotor_diameter are normalized by sigma, according to the definitions in overlap_gauss              
            delta_image(:,j) = C_j.*overlap_gauss(r_dist_image_j./sigma_j, rotor_diameter(i)/2./sigma_j);            
        end
        
        % Include ground effects (i.e. image turbines) 
        delta = [delta_real, delta_image];
        
        % Katic superposition of wake deficits from upstream turbines
        delta_total(i,:) = sqrt(sum(delta.^2,2,'omitnan'));
    end
    
    % Resulting wind speed at turbine i given the combined deficit delta_total
    u(i,:) = u0.*(1-min(delta_total(i,:),1)).*u_corr(i);
    
    % Ct value given the resulting wind speed at turbine i
    ct(i,:) = power_curve(power_curve_index(i)).interpolant_ct(u(i,:));
    
    % Power given the resulting wind speed at turbine i
    power(i,:) = power_curve(power_curve_index(i)).interpolant_power(u(i,:));
    
end

% Reorder the calculated power and wind speed values to match original unsorted turbine coordinates
[~,orig_index] = sort(x_index);
power = power(orig_index,:);
u = u(orig_index,:);

end