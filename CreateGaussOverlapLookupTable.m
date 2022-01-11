function overlap_lookup_table = CreateGaussOverlapLookupTable
%CREATEGAUSSOVERLAPLOOKUPTABLE creates the look-up table used in the TurbOPark wake model when calculating the
% rotor-averaged Gaussian deficit. More specifically the function numerically calculates the integral of
% (1/A)*exp(-rw^2/(2*sigma)) over the rotor being waked. A is the rotor area, rw is the radial distance from the
% wake centreline, and sigma is the charactestic wake width.
%
% Inputs: none
%
% Outputs:
%   gauss_lookup_table (struct):
%       - dist: array, radial distance between the upstream and downstream rotor centers, normalized to the wake width
%       - radius_down: array, radius of downstream turbine, normalized to the wake width
%       - overlap_gauss: matrix of size [length(dist),length(radius_down)], the integral evaluated at dist and
%                        radius_down 
%       NB: The struct is also saved to the folder of this function as gauss_lookup_table.mat
%
% Examples:
%   ::
%
%       gauss_lookup_table = CreateGaussOverlapLookupTable;
%
% References:
%   - https://github.com/OrstedRD/TurbOPark/blob/main/TurbOPark%20description.pdf
%
% See also:
%   - TurbOPark
%
%   Author: Benny Lassen
%   Checked by: Sidse Damgaard Hansen, Jesper Grønnegaard
%
%   Copyright (c) 2021 by Ørsted

% Radial distance between the upstream and downstream rotor centers, normalized to the characteristic wake width at the
% downstream rotor: 
dist = 0:0.01:10; 

% Radius of downstream rotor, normalized to the characteristic wake width at the downstream rotor:
radius_down = 0:0.01:20; 

% The integral of the gaussian wake over the downstream rotor area for combinations of dist and radius_down
overlap_gauss = zeros(length(dist),length(radius_down)); 

for i = 1:length(dist) % Loop over distance
    for j = 1:length(radius_down) % Loop over rotor radius       
        if radius_down(j) > 0
            % Define the integral (without the 1/area)
            fun = @(r,theta) exp(-(r.^2+dist(i)^2-2*dist(i)*r.*cos(theta))/2).*r;
            % Do the integral numerically
            out = integral2(fun, 0, radius_down(j), 0, 2*pi, 'Method', 'iterated', 'AbsTol', 0, 'RelTol',1e-10);
            % Include the 1/area
            out = out/(pi*radius_down(j)^2);
        else % In the limit of radius_down -> 0
            out = exp(-(dist(i)^2)/2); % Derived using L'Hôpital's rule and Leibniz integral rule
        end
        overlap_gauss(i,j) = out;
    end
end

% Save the result to a mat-file as a struct
overlap_lookup_table = struct('dist', dist,...
    'radius_down', radius_down,...
    'overlap_gauss', overlap_gauss);

save('gauss_lookup_table.mat', 'overlap_lookup_table');

end
