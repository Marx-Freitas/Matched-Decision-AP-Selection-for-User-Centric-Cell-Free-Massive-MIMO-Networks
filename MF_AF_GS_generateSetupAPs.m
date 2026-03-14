% This function performs the positioning of access points (APs) in the
% coverage area. The APs can be distributed in the coverage area following
% conventional geometries (such as a matrix, i.e., a grid) or following a
% random distribution. Please, note that this function assumes that the APs
% are placed in a rectangular area, where area = x_axis_size*y_axis_size

% Authors: Marx Freitas and Gilvan Soares Borges

% INPUT:
% distributionOfAPs = the distribution of APs in the coverage area (XY plane)
% x_axis_size       = the length in meters of the X-axis of the rectangular coverage area
% y_axis_size       = the length in meters of the Y-axis of the rectangular coverage area
% nbrOfAPs          = the number of APs in the newtork
% HeightAP_m        = the height of each AP in meters

% OUTPUT:
% APPositionsXYZ_m  = the position of the APs
%                     Dim: nbrOfAPs x 3 (i.e., a 3D dimension)

function APPositionsXYZ_m = MF_AF_GS_generateSetupAPs(x_axis_size, y_axis_size, HeightAP_m, nbrOfAPs, distributionOfAPs)

switch distributionOfAPs

    case 'Gaussian'
        % Author: Marx M. M. Freitas
        %=====================================================================%
        % When the distribution of the APs follows a Gaussian distribution
        %=====================================================================%
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = [ x_axis_size * randn(nbrOfAPs, 1)...
            y_axis_size * randn(nbrOfAPs, 1) HeightAP_m * ones(nbrOfAPs, 1) ];

    case 'Uniform'
        % Author: Marx M. M. Freitas
        %=====================================================================%
        % When the distribution of the APs follows a uniform distribution
        %=====================================================================%
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = [ -x_axis_size+(2*x_axis_size)*rand(nbrOfAPs, 1)...
            -y_axis_size+(2*y_axis_size)*rand(nbrOfAPs, 1) ...
            HeightAP_m * ones(nbrOfAPs, 1) ];

    case 'HCPP'
        % Author: Gilvan Soares Borges (a very intelligent man)
        %=====================================================================%
        % When the distribution of the APs follows a hard core point
        % poission process (HCPP)
        %=====================================================================%
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = GS_apsDistribution(nbrOfAPs,[x_axis_size; y_axis_size], HeightAP_m, 'regularity', 1);

        % This code line fixes a very unlikely numerical error
        while(length(APPositionsXYZ_m)> nbrOfAPs || length(APPositionsXYZ_m) < nbrOfAPs)
            APPositionsXYZ_m = GS_apsDistribution(nbrOfAPs,[x_axis_size; y_axis_size], HeightAP_m, 'regularity', 1);
        end
        
        % This code line adjusts the AP positions
        APPositionsXYZ_m(:,1) = (2.*APPositionsXYZ_m(:,1)) - x_axis_size;
        APPositionsXYZ_m(:,2) = (2.*APPositionsXYZ_m(:,2)) - y_axis_size;

        clear APPositionsXYZdim1_m APPositionsXYZdim2_m APPositionsXYZdim3_m APPositionsXYZdim4_m...
            APPositionsRecXYZ_m APPositionsXYZXdiag_m indexXdiag numberOfAPsXdiag

end

% Adjusting the APs locations at XY plane. By doing so, the coverage area
% will be contained in (x=0:x_axis_size, y=0:y_axis_size)
APPositionsXYZ_m(:,1) = (APPositionsXYZ_m(:,1) + x_axis_size)/2;
APPositionsXYZ_m(:,2) = (APPositionsXYZ_m(:,2) + y_axis_size)/2;

end