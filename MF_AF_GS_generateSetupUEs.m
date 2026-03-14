% This function performs the positioning of user-equipment (UEs) in the
% coverage area. The UEs can be distributed in the coverage area following
% conventional geometries (such as a matrix, grid) or following a random
% distribution. Please, note that the function assumes that the UEs are
% placed in a rectangular area, where area = x_axis_size*y_axis_size

% Author: Marx Freitas

% INPUT:
% distributionOfUEs = the distribution of UEs in the coverage area (XY plane)
% x_axis_size       = the length in meters of the X-axis of the rectangular coverage area 
% y_axis_size       = the length in meters of the Y-axis of the rectangular coverage area
% nbrOfUEs          = the number of UEs in the newtork
% nbrOfAPs          = the number of APs in the newtork
% HeightUE_m        = the height of the UE in meters

% OUTPUT:
% UEPositionsXYZ_m  = the position of UEs
%                     Dim: nbrOfUEs x 3 (i.e., a 3D dimension) 

function UEPositionsXYZ_m = MF_AF_GS_generateSetupUEs(distributionOfUEs, ...
    x_axis_size, y_axis_size, nbrOfUEs, nbrOfAPs, HeightUE_m)

% Deploy the UEs positions in the coverage area
switch distributionOfUEs
        
    case ('Gaussian')
        % Author: Marx Freitas
        %=====================================================================%
        % When the distribution of the UEs follows a Gaussian distribution
        %=====================================================================%
        
        % Set the UEs locations at XYZ plane
        UEPositionsXYZ_m       = zeros(nbrOfUEs, 3);
        UEPositionsXYZ_m(:, 1) = (x_axis_size/2).*randn(1, nbrOfUEs);
        UEPositionsXYZ_m(:, 2) = (y_axis_size/2).*randn(1, nbrOfUEs);
        UEPositionsXYZ_m(:, 3) = HeightUE_m;
        
        % In order to ensure the UEs remain in the area of the simulated setup
        UEPositionsXYZ_m(abs(UEPositionsXYZ_m(:,1))>x_axis_size,1) = x_axis_size;
        UEPositionsXYZ_m(abs(UEPositionsXYZ_m(:,2))>y_axis_size,2) = y_axis_size;
        
    case ('Uniform')
        % Author: Marx Freitas
        %=====================================================================%
        % When the distribution of UEs APs follows a uniform distribution
        %=====================================================================%
        
        % Set the UEs locations at XYZ plane
        UEPositionsXYZ_m       = zeros(nbrOfUEs, 3);
        UEPositionsXYZ_m(:, 1) = (2*x_axis_size).*rand(1, nbrOfUEs)-x_axis_size;
        UEPositionsXYZ_m(:, 2) = (2*y_axis_size).*rand(1, nbrOfUEs)-y_axis_size;
        UEPositionsXYZ_m(:, 3) = HeightUE_m;
        
    case ('FatherPoints')
        % Author: Marx M.M. Freitas (experimental, not used in simulations)
        %=====================================================================%
        % When UEs are clustered around a certain point (father point)
        %=====================================================================%
        UEsPerCluster = 25;
        APsUEsRatioStandard = 4;        
        APsUEsRatio  = floor(nbrOfAPs/nbrOfUEs);
        UEsPerCluster = min(UEsPerCluster*ceil(APsUEsRatio/APsUEsRatioStandard),nbrOfUEs);
        nbrOfClusters = ceil(nbrOfUEs/UEsPerCluster);        
         
        % Set the UEs locations at XYZ plane
        LengthPercent = 0.85;
        FatherPoints_m = GS_apsDistribution(nbrOfClusters,[LengthPercent*x_axis_size; LengthPercent*y_axis_size], 'regularity', 1);
        FatherPoints_m(:,1) = (2.*FatherPoints_m(:,1)) - x_axis_size;
        FatherPoints_m(:,2) = (2.*FatherPoints_m(:,2)) - y_axis_size;
        
        cluesterLocations = FatherPoints_m(:,1) + 1i*FatherPoints_m(:,2);                
        clusterPercentage = 0.1;        
        nbrUEsPerCluster = zeros(nbrOfClusters,1);
        auxUEposition    = zeros(nbrOfClusters+1,1);
        
        for indexCluster = 1:nbrOfClusters
            
            nbrUEsPerCluster(indexCluster,1) = UEsPerCluster;
            auxUEposition(indexCluster+1,1) = sum(nbrUEsPerCluster(1:indexCluster,1));
            
            indexCLocation = randi(length(cluesterLocations));
                        
            UEPositionsXYZ_m(1+auxUEposition(indexCluster,1):auxUEposition(indexCluster+1),1) = ...
                real(cluesterLocations(indexCLocation)) + ...
                (clusterPercentage*x_axis_size)*randn(nbrUEsPerCluster(indexCluster,1),1);
            
            UEPositionsXYZ_m(1+auxUEposition(indexCluster,1):auxUEposition(indexCluster+1),2) = ...
                imag(cluesterLocations(indexCLocation)) + ...
                (clusterPercentage*y_axis_size)*randn(nbrUEsPerCluster(indexCluster,1),1);
            
            cluesterLocations(indexCLocation) = [];
            
            remainingUEs = nbrOfUEs - (UEsPerCluster)*indexCluster;
            UEsPerCluster = min(UEsPerCluster,remainingUEs);
            UEsPerCluster(UEsPerCluster<0) = 0;
            
        end        
        
        UEPositionsXYZ_m(:, 3) = HeightUE_m;
                
end

UEPositionsXYZ_m(:,1) = (UEPositionsXYZ_m(:,1) + x_axis_size)/2;
UEPositionsXYZ_m(:,2) = (UEPositionsXYZ_m(:,2) + y_axis_size)/2;

end