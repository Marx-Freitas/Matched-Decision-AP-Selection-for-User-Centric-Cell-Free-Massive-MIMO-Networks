% This function implements the path loss models defined in the 3GPP
% technical report (TR) TR 138 901 V16.1.0 (2020-11), named "5G; Study on
% channel model for frequencies from 0.5 to 100 GHz (3GPP TR 38.901
% version 16.1.0 Release 16)"

% Authors: André Fernandes and Marx Freitas

% This function implements the follwing path-loss models
% OUTDOOR:
% Urban Micro (UMi): UMi-Street Canyon

% INDOOR: 
% Indoor Factory (InF): 'InF-SL' and 'InF-DL'
% Indoor Hotspot (InH): 'InH-Open' (Indoor - Open office) and 'InH-Mixed'
% (Indoor - Mixed office)

% Note: important informations and descriptions can be found throughout the
% code. For instance, the code mentions the tables and parameters from 3GPP
% TR 38.901 we utilized for writting this code

% INPUT:
% pathLossmodel    = the path loss model utilized in the simulation
% fc_GHz           = the frequency center in GHz
% noiseVariancedBm = the noise variance in dBm
% distance3D_m     = the Euclidean 3D distance (x, y, z) between the UE and
%                    and APs in meters.
%                    Dim: L x 1 (i.e., nbrOfAPs x 1)
% h_AP_m           = the height of the AP in meters
% h_UE_m           = the height of the UE in meters

% OUTPUT:
% gainOverNoisedB   = the channel gain normalized by noise variance of UE k
%                     regarding all APs normalized by the noise variance.
%                     Dim: L x 1 (i.e., nbrOfAPs x 1)

% Pr_LOS            = the line-of-sight (LOS) probability of UE k regarding
%                     all APs. Dim: L x 1 (i.e., nbrOfAPs x 1)

% selector_LOS_NLOS = indicates if the channel between the UE k and AP l is
%                     LOS or non-line-of-sight (NLOS).
%                     Dim: L x 1 (i.e., nbrOfAPs x 1)

function [gainOverNoisedB, Pr_LOS, selector_LOS_NLOS] = AF_LargeScaleGainOverNoise(pathLossmodel, fc_GHz, noiseVariancedBm, distance3D_m, h_AP_m, h_UE_m)

heigth_diference=h_AP_m-h_UE_m;

if strcmp(pathLossmodel,'Umi')==1
    % Description: UMi (Street canyon, open area) with O2O and O2I: This is
    % similar to 3D-UMi scenario, where the BSs are mounted below rooftop
    % levels of surrounding buildings. UMi open area is intended to capture
    % real-life scenarios such as a city or station square. The width of
    % the typical open area is in the order of 50 to 100 m.
    % Example: [Tx height:10m, Rx height: 1.5-2.5 m, ISD: 200m]

    % Prepare to store LOS probability (Pr_LOS) and path loss (PL)
    Pr_LOS = zeros(size(distance3D_m));
    PL     = zeros(size(distance3D_m));
    
    % Compute the breakpoint distance. See note 1 (page 29) of TR 138 901 
    he_m = 1; % Effective environment height 
    h_AP_line_m = h_AP_m - he_m; % Effective access point (AP) height
    h_UE_line_m = h_UE_m - he_m; % Effective user-equipment (UE) height

    % Breakpoint distance, see note 1, page 29 of TR 138 901 V16.1.0
    d_bp_line = 4*h_AP_line_m*h_UE_line_m*fc_GHz*10^9/(300*10^6);
    
    % Compute distance2D (d_2D) and selector_LOS_NLOS
    distance2D = sqrt(distance3D_m.^2-heigth_diference^2);
    selector_LOS_NLOS = rand(size(distance3D_m));
    
    % Computing the LOS probability (Pr_LOS), see Table 7.4.2-1
    Pr_LOS(distance2D<=18) = 1;
    if sum(Pr_LOS) ~= length(distance2D)
        Pr_LOS(distance2D>18) = ones(size(distance2D(distance2D>18)))*18./distance2D(distance2D>18) +...
            exp(-distance2D(distance2D>18)/36).*(1-(ones(size(distance2D(distance2D>18)))*18./distance2D(distance2D>18)));
    end
    
    % Set selector_LOS_NLOS to 1 for LOS, and 0 for NLOS
    selector_LOS_NLOS = (Pr_LOS>selector_LOS_NLOS);
    
    % Compute the LOS path loss (PL), see Table 7.4.1-1
    PL((selector_LOS_NLOS==1) == (distance2D<=d_bp_line)) = 32.4+21*log10(distance3D_m((selector_LOS_NLOS==1)==(distance2D<=d_bp_line)))+20*log10(fc_GHz);
    PL((selector_LOS_NLOS==1) == (distance2D>d_bp_line))  = 32.4+40*log10(distance3D_m((selector_LOS_NLOS==1)==(distance2D>d_bp_line))) +20*log10(fc_GHz)...
        -9.5*log10(d_bp_line^2+heigth_diference^2);
    
    % Computing the NLOS path loss (PL), see Table 7.4.1-1
    % Auxiliary variables to facilitate calculation of NLOS path loss
    aux=zeros(size(distance3D_m));
    aux((selector_LOS_NLOS==0) == (distance2D<=d_bp_line)) = 32.4+21*log10(distance3D_m((selector_LOS_NLOS==0)==(distance2D<=d_bp_line)))+20*log10(fc_GHz);
    aux((selector_LOS_NLOS==0) == (distance2D>d_bp_line))  = 32.4+40*log10(distance3D_m((selector_LOS_NLOS==0)==(distance2D>d_bp_line))) +20*log10(fc_GHz)...
        -9.5*log10(d_bp_line^2+heigth_diference^2);
    aux2=zeros(size(distance3D_m));
    aux2((selector_LOS_NLOS==0)) = 35.3*log10(distance3D_m((selector_LOS_NLOS==0))) + 22.4 +21.3*log10(fc_GHz) -0.3*(h_UE_m-1.5);
    
    PL((selector_LOS_NLOS==0)) = max(aux(selector_LOS_NLOS==0),aux2(selector_LOS_NLOS==0));
    % QUESTION: Should shadowing be 7.82 for NLOS?
    
    % Compute the channel gain normalized by the noise variance 
    gainOverNoisedB = -PL-noiseVariancedBm;

elseif strcmp(pathLossmodel,'InH-Mixed')==1
    % Description: 'InH-Mixed' means Indoor - Mixed office. This scenario
    % represents a mixed office environment where there might be a
    % combination of open areas and more obstructed spaces, such as walls
    % and partitions.
    
    % Details on indoor-office scenarios are listed in Table 7.2-2 and
    % presented in Figure 7.2-1. More details, if necessary, can be added
    % to Figure 7.2-1. The only difference between the open office and
    % mixed office models in this TR is the line of sight probability
    
    % Prepare to store LOS probability (Pr_LOS) and path loss (PL)
    Pr_LOS = zeros(size(distance3D_m));
    PL     = zeros(size(distance3D_m));
    
    % Compute distance2D (d_2D) and selector_LOS_NLOS
    distance2D = sqrt(distance3D_m.^2-heigth_diference^2);
    selector_LOS_NLOS = rand(size(distance3D_m));

    % Computing the LOS probability (Pr_LOS), see Table 7.4.2-1
    Pr_LOS(distance2D<=1.2) = 1;
    Pr_LOS((distance2D<=6.5) == (distance2D>1.2)) = exp(-(distance2D((distance2D<=6.5)==(distance2D>1.2))-1.2)/4.7);
    Pr_LOS(distance2D>6.5) = exp(-(distance2D(distance2D>6.5)-6.5)/32.6)*0.32;
    
    % Set selector_LOS_NLOS to 1 for LOS, and 0 for NLOS
    selector_LOS_NLOS=(Pr_LOS>selector_LOS_NLOS);
    
    % Compute the LOS path loss (PL), see Table 7.4.1-1 (InH - Office)
    PL(selector_LOS_NLOS==1) = 32.4+17.3*log10(distance3D_m(selector_LOS_NLOS==1)) + 20*log10(fc_GHz);
    
    % Compute the NLOS path loss (PL), see Table 7.4.1-1 (InH - Office)
    PL(selector_LOS_NLOS==0) = max(38.3*log10(distance3D_m(selector_LOS_NLOS==0)) + 17.3 + 24.9*log10(fc_GHz)...
        ,32.4+17.3*log10(distance3D_m(selector_LOS_NLOS==0))+20*log10(fc_GHz));
    
    % Compute the channel gain normalized by the noise variance
    gainOverNoisedB=-PL-noiseVariancedBm;

elseif strcmp(pathLossmodel,'InH-Open')==1
    % Description: 'InH-Open' means Indoor - Open office.
    % This scenario models an open office indoor environment where there is
    % relatively more open space without dense obstacles.
    
    % Details on indoor-office scenarios are listed in Table 7.2-2 and
    % presented in Figure 7.2-1. More details, if necessary, can be added
    % to Figure 7.2-1. The only difference between the open office and
    % mixed office models in this TR is the line of sight probability

    % Prepare to store LOS probability (Pr_LOS) and path loss (PL)
    Pr_LOS = zeros(size(distance3D_m));
    PL     = zeros(size(distance3D_m));
    
    % Compute distance2D (d_2D) and selector_LOS_NLOS
    distance2D = sqrt(distance3D_m.^2-heigth_diference^2);
    selector_LOS_NLOS = rand(size(distance3D_m));

    % Computing the LOS probability (Pr_LOS), see Table 7.4.2-1
    Pr_LOS(distance2D<=5) = 1;
    Pr_LOS((distance2D<=49) == (distance2D>5)) = exp(-(distance2D((distance2D<=49)==(distance2D>5))-5)/70.8);
    Pr_LOS(distance2D>49) = exp(-(distance2D(distance2D>49)-49)/211.17)*0.54;
    
    % Set selector_LOS_NLOS to 1 for LOS, and 0 for NLOS
    selector_LOS_NLOS=(Pr_LOS>selector_LOS_NLOS);
    
    % Compute the LOS path loss (PL), see Table 7.4.1-1 (InH - Office)
    PL(selector_LOS_NLOS==1) = 32.4+17.3*log10(distance3D_m(selector_LOS_NLOS==1))+20*log10(fc_GHz);

    % Compute the NLOS path loss (PL), see Table 7.4.1-1 (InH - Office)
    PL(selector_LOS_NLOS==0) = max(38.3*log10(distance3D_m(selector_LOS_NLOS==0)) + 17.3 + 24.9*log10(fc_GHz)...
        ,32.4+17.3*log10(distance3D_m(selector_LOS_NLOS==0))+20*log10(fc_GHz));
    
    % Compute the channel gain normalized by the noise variance
    gainOverNoisedB = -PL-noiseVariancedBm;

elseif strcmp(pathLossmodel,'InF-SL')==1
    % The indoor factory (InF) scenario focuses on factory halls of varying
    % sizes and with varying levels of density of "clutter", e.g.
    % machinery, assembly lines, storage shelves, etc. Details of the InF
    % scenario are listed in Table 7.2-4.

    % Description: 'InF-SL' means indoor Factory with sparse clutter and
    % low base station heigh (both Tx and Rx are below the average height
    % of the clutter average height of the clutter)
        
    % Set clutter size (d_clutter_m) and clutter density (clutter_density)
    % See Table 7.2-4 for setting d_clutter_m and clutter_density
    d_clutter_m = 10;
    clutter_density = 0.2;
    
    % Compute distance2D (d_2D) and selector_LOS_NLOS
    distance2D=sqrt(distance3D_m.^2-heigth_diference^2);
    selector_LOS_NLOS=rand(size(distance3D_m));

    % Computing the LOS probability (Pr_LOS), see Table 7.4.2-1
    Pr_LOS = exp(-distance2D/(-d_clutter_m/log(1-clutter_density)));
    
    % Set selector_LOS_NLOS to 1 for LOS, and 0 for NLOS
    selector_LOS_NLOS = (Pr_LOS>selector_LOS_NLOS);
    
    % Compute the LOS path loss (PL), see Table 7.4.1-1
    PL(selector_LOS_NLOS==1) = 31.84 + 21.5*log10(distance3D_m(selector_LOS_NLOS==1))+19*log10(fc_GHz);

    % Compute the NLOS path loss (PL), see Table 7.4.1-1
    PL(selector_LOS_NLOS==0)=max(33+25.5*log10(distance3D_m(selector_LOS_NLOS==0))+20*log10(fc_GHz)...
        ,31.84+25.5*log10(distance3D_m(selector_LOS_NLOS==0))+19*log10(fc_GHz));
    
    % Compute the channel gain normalized by the noise variance
    gainOverNoisedB = -PL-noiseVariancedBm;

elseif strcmp(pathLossmodel,'InF-DL')==1
    % The indoor factory (InF) scenario focuses on factory halls of varying
    % sizes and with varying levels of density of "clutter", e.g.
    % machinery, assembly lines, storage shelves, etc. Details of the InF
    % scenario are listed in Table 7.2-4.

    % Description: 'InF-DL' means Indoor Factory with Dense clutter and Low
    % base station height (both Tx and Rx are below the average height of
    % the clutter)
    
    % Set clutter size (d_clutter_m) and clutter density (clutter_density)
    % See Table 7.2-4 for setting d_clutter_m and clutter_density
    d_clutter_m = 2;
    clutter_density = 0.7;
    
    % Compute distance2D (d_2D) and selector_LOS_NLOS
    distance2D=sqrt(distance3D_m.^2-heigth_diference^2);
    selector_LOS_NLOS=rand(size(distance3D_m));

    % Computing the LOS probability (Pr_LOS), see Table 7.4.2-1
    Pr_LOS = exp(-distance2D/(-d_clutter_m/log(1-clutter_density)));
    
    % Set selector_LOS_NLOS to 1 for LOS, and 0 for NLOS
    selector_LOS_NLOS = (Pr_LOS>selector_LOS_NLOS);
    
    % Compute the LOS path loss (PL), see Table 7.4.1-1
    PL(selector_LOS_NLOS==1) = 31.84+21.5*log10(distance3D_m(selector_LOS_NLOS==1))+19*log10(fc_GHz);
    
     % Compute the NLOS path loss (PL), see Table 7.4.1-1
    PL(selector_LOS_NLOS==0) = max(18.6+35.7*log10(distance3D_m(selector_LOS_NLOS==0))+20*log10(fc_GHz)...
        ,31.84+25.5*log10(distance3D_m(selector_LOS_NLOS==0))+19*log10(fc_GHz));
    
    % Compute the channel gain normalized by the noise variance
    gainOverNoisedB = -PL-noiseVariancedBm;

else
    
end

end