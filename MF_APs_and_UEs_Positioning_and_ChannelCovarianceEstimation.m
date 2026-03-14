% This code places access points (APs) and user equipments (UEs) at random
% locations in the coverage area. Then, it computes the large-scale fading
% (channel gain) and channel covariance matrix

% Authors: Marx Freitas, André Fernandes, and Daynara Souza

% Set the distribution of the APs and UEs in the coverage area (XY plane)
distributionOfUEs = 'Uniform'; % Options: 'Gaussian', 'Uniform', 'FatherPoints'
distributionOfAPs = 'HCPP'; % Options: 'Gaussian', 'Uniform', 'HCPP'

% Coverage area dimensions in meters (defined as a rectangular area)
x_axis_size = 1000;
y_axis_size = 1000;

% Set the height of APs and UEs
HeightAP_m = 11.65;
HeightUE_m = 1.65;

% Propagation model parameters
% Options: 'Umi', 'InF-SL', 'InH-Mixed', 'InH-Open', 'InF-DL'
pathLossmodel  = 'Umi';

antennaSpacing = 1/2; % Half wavelength distance (in number of wavelengths)
ASD_deg = 15; % Angular standard deviation around the nominal angle (in degrees)
noiseFigure_dB = 8; % Noise figure (in dB) (from Fedex 05)

%=========================================================================%
% =========================== APs POSITIONING =========================== %
%=========================================================================%
% Places the APs in the coverage area, where area = x_axis_size*y_axis_size
APPositionsXYZ_m = MF_AF_GS_generateSetupAPs(x_axis_size, y_axis_size, ...
    HeightAP_m, nbrOfAPs, distributionOfAPs);

% Set the AP positions as a vector of complex numbers
APpositions_m = APPositionsXYZ_m(:,1) + 1i*APPositionsXYZ_m(:,2);

% Apply the wrap around technique to balance interference in distintc AP
% positions
[APpositionsWrapped, wrapLocations_m] = ...
    MF_generateSetupAPsWrapped(x_axis_size, y_axis_size, APpositions_m, ...
    nbrOfAPs);

%=========================================================================%
% =========================== UEs POSITIONING =========================== %
%=========================================================================%
% Places the APs in the coverage area, where area = x_axis_size*y_axis_size
UEPositionsXYZ_m = MF_AF_GS_generateSetupUEs(distributionOfUEs, ...
    x_axis_size, y_axis_size, nbrOfUEs, nbrOfAPs, HeightUE_m);

% Set the UE positions as a vector of complex numbers
UEpositions_m = UEPositionsXYZ_m(:,1) + 1i*UEPositionsXYZ_m(:,2);

%=========================================================================%
% ============== LARGE SCALE FADING and COVARIANCE MATRIX =============== %
%=========================================================================%
covariance_matrix = 'perfect'; % 'perfect', 'imperfect'

% For perfect and imperfect knowledge of covariance matrices
[gainOverNoise_dB, R, Hmean, HmeanFase, H, noiseVariance_dBm] = ...
    MF_CorrMatrix_Rician(pathLossmodel, nbrOfUEs, nbrOfAPs, N, ...
    APpositionsWrapped, UEpositions_m, HeightAP_m, HeightUE_m, BW_Hz, ...
    antennaSpacing, ASD_deg, fc_GHz, noiseFigure_dB, nbrOfRealizations, ...
    wrapLocations_m);

R_imp = R;

% For an imperfect knowledge of covariance matrices
if strcmp(covariance_matrix, 'imperfect')
    % It is recommended to use this function only for N = 1.
    % For accurate results, it may be necessary to consider the neta and
    % mi parameters for N > 1. It also works for N > 1, but the reference
    % paper recomends to regularize the covarance matrix. See the function
    % bellow for more details

    [gainOverNoise_dB, R_imp] = MF_CorrMatrix_Imperfect_Emil(nbrOfRealizations, ...
        nbrOfAPs, nbrOfUEs, N, poweUL_mW, tau_p, H, noiseVariance_dBm);
end