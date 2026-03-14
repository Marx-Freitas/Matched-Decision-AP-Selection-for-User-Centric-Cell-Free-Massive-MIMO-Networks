% This function computes the perfect large-scale fading of the UE

% Authors: Daynara Souza and Marx Freitas

% INPUT:
% pathLossmodel        = the path loss model utilized in the simulation
% nbrOfUEs             = the number of UEs in the newtork
% nbrOfAPs             = the number of APs in the newtork
% N                    = the number of antenna elements per AP
% APpositionsWrapped_m = the wrapped position of APs
%                        Dim: nbrOfAPs x 9
% UEpositions_m        = the position of APs
%                        Dim: nbrOfAPs x 1
% HeightAP_m           = the height of the AP in meters
% HeightUE_m           = the height of the UE in meters
% BW_Hz                = Bandwidth in Hz
% antennaSpacing       = space between antenna elements (in number of
%                        wavelengths)
% ASD_deg              = angular standard deviation (ASD) around the
%                        nominal angle (in degrees)
% fc_GHz               = the frequency center in GHz
% noiseFigure_dB       = Noise figure (in dB)
% nbrOfRealizations    = the number of realizations
% wrapLocations_m      = the center of each wrap arounded network
%                        Dim: 1 x 9

% OUTPUT:
% gainOverNoisedB   = the channel gain normalized by noise variance of UE k
%                     regarding all APs normalized by the noise variance.
%                     Dim: nbrOfAPs x 1
% R                 = the channel covariance matrix (correlation matrix)
%                     Dim: N x N x nbrOfAPs x nbrOfUEs
% Hmean             = the LOS component scaled by the Rician factor
%                     Dim: nbrOfAPs*N x nbrOfUEs
% HmeanFase         = the LOS component Hmean scaled by phase shifts
%                     Dim: nbrOfAPs*N x nbrOfRealizations x nbrOfUEs
% H                 = the real channel, comprising LOS and NLOS components
%                     Dim: nbrOfAPs*N x nbrOfRealizations x nbrOfUEs
% noiseVariance_dBm = the noise variance

% Important acronyms: line-of-sight (LOS), non-line-of-sight (NLOS)


function [gainOverNoise_dB, R, Hmean, HmeanFase, H, noiseVariance_dBm] = ...
    MF_CorrMatrix_Rician(pathLossmodel, nbrOfUEs, nbrOfAPs, N, ...
    APpositionsWrapped_m, UEpositions_m, HeightAP_m, HeightUE_m, BW_Hz, ...
    antennaSpacing, ASD_deg, fc_GHz, noiseFigure_dB, nbrOfRealizations, ...
    wrapLocations_m)

% Set the fading profile
fading_profile = 'Rayleigh'; % Options: 'Rician', 'Rayleigh'

% Standard deviation of the shadow fading in (5.43). Reference:
% Foundations of User-Centric Cell-Free Massive MIMO
sigma_sf = 4;

% Decorrelation distance of the shadow fading in (5.43). Reference:
% Foundations of User-Centric Cell-Free Massive MIMO
decorr = 9;

% Compute noise power
noiseVariance_dBm = -170 + 10*log10(BW_Hz) + noiseFigure_dB;

%% Prepare to store data
distance_AP_UE_3D_m = zeros(nbrOfAPs, nbrOfUEs);
gainOverNoise_dB    = zeros(nbrOfAPs, nbrOfUEs);
R                   = zeros(N, N, nbrOfAPs, nbrOfUEs); % Covariance matrix

% Prepare to store shadowing correlation matrix
shadowCorrMatrix     = sigma_sf^2*ones(nbrOfUEs, nbrOfUEs);
shadowAPrealizations = zeros(nbrOfUEs, nbrOfAPs);

% Prepare to store Rician factor matrices
LosVerifier  = zeros(nbrOfAPs, nbrOfUEs);
PrLOS        = zeros(nbrOfAPs, nbrOfUEs);
ricianFactor = zeros(nbrOfAPs, nbrOfUEs);

% Prepare to store the channel components, LOS and NLOS
Hmean     = zeros(nbrOfAPs*N, nbrOfUEs); % LOS component scaled by the Rician factor
HmeanFase = zeros(nbrOfAPs*N, nbrOfRealizations, nbrOfUEs); % Hmean scaled by phase shifts
H         = zeros(nbrOfAPs*N, nbrOfRealizations, nbrOfUEs); % Real channel, contains LOS and NLOS components

% Calculate the NLOS channel component before applying the channel
% covariance R and Rician factor.
Hnlos     = randn(nbrOfAPs*N, nbrOfRealizations,nbrOfUEs) + 1i*randn(nbrOfAPs*N, nbrOfRealizations, nbrOfUEs);

%% Processing

%=========================================================================%
% CHANNEL COVARIANCE MATRIX AND GENERATES THE REAL CHANNELS
%=========================================================================%

for indexUE = 1:nbrOfUEs

    % ---------------------------------------------------------------------
    % LARGE-SCALE FADING (CHANNEL GAIN) WITHOUT SHADOWING
    % ---------------------------------------------------------------------

    % Compute the 2D (X,Y) distancethe between the UE and APs.
    % NOTE:'min' is used because we adopted the wrap-around technique
    [distance_AP_UE_2D_m,whichpos] = min(abs(APpositionsWrapped_m - repmat(UEpositions_m(indexUE),size(APpositionsWrapped_m))),[],2);

    % Compute the 3D (X,Y) distancethe between the UE and APs
    % Distance between two points - linear algebra
    distance_AP_UE_3D_m(:,indexUE) = sqrt((HeightAP_m-HeightUE_m).^2 + distance_AP_UE_2D_m.^2) + eps;

    % Compute the channel gain divided by noise power

    if strcmp(fading_profile, 'Rician')
        [gainOverNoise_dB(:,indexUE), PrLOS(:,indexUE), LosVerifier(:,indexUE)] = AF_LargeScaleGainOverNoise(pathLossmodel, ...
            fc_GHz, noiseVariance_dBm, distance_AP_UE_3D_m(:,indexUE), HeightAP_m, HeightUE_m);
    
    elseif strcmp(fading_profile, 'Rayleigh')
    [gainOverNoise_dB(:,indexUE), ~, ~] = AF_LargeScaleGainOverNoise(pathLossmodel, ...
            fc_GHz, noiseVariance_dBm, distance_AP_UE_3D_m(:,indexUE), HeightAP_m, HeightUE_m);
    else
    % Do nothing
    end

    % ---------------------------------------------------------------------
    % SHADOWING
    % ---------------------------------------------------------------------
    % Calculate the conditional mean and standard deviation required to
    % obtain the new shadow fading realizations, when the shadow fading
    % realization of the previous UEs has already been generated. This
    % calculation is based on Theorem 10.2 in "Fundamentals of
    % Statistical Signal Processing: Estimation Theory" by S. Kay

    % If this is not the first UE
    if indexUE-1>0

        % Compute the distances of the new UE k regarding the other UEs
        shortestDistances = zeros(indexUE-1, 1);
        for i = 1:indexUE-1
            shortestDistances(i) = min(abs(UEpositions_m(indexUE) - UEpositions_m(i) + wrapLocations_m));
        end

        newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
        term1 = newcolumn'/shadowCorrMatrix(1:indexUE-1,1:indexUE-1);
        meanvalues = term1*shadowAPrealizations(1:indexUE-1,:);
        stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);

    else % If this is the first UE

        % Add the UE and begin to store shadow fading correlation values
        meanvalues = 0;
        stdvalue = sigma_sf;
        newcolumn = [];

    end

    % Generate the shadow fading realizations
    shadowing = meanvalues + stdvalue*randn(1,nbrOfAPs);

    % Update shadowing correlation matrix and store realizations
    shadowCorrMatrix(1:indexUE-1,indexUE) = newcolumn;
    shadowCorrMatrix(indexUE,1:indexUE-1) = newcolumn';
    shadowAPrealizations(indexUE,:)       = shadowing;

    % ---------------------------------------------------------------------
    % LARGE-SCALE FADING (CHANNEL GAIN) WITH SHADOWING
    % ---------------------------------------------------------------------
    gainOverNoise_dB(:,indexUE) = gainOverNoise_dB(:,indexUE) + transpose(shadowing);

    % ---------------------------------------------------------------------
    % RICIAN FACTOR
    % ---------------------------------------------------------------------
    ricianFactor(:,indexUE) = PrLOS(:,indexUE)./(1-PrLOS(:,indexUE));

    % ---------------------------------------------------------------------
    % COVARIANCE MATRIX (CORRELATION MATRIX)
    % ---------------------------------------------------------------------
    % Compute local scattering model of the spatial correlation matrix

    for indexAP = 1:nbrOfAPs
        %% Compute nominal angles between UE k and AP l
        angletoUE_varphi = angle(UEpositions_m(indexUE)-APpositionsWrapped_m(indexAP,whichpos(indexAP))); % azimuth angle
        angletoUE_theta  = asin((HeightAP_m-HeightUE_m)./distance_AP_UE_3D_m(indexAP,indexUE)); % elevation angle

        % Converting the angular standard deviations to radians
        ASD_varphi = ASD_deg*pi/180; % around the nominal azimuth angle
        ASD_theta  = ASD_deg*pi/180; % around the nominal elevation angle

        % Compute the local spatial covariance matrix
        R(:,:,indexAP,indexUE) = functionRlocalscatteringWithRician(N, ...
            angletoUE_varphi, angletoUE_theta, ASD_varphi, ASD_theta, ...
            antennaSpacing);

        % ---------------------------------------------------------------------
        % GENERATING THE REAL CHANNELS
        % ---------------------------------------------------------------------
        % Generate Rician channel realizations

        % Rayleigh fading channel, i.e., the LOS components are blocked
        if LosVerifier(indexAP,indexUE)==0

            % Compute the real channels
            H((indexAP-1)*N+1:indexAP*N,:,indexUE) = sqrt(0.5*db2pow(gainOverNoise_dB(indexAP,indexUE)))*sqrtm(R(:,:,indexAP,indexUE))*Hnlos((indexAP-1)*N+1:indexAP*N,:,indexUE);

            % Recompute the coraviance matrix
            R(:,:,indexAP,indexUE) = db2pow(gainOverNoise_dB(indexAP,indexUE))*R(:,:,indexAP,indexUE);

            % Rician fading channel, i.e., the LOS components are not blocked
        else
            if ricianFactor(indexAP,indexUE) == Inf
                factorNLoS = 0; % There are no NLOS components
                factorLoS  = 1; % Only LOS components
            else
                % Compute the contribution of LOS and NLOS components to
                % the channel vector
                factorNLoS = 1/(ricianFactor(indexAP,indexUE)+1);
                factorLoS  = ricianFactor(indexAP,indexUE)/(ricianFactor(indexAP,indexUE)+1);
            end

            % Apply correlation to the rician channel realizations
            % LOS channel using 3D ULA response
            dn = linspace(1, N, N)';

            % Compute the LOS component scaled by the Rician factor
            Hmean((indexAP-1)*N+1:indexAP*N,indexUE) = (sqrt(factorLoS*db2pow(gainOverNoise_dB(indexAP,indexUE)))*exp(-1j*pi*sin(angletoUE_varphi)*cos(angletoUE_theta).*(dn-1)));

            % Compute Hmean scaled by phase shifts
            HmeanFase((indexAP-1)*N+1:indexAP*N,:,indexUE) = Hmean((indexAP-1)*N+1:indexAP*N,indexUE)*exp(1j*rand(1,nbrOfRealizations)*2*pi);

            % Compute the real channels
            H((indexAP-1)*N+1:indexAP*N,:,indexUE) = HmeanFase((indexAP-1)*N+1:indexAP*N,:,indexUE) ...
                + sqrt(0.5*factorNLoS*db2pow(gainOverNoise_dB(indexAP,indexUE)))*sqrtm(R(:,:,indexAP,indexUE))*Hnlos((indexAP-1)*N+1:indexAP*N,:,indexUE);

            % Recompute the coraviance matrix
            R(:,:,indexAP,indexUE) = factorNLoS*db2pow(gainOverNoise_dB(indexAP,indexUE))*R(:,:,indexAP,indexUE);
        end

    end

    end