clear; clc; tic
%=========================================================================%
%% Input parametes
%=========================================================================%
% Choose the access point (AP) selection strategy utilized to associate the
% user-equipment (UE) with a subset of APs.
APselectionMethod = 'matchedDecision'; % Options: 'CanonicalCF', 'UCC', 'LSFB', 'ScalableCF', 'matchedDecision', 'MD_LSFB'

% Select the AP cluster adjustment strategy 
fineTunning = 'No_FineTuning'; % Options: 'No_FineTuning', 'Algorithm_3', 'Algorithm_4', 'Algorithm_5'

% Set the number of APs (L) and UEs (K) in the coverage area
nbrOfAPs = 100;
nbrOfUEs = 25;

% Set the number of antenna elements per AP (N)
N = 1;

% Set the maximum number of UEs that each AP can serve
U_max = 4;

% Set the number of monte-Carlo setups and channel realizations
nbrOfSetups = 400; nbrOfRealizations = 400;

% Set the number of samples per coherence block and for pilot signaling
tau_c = 200; tau_p = 10;

% Transmission power of uplink (UL) and downlink (DL)
poweUL_mW = 100; PowerDL_mW = 1000;

% Propagation model parameters
BW_Hz = 100e6; % Bandwidth in Hz
fc_GHz = 3.5; % Frequency center

%%
%=========================================================================%
% Processing
%=========================================================================%
% Preparing to store results
MF_PreparingToStoreResults

% Processing capacity limitation, in the access points (APs)
U_max = min(U_max, tau_p); % If tau_p < U_max, we set U_max to be equal to tau_p

% Go through all setups
for nSetup = 1:nbrOfSetups

    % Display simulation progress
    disp(['Setup ' num2str(nSetup) ' out of ' num2str(nbrOfSetups)]);

    % Places the APs and UEs in the coverage area. Then, it calculates the
    % large-scale fading (channel gain) and computes the channel covariance
    % matrix
    MF_APs_and_UEs_Positioning_and_ChannelCovarianceEstimation

    % Performs AP clustering and channel estimation
    MF_APselection_APclusterAdjustment_and_ChannelEstimation

    % Computes the computational complexity (CC), spectral efficiency (SE),
    % and energy efficiency (EE). The code below also calculates the number
    % of APs serving each UE (|M_{k}|) and the number of UEs served by each
    % AP (|D_{l}|)
    MF_Compute_SE_EE_APsperUE_and_UEsperAP

end

% Plot the achieved results
MF_PlotAndSaveResults

toc

