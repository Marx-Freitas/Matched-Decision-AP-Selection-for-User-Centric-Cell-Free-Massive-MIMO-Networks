% This code computes some performance metrics associated with user-centric
% cell-free massive MIMO systems. Specifically, it calculates the spectral
% efficiency (SE), and energy efficiency (EE). This code also calculates
% the number of APs per UE (|M_{k}|) and UEs per AP (|U_{l}|)

% Authors: Marx M. M. Freitas, André L. P. Fernandes, and Daynara D. Souza

%=========================================================================%
% ================== DOWNLINK SPECTRAL EFFICIENCY (SE) ================== %
%=========================================================================%

% ----------------------------------------------------------------------- %
% Compute the DL power allocation for distributed implementation
% ----------------------------------------------------------------------- %
% Downlink power allocation for MR and LP-MMSE precoding schemes
rho_dist_mW = DLpowerAllocation(nbrOfUEs, nbrOfAPs, D, PowerDL_mW,...
    gainOverNoise_dB);

% ----------------------------------------------------------------------- %
% Compute the DL power allocation for centralized implementation
% ----------------------------------------------------------------------- %
% The downlink power allocation for P-RZF and P-MMSE precoding schemes
% is performed inside the function "functionComputeSE_downlinkWithRician"

% ----------------------------------------------------------------------- %
% Compute SE using the capacity bound in Proposition 3
% ----------------------------------------------------------------------- %
% For perfect and imperfect knowledge of covariance matrices
[SE_P_MMSE_tot(:,nSetup), SE_P_RZF_tot(:,nSetup), SE_LP_MMSE_tot(:,nSetup), ...
    SE_MR_tot(:,nSetup), w_MR, w_LP_MMSE, ~, ~]...
    = functionComputeSE_downlinkWithRician(Hhat, H, D, C, tau_c, tau_p, ...
    nbrOfRealizations, N, nbrOfUEs, nbrOfAPs, poweUL_mW, rho_dist_mW, ...
    gainOverNoise_dB, PowerDL_mW);

% For an imperfect knowledge of covariance matrices
if strcmp(covariance_matrix, 'imperfect')
    Ts = 0.5; ts = BW_Hz*Ts/tau_c;
    alfa = (tau_p*nbrOfRealizations/2)/(ts*tau_c);
    oldPrelogFactor = (tau_c-tau_p)/tau_c;
    newPrelogFactor = (1-tau_p/tau_c-alfa);
    adjustedPrelogFactor = newPrelogFactor/oldPrelogFactor;

    SE_P_MMSE_tot(:,nSetup)  = SE_P_MMSE_tot(:,nSetup)*adjustedPrelogFactor;
    SE_LP_MMSE_tot(:,nSetup) = SE_LP_MMSE_tot(:,nSetup)*adjustedPrelogFactor;
    SE_P_RZF_tot(:,nSetup)   = SE_P_RZF_tot(:,nSetup)*adjustedPrelogFactor;
    SE_MR_tot(:,nSetup)      = SE_MR_tot(:,nSetup)*adjustedPrelogFactor;
end

%=========================================================================%
% ======================= ENERGY EFFICIENCY (EE) ======================== %
%=========================================================================%
% Compute energy efficiency
Ee_MR(nSetup,1)      = MF_EnergyEfficiency (BW_Hz, D, SE_MR_tot(:,nSetup),...
    nbrOfAPs, N, w_MR, nbrOfRealizations);
Ee_LP_MMSE(nSetup,1) = MF_EnergyEfficiency (BW_Hz, D, SE_LP_MMSE_tot(:,nSetup),...
    nbrOfAPs, N, w_LP_MMSE, nbrOfRealizations);


%=========================================================================%
% =========== APs PER UE (|M_{k}|), AND UEs PER AP (|D_{l}|) ============ %
%=========================================================================%
%% Number of active APs, APs connected to each UE, and UEs connected to each AP
nbrUEsAPisServing = sum(transpose(D(1:end,:)));
ActiveAPs(nSetup,1) = nbrOfAPs - sum(nbrUEsAPisServing == 0);
APsConnectedToUE(nSetup,:) = sum(D(:,1:end));
UEsConnectedperAP(nSetup,:) = sum(D,2);