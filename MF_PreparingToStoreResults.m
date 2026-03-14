% Prepare to store results
% Author: Marx M. M. Freitas

%=========================================================================%
% ====================== SPECTRAL EFFICIENCY (EE) ======================= %
%=========================================================================%
% Prepare to save downlink simulation results
SE_LP_MMSE_tot = zeros(nbrOfUEs, nbrOfSetups);
SE_P_MMSE_tot  = zeros(nbrOfUEs, nbrOfSetups);
SE_MR_tot      = zeros(nbrOfUEs, nbrOfSetups);
SE_P_RZF_tot   = zeros(nbrOfUEs, nbrOfSetups);

%=========================================================================%
% ======================= ENERGY EFFICIENCY (EE) ======================== %
%=========================================================================%
% Prepare to store energy efficiency results
Ee_MR      = zeros(nbrOfSetups, 1); % MR Precoding
Ee_LP_MMSE = zeros(nbrOfSetups, 1); % LP_MMSE Precoding
Ee_P_MMSE  = zeros(nbrOfSetups, 1); % P_MMSE Precoding
Ee_P_RZF   = zeros(nbrOfSetups, 1); % P_RZF Precoding

%=========================================================================%
% =========== APs PER UE (|M_{k}|), AND UEs PER AP (|D_{l}|) ============ %
%=========================================================================%
% Prepare to store the number of active APs connected to each UE and auxiliar variable for Cmin
ActiveAPs        = zeros(nbrOfSetups,1);
APsConnectedToUE = zeros(nbrOfSetups,nbrOfUEs);