% This function performs the Scalable CF AP selection scheme
% Author: Marx Freitas

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% nbrOfUEs          = the number of UEs (K) in the newtork
% nbrOfAPs          = the number of APs (L) in the newtork
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% mastersOfUEs      = the index of the master AP of each UE
%                     Dim: nbrOfUEs x 1
% U_max             = the maximum number of UEs that each AP can serve
% tau_p             = the number of orthogonal pilot signals
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs
% pilotIndex = the pilot assigned to each UE
%              Dim: nbrOfUEs x 1
% masterAPs  = the master AP of each UE
%              Dim: nbrOfAPs x nbrOfUEs

function [D, pilotIndex, masterAPs] = MF_ScalableCF(pilotIndex, nbrOfUEs, ...
    nbrOfAPs, masterAPs, mastersOfUEs, U_max, tau_p, gainOverNoise_dB)

% Prepare to store the AP cluster
D = zeros(nbrOfAPs,nbrOfUEs);

% Set the threshold value in dB of non-master AP selection
threshold_dB = -40;

% Set one master AP per UE
nbrOfmasterAPs = ones(nbrOfUEs,1);

for indexUE = 1:nbrOfUEs
    %=====================================================================%
    % ========================== MASTER AP ============================== %
    %=====================================================================%
    % Associate the UE with a master AP (or masters APs)
    [masterAPs, mastersOfUEs] = MF_mastersAPs(masterAPs, mastersOfUEs,...
        nbrOfmasterAPs, U_max, indexUE, gainOverNoise_dB);
    D = masterAPs;

    %=====================================================================%
    % ======================= PILOT ASSIGNMENT ========================== %
    %=====================================================================%
    pilotIndex = PilotAssignment(pilotIndex, indexUE, tau_p, ...
        gainOverNoise_dB, mastersOfUEs(indexUE));

end

%=========================================================================%
% ========================= NON MASTER APs ============================== %
%=========================================================================%
D = nonMasterSelection(nbrOfAPs, nbrOfUEs, tau_p, U_max, nbrOfAPs, ...
    threshold_dB, pilotIndex, mastersOfUEs, gainOverNoise_dB, D);

% REFERENCES
% [4] E. Björnson and L. Sanguinetti, "Scalable cell-free massive MIMO
% systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247–4261, Jul. 2020.