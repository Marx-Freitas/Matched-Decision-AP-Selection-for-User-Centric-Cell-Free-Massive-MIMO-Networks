% This function performs the Canonical Cell-Free (CF) AP selection scheme
% Author: Marx Freitas

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% tau_p             = the number of orthogonal pilot signals
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs
% pilotIndex = the pilot assigned to each UE
%              Dim: nbrOfUEs x 1
% masterAPs  = the master AP of each UE
%              Dim: nbrOfAPs x nbrOfUEs

function [D, pilotIndex, masterAPs] = MF_canonicalCF(pilotIndex, nbrOfAPs, nbrOfUEs, gainOverNoise_dB, tau_p, masterAPs)

%=========================================================================%
% ====================== CANONICAL CF AP SELECTION ====================== %
%=========================================================================%
D = ones(nbrOfAPs,nbrOfUEs);
gainOverNoise_Lin = db2pow(gainOverNoise_dB);

%=========================================================================%
% ========================== PILOT ASSIGNMENT =========================== %
%=========================================================================%
for indexUE = 1:nbrOfUEs

    servingAPs = D(:, indexUE) == 1;
    [~, strongestAP_id]  = max(gainOverNoise_Lin(servingAPs, indexUE));

    % Assign orthogonal pilots to the first tau_p UEs
    pilotIndex = PilotAssignment(pilotIndex,indexUE,tau_p,gainOverNoise_dB,strongestAP_id);

    % Master AP assignment
    masterAPs(strongestAP_id, indexUE) = 1;

end

% REFERENCES
% [1] H. Q. Ngo, A. Ashikhmin, H. Yang, E. G. Larsson, and T. L. Marzetta,
% "Cell-free massive MIMO versus small cells," IEEE Trans. Wireless
% Commun., vol. 16, no. 3, pp. 1834–1850, Mar. 2017.