% This function performs the user-centric clustering (UCC) [2] AP selection
% scheme. It is a modified version of [2], since [2] considers the
% estimated channel in each coherence block to perform AP selection. In our
% case, we have considered that [2] utilizes only the large-scale fading
% coefficients.

% Author: Marx Freitas

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% tau_p             = the number of orthogonal pilot signals
% U_max             = the maximum number of UEs that each AP can serve
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs
% pilotIndex = the pilot assigned to each UE
%              Dim: nbrOfUEs x 1
% masterAPs  = the master AP of each UE
%              Dim: nbrOfAPs x nbrOfUEs


function [D, pilotIndex, masterAPs] = MF_UCC(pilotIndex, nbrOfAPs, nbrOfUEs, tau_p, U_max, gainOverNoise_dB, masterAPs)
% Prepare to store the AP cluster of the UE
D = zeros(nbrOfAPs,nbrOfUEs);

% Compute the maximum number of UEs that each AP can serve
nbrUEstoSelect = min(U_max,nbrOfUEs);

% Convert the channel gain from dB to linear
gainOverNoise_lin = db2pow(gainOverNoise_dB);

%=========================================================================%
% ========================== UCC AP selection =========================== %
%=========================================================================%
for indexAP = 1:nbrOfAPs
    [~,indexUE] = sort(gainOverNoise_lin(indexAP,:));
    indexUE = indexUE(end-nbrUEstoSelect+1:end);
    D(indexAP,indexUE) = 1;
end

%=========================================================================%
% ========================== PILOT ASSIGNMENT =========================== %
%=========================================================================%
gainOverNoise_lin(D == 0) = 0; % It consider only the APs serving the UE

for indexUE = 1:nbrOfUEs
    
    % Identify the AP serving the UE presenting the strongest channel gain 
    [~,strongestAP_id] = max(gainOverNoise_lin(:,indexUE));
    
    % Assign orthogonal pilots to the first tau_p UEs
    pilotIndex = PilotAssignment(pilotIndex,indexUE,tau_p,gainOverNoise_dB,strongestAP_id);

    % The AP with strongest channel gain is considered the "master" AP of
    % the UE. Note that there is no master AP in this scheme. This line is
    % used only to allow the code to run in different AP selection schemes
    masterAPs(strongestAP_id, indexUE) = 1;

end

% References
% [2] S. Buzzi and C. D'Andrea, "Cell-free massive MIMO: User-centric
% approach," IEEE Wireless Commun. Lett., vol. 6, no. 6, pp. 706–709,
% Dec. 2017.