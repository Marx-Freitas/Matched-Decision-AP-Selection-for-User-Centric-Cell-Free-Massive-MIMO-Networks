% This function performs the Largest-large-scale-fading-based (LSFB) AP
% selection scheme
% Author: Marx Freitas

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% nbrOfUEs          = the number of UEs (K) in the newtork
% nbrOfAPs          = the number of APs (L) in the newtork
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

function [D, pilotIndex, masterAPs] = MF_LSFB(pilotIndex, nbrOfUEs, nbrOfAPs, gainOverNoise_dB, tau_p, masterAPs)


% Set the threshold gain of LSFB scheme
sigma_percent = 0.999;

% Prepare to store auxiliar variables
APtoUE_normGain   = zeros(nbrOfAPs,nbrOfUEs);
normGain_sorted   = zeros(nbrOfAPs,nbrOfUEs); % Normalized sorted gain
normGain_cumula   = zeros(nbrOfAPs,nbrOfUEs); % Normalized sorted gain,
                                              % cumulative sum

% Prepare to store auxiliar variables
indexAP_sorted    = zeros(nbrOfAPs,nbrOfUEs);
selectedAP        = zeros(nbrOfAPs,nbrOfUEs);
indexAP_flip      = zeros(nbrOfAPs,nbrOfUEs);

% Convert the channel gain from dB to linear
gainOverNoise_lin = db2pow(gainOverNoise_dB);

% Prepare to store the AP cluster
D = zeros(nbrOfAPs,nbrOfUEs);

%=========================================================================%
% ========================== LSFB AP selection ========================== %
%=========================================================================%
for indexUE = 1:nbrOfUEs
    
    % Compute the channel gain of each AP normalized by the sum of the UE channel gain
    APtoUE_normGain(:,indexUE) = gainOverNoise_lin(:,indexUE).*(1./sum(gainOverNoise_lin(:,indexUE))+eps);
    [normGain_sorted(:,indexUE), indexAP_sorted(:,indexUE)] = sort(APtoUE_normGain(:,indexUE));
    
    % Sort the normalized channel gains in ascending order
    indexAP_flip(:,indexUE)    = flip(indexAP_sorted(:,indexUE));
    normGain_cumula(:,indexUE) = cumsum(flip(normGain_sorted(:,indexUE)));

    % Select the APs with the largest contribution to the total channel gain
    selectedAP(:,indexUE) = normGain_cumula(:,indexUE) < sigma_percent;
    lastSelectedAP_id = sum(selectedAP,1);
    
    % Set the AP cluster
    D(indexAP_flip(1:lastSelectedAP_id(1,indexUE),indexUE),indexUE) = 1;

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

% REFERENCES
% [3] H. Q. Ngo, L. Tran, T. Q. Duong, M. Matthaiou, and E. G. Larsson,
% "On the total energy efficiency of cell-free massive MIMO," IEEE Trans.
% Green Commun. Netw, vol. 2, no. 1, pp. 25–39, 2018.