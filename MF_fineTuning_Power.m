% This function performs the proposed fine-tuning scheme based on power
% allocation, i.e., the algorithm 3
% This function is based in the paper: "Matched-Decision AP Selection for
% User-Centric Cell-Free Massive MIMO Networks"

% Author: Marx Freitas

% INPUT:
% D                 = the AP cluster of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% PowerDL_mW        = the total transmission power in each AP
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% Gamma             = the Γ% of the cumulative sum in (24) 

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs

function D = MF_fineTuning_Power(D, PowerDL_mW, nbrOfAPs, nbrOfUEs, gainOverNoise_dB, masterAPs, Gamma_kl)

% Auxiliary AP cluster
D_aux = D;

% Prepare to store the power allocation and final AP clusters
rho_dist_mW       = zeros(nbrOfAPs,nbrOfUEs);
FinalAPCluster    = zeros(nbrOfAPs,nbrOfUEs);

% Converting the channel gain from dB to linear
gainOverNoise_Lin = db2pow(gainOverNoise_dB);

% Fine-tuning the AP cluster
for indexAP =1:nbrOfAPs

    % Find the UEs served by the AP
    servedUEs_id = find(D_aux(indexAP,:) == 1);
    normalizationAPl = sum(sqrt(gainOverNoise_Lin(indexAP,servedUEs_id)));
    
    % Compute the power that the AP allocates to each served UE
    rho_dist_mW(indexAP, servedUEs_id) = PowerDL_mW*sqrt(gainOverNoise_Lin(indexAP,servedUEs_id))/normalizationAPl;
    
    % Sort the UEs served in descending order according to allocated power
    [rho_prime_kl_mW, id_k] = sort(rho_dist_mW(indexAP,:)/PowerDL_mW,'descend');
    
    % The AP serves only the UEs that contribute to at least Γ% of the
    % cumulative sum in (24)
    id_aux = sum(cumsum(rho_prime_kl_mW) <= Gamma_kl);
    id_k = id_k(1:id_aux);
    FinalAPCluster(indexAP, id_k) = 1;
end

% Final AP cluster
D = or(FinalAPCluster,masterAPs);