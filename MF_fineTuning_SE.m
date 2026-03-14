% This function performs the proposed fine-tuning scheme based on spectral
% efficiency, i.e., the algorithm 4
% This function is based in the paper: "Matched-Decision AP Selection for
% User-Centric Cell-Free Massive MIMO Networks"

% Author: Marx Freitas

% INPUT:
% D                 = the AP cluster of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% R                 = the channel covariance matrix (correlation matrix)
%                     Dim: N x N x nbrOfAPs x nbrOfUEs
% noiseVariance_dBm = the noise variance in dBm
% PowerDL_mW        = the total transmission power in each AP
% U_max             = the maximum number of UEs that each AP can serve
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% varepsilon        = threshold for SE losses
% tau_p             = the number of orthogonal pilot signals
% tau_c             = the number of samples per coherence block
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs

function D = MF_fineTuning_SE(D, R, noiseVariance_dBm, PowerDL_mW, U_max, ...
    nbrOfAPs, nbrOfUEs, varepsilon, tau_p, tau_c, masterAPs)

% Prepare to store results
DS_k   = zeros(1,nbrOfUEs); % Desidered signal of UE k
IS_P_k = zeros(nbrOfUEs,nbrOfUEs); % Interference from interfering UEs
DS_kl  = zeros(nbrOfAPs,nbrOfUEs); % Contribution of each serving AP to DS_k
D_aux  = zeros(nbrOfAPs,nbrOfUEs); % Auxiliary AP cluster

% Power allocated to each UE. The equal power allocation is assumed in this
% example for simplicity
varrho_k_W = (PowerDL_mW/U_max(end))*10^-3;
varrho_i_W = varrho_k_W;

%=========================================================================%
% ================= DESIRED SIGNAL AND INTERFERING UEs ================== %
%=========================================================================%

for indexUE = 1:nbrOfUEs

    % Find the APs that are serving the UE k
    servingAPs_id_k = find(D(:,indexUE) == 1);

    % Subset of UEs that are partially served by the same APs as UE k
    P_k = find(sum(D(servingAPs_id_k,:),1)>=1);

    %---------------------------------------------------------------------%
    % DESIRED SIGNAL                                                      
    %---------------------------------------------------------------------%
    % Desired signal of UE k to each AP that is in its AP cluster

    for id_k = 1: length(servingAPs_id_k)

        % Compute the desired signal in (25)
        DS_k(1,indexUE) = DS_k(1,indexUE) + sqrt(varrho_k_W*trace(R(:,:,servingAPs_id_k(id_k), indexUE)));

        % Compute the contribution of each serving AP to the desired signal
        DS_kl(servingAPs_id_k(id_k),indexUE) = sqrt(varrho_k_W*trace(R(:,:,servingAPs_id_k(id_k), indexUE)));
    end

    %---------------------------------------------------------------------%
    % INTERFERENCE FROM INTERFERING UEs
    %---------------------------------------------------------------------%
    % Interfering signals: UEs partially served by the same APs as UE k

    for index_i = 1:length(P_k)

        % Find the APs that serve the UE i
        servingAPs_id_i = find(D(:,P_k(index_i)) == 1);

        for id_i = 1:length(servingAPs_id_i)
            if D(servingAPs_id_i(id_i),indexUE)  == 1
                IS_P_k(indexUE,P_k(index_i)) = IS_P_k(indexUE,P_k(index_i)) + ...
                    varrho_i_W*trace(R(:,:,servingAPs_id_i(id_i),index_i)...
                    *R(:,:,servingAPs_id_i(id_i),indexUE))/trace(R(:,:,servingAPs_id_i(id_i),index_i));
            end
        end
    end
end


%=========================================================================%
% ===================== SPECTRAL EFFICIENCY (SE) ======================== %
%=========================================================================%

% Compute the noise variance
noiseVariance_mW = 10^(noiseVariance_dBm/10);

% Compute the interference signal in (25)
IS_k = transpose(sum(IS_P_k,2)); 
aux_IS_k = IS_k-(abs(DS_k)).^2; % Auxiliary variable

% To fix numerical erros
numErrors_id = aux_IS_k<0;
aux_IS_k(numErrors_id) = IS_k(numErrors_id);

% Removing numerical erros from (25)
IS_k = aux_IS_k;

% Compute the SINR in (13) and the SE in (12)
SINR_k = abs(DS_k).^2./(IS_k+noiseVariance_mW);
SE_k = ((tau_c-tau_p)/tau_c)*log2(1+SINR_k);


%=========================================================================%
% ============================== FINE-TUNING ============================ %
%=========================================================================%
% Fine-tuning spectral efficiency
% Set the vector q_kl
q_kl = DS_kl;

for indexUE = 1:nbrOfUEs
    % Sort each AP (in the AP cluster of UE k)according to its contribution to DS_kl
    [DS_kl_sort, index_q_kl] = sort(q_kl(:,indexUE),'descend');
    % Computes the contribution of each AP to SE_k in percentage
    q_trace_kl_prime_sum = cumsum(DS_kl_sort); % Sums cumulatively the sorted DS_kl
    SINR_kl_prime_sum = abs(q_trace_kl_prime_sum).^2/IS_k(indexUE); % Cumulative and sorted SINR_kl
    SE_kl_prime_sum = ((tau_c-tau_p)/tau_c)*log2(1+SINR_kl_prime_sum); % Cumulative and sorted SE_kl
    index_q_kl_prime_end = (SE_k(indexUE)-SE_kl_prime_sum)<=varepsilon; % Identify the APs with a marginal contribution to the SE_k
    index_q_kl(index_q_kl_prime_end) = []; % Exclude the APs with a marginal contribution to the SE_k
    D_aux(index_q_kl,indexUE) = 1; % Keeps the APs with significant contribution to the SE_k
end

% Final AP clusters
D = or(D_aux, masterAPs);