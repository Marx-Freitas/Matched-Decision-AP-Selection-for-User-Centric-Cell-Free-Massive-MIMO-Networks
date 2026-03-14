% This function performs the proposed fine-tuning scheme based on energy
% efficiency, i.e., the algorithm 5
% This function is based in the paper: "Matched-Decision AP Selection for
% User-Centric Cell-Free Massive MIMO Networks"

% INPUT:
% D                 = the AP cluster of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% R                 = the channel covariance matrix (correlation matrix)
%                     Dim: N x N x nbrOfAPs x nbrOfUEs
% noiseVariance_dBm = the noise variance
% PowerDL_mW        = the total transmission power in each AP
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% tau_p             = the number of orthogonal pilot signals
% tau_c             = the number of samples per coherence block
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% gainOverNoisedB   = the channel gain normalized by noise variance of UE k
%                     regarding all APs normalized by the noise variance.
%                     Dim: nbrOfAPs x 1
% zeta              = threshold for algorithm 5

% OUTPUT:
% D                 = the AP cluster of each UE
%                     Dim: nbrOfAPs x nbrOfUEs

function D = MF_fineTuning_EE(D ,R, noiseVariance_dBm, PowerDL_mW, ...
    nbrOfAPs, nbrOfUEs, tau_p, tau_c, masterAPs, gainOverNoise_dB, zeta)

% Prepare to store results
DS_k   = zeros(1,nbrOfUEs); % Desidered signal of UE k
DS_kl  = zeros(nbrOfAPs,nbrOfUEs); % Interference from interfering UEs
IS_P_k = zeros(nbrOfUEs,nbrOfUEs); % Contribution of each serving AP to DS_k
D_aux  = zeros(nbrOfAPs,nbrOfUEs); % Auxiliary AP cluster

% Power allocated to each UE. The equal power allocation is assumed in this
% example for simplicity
varrho_kl = DLpowerAllocation(nbrOfUEs, nbrOfAPs, D, PowerDL_mW,...
    gainOverNoise_dB);

%=========================================================================%
% ================= DESIRED SIGNAL AND INTERFERING UEs ================== %
%=========================================================================%
% Computing the desired signal and interference

for indexUE = 1:nbrOfUEs

    % Find the APs that serve the UE k
    servingAPs_id_k = find(D(:,indexUE) == 1);

    % Subset of UEs that are partially served by the same APs as UE k
    P_k = find(sum(D(servingAPs_id_k,:),1)>=1);

    %---------------------------------------------------------------------%
    % DESIRED SIGNAL
    %---------------------------------------------------------------------%
    % Desired signal of UE k to each AP that is in its AP cluster

    for id_k = 1: length(servingAPs_id_k)
        % Compute the desired signal in (25)
        DS_k(1,indexUE) = DS_k(1,indexUE) + sqrt(varrho_kl(servingAPs_id_k(id_k),indexUE)*trace(R(:,:,servingAPs_id_k(id_k),indexUE)));

        % Compute the contribution of each serving AP to the desired signal
        DS_kl(servingAPs_id_k(id_k),indexUE) = sqrt(varrho_kl(servingAPs_id_k(id_k),indexUE)*trace(R(:,:,servingAPs_id_k(id_k),indexUE)));
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
                    varrho_kl(servingAPs_id_k(id_k),index_i).*trace(R(:,:,servingAPs_id_i(id_i),index_i)...
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
% Fine-tuning energy efficiency

% Identify the SE of each UE served by the AP
SE_k_l = repmat(SE_k, [nbrOfAPs 1]);
SE_k_l = SE_k_l.*D;

% Compute the ratio between the SE and power allocated
rk_l = SE_k_l./(varrho_kl+eps);

for indexAP = 1:nbrOfAPs

    % Find the UEs that the AP serves
    servedUEs_id   = find(D(indexAP,:) == 1);
    masteredUEs = find(masterAPs(indexAP, :) == 1);
    servedUEs_id   = setdiff(servedUEs_id, masteredUEs);

    if isempty (servedUEs_id)
        % Do nothing
    else
        % Keeping and dropping UEs
        rk_l_max = max(rk_l(indexAP, servedUEs_id));
        tolerance = zeta/rk_l_max;
        rk_l_minus_1 = 1./rk_l(indexAP, servedUEs_id);
        RemainingUEs_id = rk_l_minus_1 > tolerance;

        % Keep or drop some UEs
        D_aux(indexAP, servedUEs_id) = 0;
        D_aux(indexAP, servedUEs_id(RemainingUEs_id)) = 1;

    end
end

% Final AP clusters
D_aux = or(D_aux, masterAPs);
D = D_aux;