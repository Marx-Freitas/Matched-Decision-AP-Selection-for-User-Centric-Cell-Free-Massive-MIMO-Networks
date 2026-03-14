% This function performs the Matched Decision (MD) AP selection scheme
% Author: Marx Freitas

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% mastersOfUEs      = the index of the master AP of each UE
%                     Dim: nbrOfUEs x 1
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% U_max             = the maximum number of UEs that each AP can serve
% tau_p             = the number of orthogonal pilot signals
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% threshold_dB      = A threshold value in dB utilized by MD scheme

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs
% pilotIndex = the pilot assigned to each UE
%              Dim: nbrOfUEs x 1
% masterAPs  = the master AP of each UE
%              Dim: nbrOfAPs x nbrOfUEs


function [D, pilotIndex, masterAPs] = MF_matchedDecision(pilotIndex, ...
    masterAPs, mastersOfUEs, nbrOfAPs, nbrOfUEs, U_max, tau_p, ...
    gainOverNoise_dB, threshold_dB)

% Prepare to store data
E = zeros(nbrOfAPs,nbrOfUEs); % Decision taken by the UE
F = zeros(nbrOfAPs,nbrOfUEs); % Decision taken by the AP

% Prepare to store the decisions taken by the UE and APs
D = zeros(nbrOfAPs,nbrOfUEs); % AP clusters

gainOverNoise_lin    = db2pow(gainOverNoise_dB); % From dB to linear
auxGainOverNoise_lin = db2pow(gainOverNoise_dB); % Auxiliar variable


for indexUE = 1:nbrOfUEs
    %=====================================================================%
    % =================== DECISION TAKEN BY THE UE ====================== %
    %=====================================================================%
    
    %---------------------------------------------------------------------%
    % THE DECISIONS TAKEN BY THE UEs
    %---------------------------------------------------------------------%
    % IT IS A PART OF THE INTERMEDIATE STAGE

    % The UE requests connection to a subset of APs following decision
    % criterion. In this case, it is considered beta > threshold_dB

    % Sort the the large-scale fading gains in ascending order
    [~, APsSelectedByTheUE_id] = sort(gainOverNoise_lin(:,indexUE));

    % Organize the large-scale fading gains in descending order
    APsSelectedByTheUE_id = flip(APsSelectedByTheUE_id);
    APsSelectedByTheUE_id = APsSelectedByTheUE_id(1:end);

    % Compute vector E (decision taken by the UE)
    E(APsSelectedByTheUE_id,indexUE) = 1;

    % Drop the APs with the weakest channel gains
    indexAPsRejected = gainOverNoise_dB(:,indexUE) < threshold_dB;
    E(indexAPsRejected,indexUE) = 0;
    E_id = find(E(:,indexUE) == 1);

    %=====================================================================%
    % ========================== MASTER AP ============================== %
    %=====================================================================%
    % Set the number of master APs per UE. It is considered one per UE
    nbrOfmasterAPs = ones(nbrOfUEs,1);
    [masterAPs,mastersOfUEs] = MF_mastersAPs(masterAPs, mastersOfUEs, ...
        nbrOfmasterAPs, U_max, indexUE, gainOverNoise_dB);

    % Ensure connection to at least one AP, i.e., the master AP
    if sum(F(mastersOfUEs(indexUE),:),2) >= U_max(end)
        servedUEs_id  = find(F(mastersOfUEs(indexUE),:) == 1);
        weakerUE_gain = min(auxGainOverNoise_lin(mastersOfUEs(indexUE),servedUEs_id));
        weakerUE_pos  = auxGainOverNoise_lin(mastersOfUEs(indexUE),servedUEs_id) == weakerUE_gain;
        weakerUE_id   = servedUEs_id(weakerUE_pos);
        F(mastersOfUEs(indexUE),weakerUE_id) = 0;
    end

    F(mastersOfUEs(indexUE),indexUE) = 1;
    E(mastersOfUEs(indexUE),indexUE) = 1;

    %=====================================================================%
    % ======================= PILOT ASSIGNMENT ========================== %
    %=====================================================================%
    pilotIndex = PilotAssignment(pilotIndex, indexUE, tau_p, ...
        gainOverNoise_dB, mastersOfUEs(indexUE));

    %=====================================================================%
    % =================== DECISION TAKEN BY THE APs ===================== %
    %=====================================================================%

    for indexAP = 1:length(E_id)
        
        % Compute the number of UEs the AP serves as master
        UEsdiscard_id  = masterAPs(E_id(indexAP),:) == 1;

        % Compute the served UEs
        servedUEs_id   = find(F(E_id(indexAP),:) == 1); % ID
        nbrOfservedUEs = sum(F(E_id(indexAP),:),2); % Number of served UEs

        % Do not consider the UEs that use the AP as master
        auxGainOverNoise_lin(E_id(indexAP),UEsdiscard_id) = 10^100;
    
    %---------------------------------------------------------------------%
    % THESE LINES INCLUDE BOTH THE INTERMEDIATE AND FINAL STAGE IN THE APs
    %---------------------------------------------------------------------%
        if nbrOfservedUEs < U_max(end)
            F(E_id(indexAP),indexUE) = 1;
        else
            weakerUE_gain = min(auxGainOverNoise_lin(E_id(indexAP),servedUEs_id));
            weakerUE_pos  = auxGainOverNoise_lin(E_id(indexAP),servedUEs_id) == weakerUE_gain;
            weakerUE_id   = servedUEs_id(weakerUE_pos);
            newUE_gain    = auxGainOverNoise_lin(E_id(indexAP),indexUE);

            if newUE_gain > weakerUE_gain
                F(E_id(indexAP),weakerUE_id) = 0;
                F(E_id(indexAP),indexUE) = 1;
            end
        end
    end
end

%=========================================================================%
% ========================== MATCHED DECISION =========================== %
%=========================================================================%
G = and(E,F); G = or(masterAPs,G); D(G) = 1;

end