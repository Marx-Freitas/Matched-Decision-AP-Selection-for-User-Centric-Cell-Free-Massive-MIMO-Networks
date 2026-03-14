% This function performs the MD LSFB AP selection scheme
% Author: Marx Freitas

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% mastersOfUEs      = the index of the master AP of each UE
%                     Dim: nbrOfUEs x 1
% U_max             = the maximum number of UEs that each AP can serve
% nbrOfUEs          = the number of UEs (K) in the newtork
% nbrOfAPs          = the number of APs (L) in the newtork
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% tau_p             = the number of orthogonal pilot signals


% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs
% pilotIndex = the pilot assigned to each UE
%              Dim: nbrOfUEs x 1
% masterAPs  = the master AP of each UE
%              Dim: nbrOfAPs x nbrOfUEs

function [D, pilotIndex, masterAPs] = MF_MD_LSFB(pilotIndex, masterAPs, ...
    mastersOfUEs, U_max, nbrOfUEs, nbrOfAPs, gainOverNoise_dB, tau_p)

% A threshold value in dB utilized by MD LSFB
sigma_percent = 0.999;

% Prepare to store the decisions taken by the UE and APs
E = zeros(nbrOfAPs,nbrOfUEs); % Decision taken by the UE
F = zeros(nbrOfAPs,nbrOfUEs); % Decision taken by the AP

% Prepare to store the AP cluster
D = zeros(nbrOfAPs,nbrOfUEs); % AP clusters

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

for indexUE = 1:nbrOfUEs

    %=====================================================================%
    % =================== DECISION TAKEN BY THE UE ====================== %
    %=====================================================================%
    
    %---------------------------------------------------------------------%
    % THE DECISIONS TAKEN BY THE UE ARE A PART OF INTERMEDIATE STAGE
    %---------------------------------------------------------------------%

    % The UE requests connection to a subset of APs following decision
    % criterion. In this case, the UE selects the APs with the largest
    % contribution to the total channel gain, i.e., the LSFB scheme
    
    % Compute the channel gain of each AP normalized by the sum of the UE channel gain
    APtoUE_normGain(:,indexUE) = gainOverNoise_lin(:,indexUE).*(1./sum(gainOverNoise_lin(:,indexUE))+eps);
    [normGain_sorted(:,indexUE),indexAP_sorted(:,indexUE)] = sort(APtoUE_normGain(:,indexUE));
    
    % Sort the normalized channel gains in ascending order
    indexAP_flip(:,indexUE) = flip(indexAP_sorted(:,indexUE));
    normGain_cumula(:,indexUE) = cumsum(flip(normGain_sorted(:,indexUE)));

    % Select the APs with the largest contribution to the total channel gain
    selectedAP(:,indexUE) = normGain_cumula(:,indexUE)<sigma_percent;
    lastSelectedAP_id = sum(selectedAP,1);
    
    % Decision taken by the UE
    E(indexAP_flip(1:lastSelectedAP_id(1,indexUE),indexUE),indexUE) = 1;
    E_id = find(E(:,indexUE) == 1);

    %=====================================================================%
    % ========================== MASTER AP ============================== %
    %=====================================================================%
    % Set the number of master APs per UE. It is considered one per UE
    nbrOfmasterAPs = ones(nbrOfUEs,1);
    [masterAPs, mastersOfUEs] = MF_mastersAPs(masterAPs, mastersOfUEs, ...
        nbrOfmasterAPs, U_max, indexUE, gainOverNoise_dB);

    % Ensure connection to at least one AP, i.e., the master AP
    if sum(F(mastersOfUEs(indexUE),:),2) >= U_max(end)
        servedUEs_id  = find(F(mastersOfUEs(indexUE),:) == 1);
        weakerUE_gain = min(gainOverNoise_lin(mastersOfUEs(indexUE),servedUEs_id));
        weakerUE_pos  = gainOverNoise_lin(mastersOfUEs(indexUE),servedUEs_id) == weakerUE_gain;
        weakerUE_id   = servedUEs_id(weakerUE_pos);
        F(mastersOfUEs(indexUE),weakerUE_id) = 0;
    end

    F(mastersOfUEs(indexUE),indexUE) = 1; 
    E(mastersOfUEs(indexUE),indexUE) = 1;

    %=====================================================================%
    % ======================== PILOT ASSIGNMENT ========================= %
    %=====================================================================%
    [~,indexSortedAPs] = max(gainOverNoise_lin(:,indexUE));
    pilotIndex = PilotAssignment(pilotIndex, indexUE, tau_p, ...
        gainOverNoise_dB, indexSortedAPs);

    %=====================================================================%
    % =================== DECISION TAKEN BY THE APs ===================== %
    %=====================================================================%

    for indexAP = 1:length(E_id)

        % Compute the number of UEs the AP serves as master
        UEsdiscard_id  = masterAPs(E_id(indexAP),:) == 1;

        % Compute the served UEs
        servedUEs_id   = find(F(E_id(indexAP),:) == 1); % ID
        nbrOfServedUEs = sum(F(E_id(indexAP),:),2); % Number of served UEs
        
        % Do not consider the UEs that use the AP as master
        gainOverNoise_lin(E_id(indexAP),UEsdiscard_id) = 10^100;

    %---------------------------------------------------------------------%
    % THESE LINES INCLUDE BOTH THE INTERMEDIATE AND FINAL STAGE IN THE APs
    %---------------------------------------------------------------------%

        if nbrOfServedUEs < U_max(end)
            F(E_id(indexAP),indexUE) = 1;
        else
            weakerUE_gain = min(gainOverNoise_lin(E_id(indexAP),servedUEs_id));
            weakerUE_pos  = gainOverNoise_lin(E_id(indexAP),servedUEs_id) == weakerUE_gain;
            weakerUE_id   = servedUEs_id(weakerUE_pos);
            newUE_gain    = gainOverNoise_lin(E_id(indexAP),indexUE);

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
