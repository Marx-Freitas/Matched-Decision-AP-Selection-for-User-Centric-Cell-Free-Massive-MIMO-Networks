% This function associates the UE with non-master APs It is a modified
% version of the method proposed in [4]

% INPUT:

% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork

% tau_p             = the number of orthogonal pilot signals
% U_max             = the maximum number of UEs that each AP can serve
% Cmax              = The maximum number of APs per UE
% threshold_dB      = A threshold value in dB utilized by non-master APs
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% mastersOfUEs      = the index of the master AP of each UE
%                     Dim: nbrOfUEs x 1
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs

% Note: Cmax is a modification we have performed in this code. However,
% it is not utilized in these simulations. Thus, we have assumed that
% Cmax = nbrOfAPs
                                
function D = nonMasterSelection(nbrOfAPs, nbrOfUEs, tau_p, U_max, Cmax, ...
    threshold, pilotIndex, mastersOfUEs, gainOverNoise_dB,D)

for indexAP = 1:nbrOfAPs
    tmax = min(tau_p, nbrOfUEs); % Helpful when nbrOfUEs < tau_p
    
    for t = 1:tmax
        
        % Find which UEs are using the pilot "t" in AP L.
        pilotUEs = find(t==pilotIndex(:));        
        
        % If at least one UE is served by the AP in the pilot "t",
        % this condition will not be satisfied
        if (sum(D(indexAP,pilotUEs)) == 0) && (sum(D(indexAP,:)) < U_max)
            
            % Assigns the pilot "t" to the UE presenting the strongest
            % channel gain
            [gainValue,UEindex] = max(gainOverNoise_dB(indexAP,pilotUEs));
            
            % Serve the UE if the channel gain is at most "threshold"
            % weaker than the master AP's channel
            if (gainValue -gainOverNoise_dB(mastersOfUEs(pilotUEs(UEindex)), pilotUEs(UEindex)) >= threshold) && (sum(D(:,pilotUEs(UEindex)))<Cmax)                
                D(indexAP,pilotUEs(UEindex)) = 1;                       
            
            else
                D(indexAP,pilotUEs(UEindex)) = 0;
                
            end

        end

    end

end

% REFERENCES
% [4] E. Björnson and L. Sanguinetti, "Scalable cell-free massive MIMO
% systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247–4261, Jul. 2020.