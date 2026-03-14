% This function performs pilot assignment. It follows the approach proposed
% in [4]

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% indexUE           = the index of UE k 
% tau_p             = number of ulplink pilot samples per coherence block
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% indexSortedAPs    = index sorted APs

% OUTPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1


%Pilot assignment
function pilotIndex = PilotAssignment(pilotIndex,indexUE,tau_p,gainOverNoise_dB,indexSortedAPs)
gainOverNoise_Lin = db2pow(gainOverNoise_dB);

if indexUE <= tau_p
    pilotIndex(indexUE) = indexUE;
else %Assign pilot for remaining UEs
    %Compute received power from to the master AP from each pilot
    pilotinterference = zeros(tau_p,1);    
    for t = 1:tau_p
        pilotinterference(t) = sum(gainOverNoise_Lin(indexSortedAPs(end), pilotIndex(1:indexUE-1)==t));
    end
    %Find the pilot with the least receiver power
    [~,bestpilot] = min(pilotinterference);
    pilotIndex(indexUE) = bestpilot;    
end


% REFERENCES
% [4] E. Björnson and L. Sanguinetti, "Scalable cell-free massive MIMO
% systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247–4261, Jul. 2020.