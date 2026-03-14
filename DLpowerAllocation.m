% Author: Marx Freitas
% This function performs the downlink (DL) power allocation for distributed
% implementations. Thus, it is used for MR and LP-MMSE precoding schemes.

% INPUT:
% nbrOfUEs        = the number of user equipment (UEs) in the network
% nbrOfAPs        = the number of access points (APs) in the network
% D               = the AP cluster of each UE
%                   Dim: nbrOfAPs x nbrOfUEs
% PowerDL_mW      = the total transmission power on each AP in mW
% gainOverNoisedB = the imperfect channel gain normalized by noise
%                   variance of UE k regarding all APs normalized by the
%                   noise variance.

% OUTPUT:
% rho_dist_mW     = the power allocated to each UE
%                   Dim: nbrOfAPs x nbrOfUEs

function rho_dist_mW = DLpowerAllocation(nbrOfUEs, nbrOfAPs, D, PowerDL_mW, gainOverNoise_dB)

% Prepare to store the power allocation for distributed precoding
rho_dist_mW = zeros(nbrOfAPs,nbrOfUEs);

% Compute the channel gain of the masters and non-masters APs
gainOverNoise_lin = db2pow(gainOverNoise_dB);

for indexAP = 1:nbrOfAPs

    % Identifies which UEs are served by AP l
    servedUEs = find(D(indexAP,:)==1);

    for indexServedUEs = 1:length(servedUEs)

        % Compute the power allocation for distributed precoding
        rho_dist_mW(indexAP,servedUEs(indexServedUEs)) = PowerDL_mW*sqrt(gainOverNoise_lin(indexAP,servedUEs(indexServedUEs)))/sum(sqrt(gainOverNoise_lin(indexAP,servedUEs)));
    end
end

end

