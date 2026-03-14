% This function associates UEs with master APs. It is a modified version of
% the method proposed in [4]. The main difference between this function and
% [4] is that this function allows UEs to connect to multiple master APs.
% In contrast, [4] allows the UE to only connect to a single master AP.
% Furthermore, [4] does not consider that an AP may be unable to serve a
% new UE if the AP is operating at its maximum capacity to serve UEs and if
% all these UEs use the AP as master. However, our simulations usually also
% consider a single master AP per UE

% This function performs the  AP selection scheme
% Author: Marx Freitas

% INPUT:
% masterAPs        = the master AP of each UE
%                    Dim: nbrOfAPs x nbrOfUEs
% mastersOfUEs     = the index of the master AP of each UE
%                    Dim: nbrOfUEs x 1
% Cmin_nbrAPs      = the minimum number of master APs per UE
% U_max            = the maximum number of UEs that each AP can serve
% indexUE          = the index of the UE
% gainOverNoise_dB = the channel gain normalized by noise
%                    Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% masterAPs    = the master AP of each UE
%                Dim: nbrOfAPs x nbrOfUEs
% mastersOfUEs = the index of the master AP of each UE
%                Dim: nbrOfUEs x 1


function [masterAPs, mastersOfUEs] = MF_mastersAPs(masterAPs, mastersOfUEs, Cmin_nbrAPs, U_max, indexUE, gainOverNoise_dB)

% Convert gainOverNoise_dB from dB to linear
gainOverNoise_Lin = db2pow(gainOverNoise_dB);

% Auxiliar variable
auxGainOverNoiseLin = gainOverNoise_Lin;

% Set auxGainOverNoiseLin to zero for all APs serving U_max UEs as master
auxGainOverNoiseLin((sum(masterAPs(:,:),2)>=U_max),:,:) = 0;

% Sort the channel gains of the APs in descending order
[gainSortedAPs,indexSortedAPs] = sort(auxGainOverNoiseLin(:, indexUE));
gainSortedAPs  = gainSortedAPs(end-Cmin_nbrAPs(indexUE,1)+1:end); 
indexSortedAPs = indexSortedAPs(end-Cmin_nbrAPs(indexUE,1)+1:end);

% Excludes the APs presenting auxGainOverNoiseLin = 0
indexSortedAPs(gainSortedAPs(:,1) == 0) = []; 

% Set the master APs of the UE
masterAPs(indexSortedAPs,indexUE) = 1;

% Choose the UE's main master AP
mastersOfUEs(indexUE) = indexSortedAPs(end);

end

% REFERENCES
% [4] E. Björnson and L. Sanguinetti, "Scalable cell-free massive MIMO
% systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247–4261, Jul. 2020.