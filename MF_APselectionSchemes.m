% This code performs the following AP selection strategies:
% Canonical Cell-Free (CF) [1]
% User-centric clustering (UCC) [2]
% Largest-large-scale-fading-based (LSFB) [3]
% Scalable CF [4]
% Matched-Decision (MD), i.e., the proposed framework
% MD LSFB - When the MD scheme behaves as a scalable version of LSFB
% MD UCC  - When the MD scheme behaves as a scalable version of UCC. For
% achieving the MD-UCC scheme, set threshold_dB <= -50

% INPUT:
% APselectionMethod = the AP selection scheme
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% tau_p             = the number of orthogonal pilot signals
% U_max             = the maximum number of UEs that each AP can serve
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs
% pilotIndex = the pilot assigned to each UE
%              Dim: nbrOfUEs x 1
% masterAPs  = the master AP of each UE
%              Dim: nbrOfAPs x nbrOfUEs

% Note: this code uses the terms master AP and coordination AP as synonyms
% Author: Marx M. M. Freitas

function [D, pilotIndex, masterAPs] = MF_APselectionSchemes(APselectionMethod, ...
    nbrOfAPs, nbrOfUEs, tau_p, U_max, gainOverNoise_dB)

% Prepare to store variables
masterAPs            = zeros(nbrOfAPs,nbrOfUEs); % coordination APs
mastersOfUEs         = zeros(nbrOfUEs,1);
pilotIndex           = zeros(nbrOfUEs,1);

%=========================================================================%
% ================= NON-SCALABLE AP SELECTION SCHEMES =================== %
%=========================================================================%
% Canonical Cell-Free (CF) [1]:
if strcmp(APselectionMethod,'CanonicalCF')
    [D, pilotIndex, masterAPs] = MF_canonicalCF(pilotIndex, nbrOfAPs, ...
        nbrOfUEs, gainOverNoise_dB, tau_p, masterAPs);

    % User-centric clustering (UCC) [2]:
elseif strcmp(APselectionMethod,'UCC')
    [D, pilotIndex, masterAPs] = MF_UCC(pilotIndex, nbrOfAPs, nbrOfUEs, ...
        tau_p, U_max, gainOverNoise_dB, masterAPs);

    % Largest-large-scale-fading-based (LSFB) [3]:
elseif strcmp(APselectionMethod,'LSFB')
    [D, pilotIndex, masterAPs] = MF_LSFB(pilotIndex, nbrOfUEs, ...
        nbrOfAPs, gainOverNoise_dB, tau_p);

    %=====================================================================%
    % ================= SCALABLE AP SELECTION SCHEMES =================== %
    %=====================================================================%
    % Scalable CF [4]:
elseif strcmp(APselectionMethod,'ScalableCF')
    [D, pilotIndex, masterAPs] = MF_ScalableCF(pilotIndex, nbrOfUEs, ...
        nbrOfAPs, masterAPs, mastersOfUEs, U_max, tau_p, gainOverNoise_dB);

    % Matched Decision (MD), PROPOSED FRAMEWORK
elseif strcmp(APselectionMethod,'matchedDecision')
    threshold_dB = -50;
    [D, pilotIndex, masterAPs] = MF_matchedDecision(pilotIndex, masterAPs, ...
        mastersOfUEs, nbrOfAPs, nbrOfUEs, U_max, tau_p, gainOverNoise_dB, ...
        threshold_dB);

    % MD LSFB
elseif strcmp(APselectionMethod,'MD_LSFB')
    [D, pilotIndex, masterAPs] = MF_MD_LSFB(pilotIndex, masterAPs, ...
        mastersOfUEs, U_max, nbrOfUEs, nbrOfAPs, gainOverNoise_dB, tau_p);

end

end

% REFERENCES
% [1] H. Q. Ngo, A. Ashikhmin, H. Yang, E. G. Larsson, and T. L. Marzetta,
% "Cell-free massive MIMO versus small cells," IEEE Trans. Wireless
% Commun., vol. 16, no. 3, pp. 1834–1850, Mar. 2017.

% [2] S. Buzzi and C. D'Andrea, "Cell-free massive MIMO: User-centric
% approach," IEEE Wireless Commun. Lett., vol. 6, no. 6, pp. 706–709,
% Dec. 2017.

% [3] H. Q. Ngo, L. Tran, T. Q. Duong, M. Matthaiou, and E. G. Larsson,
% "On the total energy efficiency of cell-free massive MIMO," IEEE Trans.
% Green Commun. Netw, vol. 2, no. 1, pp. 25–39, 2018.

% [4] E. Björnson and L. Sanguinetti, "Scalable cell-free massive MIMO
% systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247–4261, Jul. 2020.