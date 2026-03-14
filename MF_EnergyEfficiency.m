% This function reproduces the energy efficiency (EE) model of the paper
% "On the total energy efficiency of cell-free massive MIMO,"  [3]
% This model is only suitable for distributed implementation, i.e., for the
% MR and LP-MMMSE precoding techniques

% INPUT:
% BW_Hz             = system bandwidth
% D                 = the AP cluster of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% Se_precoding      = the spectral (SE) achieved with the precoding
% technique (MR or LP-MMSE)
% nbrOfAPs          = the number of APs (L) in the newtork
% N                 = the number of antenna elements per AP
% w_precoding       = the precoding matrix
%                     Dim: L x K x nbrOfRealizations
% nbrOfRealizations = the number of realizations

% OUTPUT:
% Ee           = the energy efficiency (EE)


function Ee = MF_EnergyEfficiency (BW_Hz, D, Se_precoding, nbrOfAPs, N, w_precoding, nbrOfRealizations)

%  Input Parameters (from table II)
alpha_m = 0.4;
Ptc_m   = 0.2; % in W
Pbt_m   = 0.25/10^9; % 0.25 W per Gb/s, but the unit is in W/Gb/s
P0_m    = 0.825;

alpha_m = alpha_m*ones(nbrOfAPs,1);
Ptc_m   = Ptc_m.*ones(nbrOfAPs,1);
Pbt_m   = Pbt_m.*ones(nbrOfAPs,1);
Se_m    = zeros(nbrOfAPs,1); % SE of the users served by the m-th AP
aux = zeros(nbrOfAPs,1);

w_precoding       = abs(w_precoding).^2;
mean_w_precoding  = sum(w_precoding,3)/nbrOfRealizations;
power_w_precoding = sum(mean_w_precoding,2);

% Computing Eqs. (18), (39) - fronthaul power consumption
for indexAP = 1:nbrOfAPs
  aux(indexAP,1) = sum(power_w_precoding((indexAP-1)*N+1:indexAP*N));
end

power_w_precoding = aux;
power   = (power_w_precoding./1e3); % Converting from mW to W

P_m = (1./alpha_m).*power + N.*Ptc_m; % Eq. (18)

% If the AP is off it is not consuming energy
APs_off_id = sum(D,2) == 0;
P_m(APs_off_id) = eps;

for indexAP = 1:nbrOfAPs
    ServedUEs_id = D(indexAP,:) == 1;
    Se_m(indexAP) = sum(Se_precoding(ServedUEs_id));
end

% Computing Eq. (19)
P0_m = P0_m.*ones(nbrOfAPs,1); % In Watts
Pbh_m = P0_m + (BW_Hz.*Se_m.*Pbt_m); % Eq.(19)
Ptotal = sum(P_m) + sum(Pbh_m); % Eq.(17)
Se = sum(Se_precoding);
Ee = (BW_Hz.*Se)/Ptotal; % Eq.(21)

% REFERENCE
% [3] H. Q. Ngo, L. Tran, T. Q. Duong, M. Matthaiou, and E. G. Larsson,
% "On the total energy efficiency of cell-free massive MIMO," IEEE Trans.
% Green Commun. Netw, vol. 2, no. 1, pp. 25–39, 2018

% Some information from the reference paper
% Ptotal = sum(P_m) + sum(Pbh_m), for l varying from 1 to L

% ("P_m") is the power consumption at the m-th AP due to the amplifier and the circuit power consumption part
% (including the power consumption of the transceiver chains and the power consumed for signal processing)
% ("Pbh_m") is the power consumed by the backhaul link connecting the CPU and the m-th AP.

% P_m = (1/alpha_m)*power_w_P_MMSE + N*Ptc_m; ou P_m = (1/alpha_m)*E[Xl] + N*Ptc_m;
% ("alpha_m") is such that: 0 < alpha_m ? 1 is the power amplifier efficiency
% ("Ptc_m") is the internal power required to run the circuit components (e.g. converters, mixers, and filters)
% related to each antenna of the m-th AP


% Pbh_m = P0_m + (B*Se*Pbt_m), where:
% ("P0_m") is a fixed power consumption of each backhaul (traffic-independentpower) which may depend on the
% distances between the APs and the CPU and the system topology
% ("Pbt_m") is the traffic-dependent power (in Watt per bit/s)
% ("B") is the system bandwidth
% ("Se") is the sum spectral efficiency of each user "k"
% ("Ee") is the energy efficiency