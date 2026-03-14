function [SE_P_MMSE, SE_P_RZF, SE_LP_MMSE, SE_MR, w_MR, w_LP_MMSE, w_P_MMSE, w_PRZF]...
    = functionComputeSE_downlinkWithRician(Hhat, H, D, C, tau_c, tau_p, nbrOfRealizations, N, K, L, p,...
    rho_dist_mW, gainOverNoisedB, PowerDL_mW)

% Compute downlink SE for different transmit precoding schemes using the
% capacity bound in Theorem 6.1 for the centralized schemes and the
% capacity bound in Corollary 6.3 for the distributed schemes. Compute the
% genie-aided downlink SE from Corollary 6.6 for the centralized and the
% distributed operations.

% Although the function name contains the term Rician, the small-scale
% fading in this code can be simulated as both Rician and Rayleigh fading


% INPUT:
% Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
% H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                     with the true channel realizations. The matrix is
%                     organized in the same way as Hhat.
% D                 = DCC matrix for cell-free setup with dimension L x K
%                     where (l,k) is one if AP l serves UE k and zero otherwise
% C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                     the spatial correlation matrix of the channel
%                     estimation error between AP l and UE k,
%                    normalized by noise variance
% tau_c             = Length of coherence block
% tau_p             = Length of pilot sequences
% nbrOfRealizations = Number of channel realizations
% N                 = Number of antennas per AP
% K                 = Number of UEs
% L                 = Number of APs
% p                 = Uplink transmit power per UE (same for everyone)
% rho_dist          = Matrix with dimension L x K where (l,k) is the power
%                     allocated to UE k by AP l in the distributed downlink
%                     operation
% gainOverNoisedB   = Matrix with dimension L x K where (l,k) is the channel
%                     gain (normalized by the noise variance) between AP l
%                     and UE k
% PowerDL_mW           = Maximum allowed transmit power for each AP


% OUTPUT:
% SE_P_MMSE         = SEs achieved with P-MMSE precoding in (6.17)
% SE_P_RZF          = SEs achieved with P-RZF precoding in (6.18)
% SE_LP_MMSE        = SEs achieved with LP-MMSE precoding in (6.33)
% SE_MR             = SEs achieved with MR precoding in (6.26)

% w_MR              = The precoding matrix of MR scheme.
%                     Dim: L x K x nbrOfRealizations
% w_LP_MMSE         = The precoding matrix of LP-MMSE scheme.
%                     Dim: L x K x nbrOfRealizations
% w_P_MMSE          = The precoding matrix of P-MMSE scheme.
%                     Dim: L x K x nbrOfRealizations
% w_PRZF            = The precoding matrix of P-RZF scheme.
%                     Dim: L x K x nbrOfRealizations


% This Matlab function was developed to generate simulation results to:
%
% Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
% "Foundations of User-Centric Cell-Free Massive MIMO",
% Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
% pp 162-472. DOI: 10.1561/2000000109
%
% This is version 1.0 (Last edited: 2021-01-31)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.

% Some comments: Marx Freitas

% Store the N x N identity matrix
eyeN = eye(N);

% Compute the prelog factor assuming only downlink data transmission
prelogFactor = (1-tau_p/tau_c);

% Prepare to store the terms that appear in SEs
scaling_MR = zeros(L,K);
interUserGains_MR = zeros(K,K,nbrOfRealizations);
signal_MR1 = zeros(K,1);
interf_MR1 = zeros(K,1);

signal_P_MMSE = zeros(K,1);
interf_P_MMSE = zeros(K,1);
scaling_P_MMSE = zeros(K,1);
portionScaling_PMMSE = zeros(L,K);
interUserGains_P_MMSE = zeros(K,K,nbrOfRealizations);

signal_P_RZF = zeros(K,1);
interf_P_RZF = zeros(K,1);
scaling_P_RZF = zeros(K,1);
portionScaling_PRZF = zeros(L,K);
interUserGains_P_RZF = zeros(K,K,nbrOfRealizations);

signal_LP_MMSE = zeros(K,1);
interf_LP_MMSE = zeros(K,1);
scaling_LP_MMSE = zeros(L,K);
interUserGains_LP_MMSE = zeros(K,K,nbrOfRealizations);

% MF editing
% Prepare to store the power consumed in the APs in each precoding scheme
w_MR      = zeros(L*N,K,nbrOfRealizations);
w_LP_MMSE = zeros(L*N,K,nbrOfRealizations);
w_P_MMSE  = zeros(L*N,K,nbrOfRealizations);
w_PRZF    = zeros(L*N,K,nbrOfRealizations);

%%
% ======================================================================= %
% ================== COMBINING AND SCALING FACTORS ====================== %
% ======================================================================= %
% First generate the combining, then the scaling factors

%Go through all channel realizations
for n = 1:nbrOfRealizations

    % ------------------------------------------------------------------- %
    % MR and LP-MMSE schemes
    % ------------------------------------------------------------------- %
    % Go through all APs
    for l = 1:L

        % Extract which UEs are served by the AP
        servedUEs = find(D(l,:)==1);

        % Compute sum of estimation error covariance matrices of the UEs
        % served by AP l
        Cserved = sum(C(:,:,l,servedUEs),4);

        % Compute MR and LP-MMSE combining
        V_MR = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        V_LP_MMSE = p*((p*(V_MR(:,servedUEs)*V_MR(:,servedUEs)')+p*Cserved+eyeN)\V_MR(:,servedUEs));

        % Compute the precoding scaling factor by Monte Carlo methods
        scaling_MR(l,servedUEs) = scaling_MR(l,servedUEs) + sum(abs(V_MR(:,servedUEs)).^2,1)/nbrOfRealizations;
        scaling_LP_MMSE(l,servedUEs) = scaling_LP_MMSE(l,servedUEs) + sum(abs(V_LP_MMSE).^2,1)/nbrOfRealizations;

    end

    % ------------------------------------------------------------------- %
    % P-MMSE and P-RZF schemes
    % ------------------------------------------------------------------- %
    % Go through all UEs
    for k = 1:K

        % Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);

        % Determine which UEs that are served by partially the same set
        % of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;

        % Extract channel realizations and estimation error correlation
        % matrices for the APs that involved in the service of UE k
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);

        for l = 1:La
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end

        % Compute MMSE, P-MMSE and P-RZF combining
        V_P_MMSE = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        V_P_RZF = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));

        % Compute the precoding scaling factor by Monte Carlo methods
        scaling_P_MMSE(k) = scaling_P_MMSE(k) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_RZF(k) = scaling_P_RZF(k) + sum(abs(V_P_RZF).^2,1)/nbrOfRealizations;

        % Go through all the serving APs
        for l=1:La
            V_P_MMSE2 = V_P_MMSE((l-1)*N+1:l*N,:);
            V_P_RZF2 = V_P_RZF((l-1)*N+1:l*N,:);

            portionScaling_PMMSE(servingAPs(l),k) = portionScaling_PMMSE(servingAPs(l),k) ...
                + sum(abs(V_P_MMSE2).^2,1)/nbrOfRealizations;

            portionScaling_PRZF(servingAPs(l),k) = portionScaling_PRZF(servingAPs(l),k) ...
                + sum(abs(V_P_RZF2).^2,1)/nbrOfRealizations;
        end
    end

end

% Normalize the norm squares of the portions for the normalized
% centralized precoders
portionScaling_PMMSE = portionScaling_PMMSE./repmat(scaling_P_MMSE.',[L 1]);
portionScaling_PRZF = portionScaling_PRZF./repmat(scaling_P_RZF.',[L 1]);

% ======================================================================= %
% =================== CENTRALIZED POWER ALLOCATION ====================== %
% ======================================================================= %
% The parameters for the scalable centralized downlink power allocation in
% (7.43)
upsilon = -0.5;
kappa = 0.5;

% Compute the power allocation coefficients for centralized precoding
% according to (7.43), fractional power allocation
rho_PMMSE = functionCentralizedPowerAllocation(K, gainOverNoisedB, D, PowerDL_mW, portionScaling_PMMSE, upsilon, kappa);
rho_PRZF  = functionCentralizedPowerAllocation(K, gainOverNoisedB, D, PowerDL_mW,  portionScaling_PRZF, upsilon, kappa);

%%
% ======================================================================= %
% ======================= PRECODING VECTORS ============================= %
% ======================================================================= %
% Generates the precoding from combining vectors

% Go through all channel realizations
for n = 1:nbrOfRealizations
    % Matrix to store Monte-Carlo results for this realization
    interf_MR_n = zeros(K,K);
    interf_P_MMSE_n = zeros(K,K);
    interf_P_RZF_n = zeros(K,K);
    interf_LP_MMSE_n = zeros(K,K);

    % ------------------------------------------------------------------- %
    % MR and LP-MMSE schemes
    % ------------------------------------------------------------------- %
    % Go through all APs
    for l = 1:L

        % Extract channel realizations from all UEs to AP l
        Hallj = reshape(H((l-1)*N+1:l*N,n,:),[N K]);

        % Extract channel estimates from all UEs to AP l
        Hhatallj = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);

        % Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);

        % Compute sum of estimation error covariance matrices of the UEs
        % served by AP l
        Cserved = sum(C(:,:,l,servedUEs),4);

        % Compute MR combining
        V_MR = Hhatallj(:,servedUEs);

        % Compute LP-MMSE combining
        V_LP_MMSE = p*((p*(V_MR*V_MR')+p*Cserved+eyeN)\V_MR);

        % Go through all UEs served by the AP
        for ind = 1:length(servedUEs)

            % Extract UE index
            k = servedUEs(ind);

            % Normalize MR precoding
            w = V_MR(:,ind)*sqrt(rho_dist_mW(l,k)/scaling_MR(l,k));
            w_MR((l-1)*N+1:l*N,k,n) = w;

            % Compute gain of the signal from UE that arrives at other UEs
            interUserGains_MR(:,k,n) = interUserGains_MR(:,k,n) + Hallj'*w;

            signal_MR1(k) = signal_MR1(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            interf_MR_n(:,k) = interf_MR_n(:,k) + Hallj'*w;

            % Normalize LP-MMSE precoding
            w = V_LP_MMSE(:,ind)*sqrt(rho_dist_mW(l,k)/scaling_LP_MMSE(l,k));
            w_LP_MMSE((l-1)*N+1:l*N,k,n) = w;

            % Compute realizations of the terms inside the expectations
            % of the signal and interference terms in Corollary 6.3
            signal_LP_MMSE(k) = signal_LP_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            interf_LP_MMSE_n(:,k) = interf_LP_MMSE_n(:,k) + Hallj'*w;

            % Compute gain of the signal from UE that arrives at other UEs
            interUserGains_LP_MMSE(:,k,n) = interUserGains_LP_MMSE(:,k,n) + Hallj'*w;


        end

    end

    % ------------------------------------------------------------------- %
    % P-RZF and P-MMSE schemes
    % ------------------------------------------------------------------- %
    % Consider the centralized schemes
    % Go through all UEs
    for k = 1:K

        % Determine the set of serving APs
        servingAPs = find(D(:,k)==1);

        La = length(servingAPs);

        % Determine which UEs that are served by partially the same set
        % of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;

        % Extract channel realizations and estimation error correlation
        % matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);

        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);

        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        % Compute P-MMSE precoding
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));

        % Apply power allocation
        w = w*sqrt(rho_PMMSE(k)/scaling_P_MMSE(k));

        % MF Modification
        for index = 1:length(servingAPs)
            l = servingAPs(index);
            w_P_MMSE(((l-1)*N+1):(l*N),k,n) = w((index-1)*N+1:index*N);
        end

        % Compute realizations of the terms inside the expectations
        % of the signal and interference terms in Theorem 6.1
        signal_P_MMSE(k) = signal_P_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_MMSE_n(:,k) = interf_P_MMSE_n(:,k) + Hallj_active'*w;

        % Compute gain of the signal from UE that arrives at other UEs
        interUserGains_P_MMSE(:,k,n) = interUserGains_P_MMSE(:,k,n) + Hallj_active'*w;

        % Compute P-RZF combining
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));

        % Apply power allocation
        w = w*sqrt(rho_PRZF(k)/scaling_P_RZF(k));

        % MF Modification
        for index = 1:length(servingAPs)
            l = servingAPs(index);
            w_PRZF(((l-1)*N+1):(l*N),k,n) = w((index-1)*N+1:index*N);
        end

        % Compute realizations of the terms inside the expectations
        % of the signal and interference terms in Theorem 6.1
        signal_P_RZF(k) = signal_P_RZF(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_RZF_n(:,k) = interf_P_RZF_n(:,k) + Hallj_active'*w;

        % Compute gain of the signal from UE that arrives at other UEs
        interUserGains_P_RZF(:,k,n) = interUserGains_P_RZF(:,k,n) + Hallj_active'*w;

    end

    % Compute interference power in one realization
    interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;
    interf_P_RZF = interf_P_RZF + sum(abs(interf_P_RZF_n).^2,2)/nbrOfRealizations;

    interf_MR1 = interf_MR1 + sum(abs(interf_MR_n).^2,2)/nbrOfRealizations;
    interf_LP_MMSE = interf_LP_MMSE + sum(abs(interf_LP_MMSE_n).^2,2)/nbrOfRealizations;

end

% ======================================================================= %
% =================== DOWNLINK SPECTRAL EFFICIENCY ====================== %
% ======================================================================= %

% Compute SE in Corollary 6.3 with MR  using the closed-form expressions in Corollary 6.4
SE_MR = prelogFactor*real(log2(1+(abs(signal_MR1).^2) ./ (interf_MR1 - abs(signal_MR1).^2 + 1)));

% Compute SE  in Corollary 6.3 with LP-MMSE
SE_LP_MMSE = prelogFactor*real(log2(1+(abs(signal_LP_MMSE).^2) ./ (interf_LP_MMSE - abs(signal_LP_MMSE).^2 + 1)));

% Compute SE in Theorem 6.1 with P-MMSE
SE_P_MMSE = prelogFactor*real(log2(1+(abs(signal_P_MMSE).^2) ./ (interf_P_MMSE - abs(signal_P_MMSE).^2 + 1)));

% Compute SE in Theorem 6.1 with P-RZF
SE_P_RZF = prelogFactor*real(log2(1+(abs(signal_P_RZF).^2) ./ (interf_P_RZF - abs(signal_P_RZF).^2 + 1)));
