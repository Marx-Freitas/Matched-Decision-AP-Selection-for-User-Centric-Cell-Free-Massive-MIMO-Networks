% function [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
%     SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
%     Rate_MR, Rate_L_MMSE, Rate_LP_MMSE,...
%     Rate_MMSE, Rate_P_MMSE, Rate_P_RZF, rho_AP_mmse, rho_AP_pmmse, rho_AP_przf] ...
%     = MF_functionComputeSE_downlinkWithRician(Hhat,H,D,C,Hhat_nfl,C_nfl,...
%     tau_c,tau_p,nbrOfRealizations,N,K,L,p,rho_dist,gainOverNoisedB,rho_tot,BW,...
%     bits_per_sample_up,acceptable_deg_data, bits_downlink_dist, bits_downlink_cent)

function [SE_P_MMSE, SE_P_RZF,  ...
    SE_LP_MMSE, SE_MR, ...
    Rate_MR, Rate_LP_MMSE,...
    Rate_P_MMSE, Rate_P_RZF, rho_AP_pmmse, rho_AP_przf] ...
    = MF_functionComputeSE_downlinkWithRician(Hhat,H,D,C,Hhat_nfl,C_nfl,...
    tau_c,tau_p,nbrOfRealizations,N,K,L,p,rho_dist,gainOverNoisedB,rho_tot,BW,...
    bits_per_sample_up,acceptable_deg_data, bits_downlink_dist, bits_downlink_cent)
%Compute downlink SE for different transmit precoding schemes using the capacity
%bound in Theorem 6.1 for the centralized schemes and the capacity bound
%in Corollary 6.3 for the distributed schemes. Compute the genie-aided
%downlink SE from Corollary 6.6 for the centralized and the distributed operations.
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
%D                 = DCC matrix for cell-free setup with dimension L x K
%                    where (l,k) is one if AP l serves UE k and zero otherwise
%B                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel estimate
%                    between AP l and UE k, normalized by noise variance
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k,
%                    normalized by noise variance
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Number of UEs
%L                 = Number of APs
%p                 = Uplink transmit power per UE (same for everyone)
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k,
%                    normalized by noise
%pilotIndex        = Vector containing the pilot assigned to each UE
%rho_dist          = Matrix with dimension L x K where (l,k) is the power
%                    allocated to UE k by AP l in the distributed downlink
%                    operation
%gainOverNoisedB   = Matrix with dimension L x K where (l,k) is the channel
%                    gain (normalized by the noise variance) between AP l
%                    and UE k
%rho_tot           = Maximum allowed transmit power for each AP
%
%OUTPUT:
%SE_MMSE           = SEs achieved with MMSE precoding in (6.16)
%SE_P_MMSE         = SEs achieved with P-MMSE precoding in (6.17)
%SE_P_RZF          = SEs achieved with P-RZF precoding in (6.18)
%SE_L_MMSE         = SEs achieved with L-MMSE precoding in (6.25)
%SE_LP_MMSE        = SEs achieved with LP-MMSE precoding in (6.33)
%SE_MR             = SEs achieved with MR precoding in (6.26)
%Gen_SE_P_MMSE     = Genie-aided SEs achieved with P-MMSE precoding in (6.17)
%Gen_SE_P_RZF      = Genie-aided SEs achieved with P-RZF precoding in (6.18)
%Gen_SE_LP_MMSE    = Genie-aided SEs achieved with LP-MMSE precoding in (6.33)
%Gen_SE_MR         = Genie-aided SEs achieved with MR precoding in (6.26)
%
%
%This Matlab function was developed to generate simulation results to:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO",
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Store the N x N identity matrix
eyeN = eye(N);

%Compute the prelog factor assuming only downlink data transmission
prelogFactor = (1-tau_p/tau_c);


%Prepare to store the terms that appear in SEs

scaling_MR = zeros(L,K);
signal_MR1 = zeros(K,1);
interf_MR1 = zeros(K,1);
interUserGains_MR = zeros(K,K,nbrOfRealizations);

%signal_MMSE = zeros(K,1);
%interf_MMSE = zeros(K,1);
scaling_MMSE = zeros(K,1);
portionScaling_MMSE = zeros(L,K);
%interUserGains_MMSE = zeros(K,K,nbrOfRealizations);

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

%signal_MMSE_nfl = zeros(K,1);
%interf_MMSE_nfl = zeros(K,1);
scaling_MMSE_nfl = zeros(K,1);
portionScaling_MMSE_nfl = zeros(L,K);

signal_P_MMSE_nfl = zeros(K,1);
interf_P_MMSE_nfl = zeros(K,1);
scaling_P_MMSE_nfl = zeros(K,1);
portionScaling_PMMSE_nfl = zeros(L,K);

signal_P_RZF_nfl = zeros(K,1);
interf_P_RZF_nfl = zeros(K,1);
scaling_P_RZF_nfl = zeros(K,1);
portionScaling_PRZF_nfl = zeros(L,K);

%signal_L_MMSE = zeros(K,1);
%interf_L_MMSE = zeros(K,1);
scaling_L_MMSE = zeros(L,K);
%interUserGains_L_MMSE = zeros(K,K,nbrOfRealizations);

signal_LP_MMSE = zeros(K,1);
interf_LP_MMSE = zeros(K,1);
scaling_LP_MMSE = zeros(L,K);
interUserGains_LP_MMSE = zeros(K,K,nbrOfRealizations);

%Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %L-MMSE, LP-MMSE precoding
    for l = 1:L
        
        %Extract which UEs are served by the AP
        servedUEs = find(D(l,:)==1);
        
        %Compute sum of estimation error covariance matrices of the UEs
        %served by AP l
        Cserved = sum(C_nfl(:,:,l,servedUEs),4);
        
        %Compute L-MMSE and LP-MMSE precoding
        V_MR = reshape(Hhat_nfl((l-1)*N+1:l*N,n,:),[N K]);
        V_L_MMSE = p*((p*(V_MR*V_MR')+p*sum(C_nfl(:,:,l,:),4)+eyeN)\V_MR(:,servedUEs));
        V_LP_MMSE = p*((p*(V_MR(:,servedUEs)*V_MR(:,servedUEs)')+p*Cserved+eyeN)\V_MR(:,servedUEs));
        
        %Compute scaling factor by Monte Carlo methods
        scaling_MR(l,servedUEs) = scaling_MR(l,servedUEs) + sum(abs(V_MR(:,servedUEs)).^2,1)/nbrOfRealizations;
        scaling_L_MMSE(l,servedUEs) = scaling_L_MMSE(l,servedUEs) + sum(abs(V_L_MMSE).^2,1)/nbrOfRealizations;
        scaling_LP_MMSE(l,servedUEs) = scaling_LP_MMSE(l,servedUEs) + sum(abs(V_LP_MMSE).^2,1)/nbrOfRealizations;
        
    end
    
    %MMSE, P-MMSE, and P-RZF precoding
    
    
    %Go through all UEs
    for k = 1:K
        
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);
        
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        
        %Compute MMSE, P-MMSE, and P-RZF precoding
        V_MMSE = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        V_P_MMSE = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        V_P_RZF = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));
        
        %Compute scaling factor by Monte Carlo methods
        scaling_MMSE(k) = scaling_MMSE(k) + sum(abs(V_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_MMSE(k) = scaling_P_MMSE(k) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_RZF(k) = scaling_P_RZF(k) + sum(abs(V_P_RZF).^2,1)/nbrOfRealizations;
        
        %_nfl
        Hhatallj_active_nfl = zeros(N*La,K);
        C_tot_blk_nfl = zeros(N*La,N*La);
        C_tot_blk_partial_nfl = zeros(N*La,N*La);
        
        for l = 1:La
            Hhatallj_active_nfl((l-1)*N+1:l*N,:) = reshape(Hhat_nfl((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk_nfl((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C_nfl(:,:,servingAPs(l),:),4);
            C_tot_blk_partial_nfl((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C_nfl(:,:,servingAPs(l),servedUEs),4);
        end
        
        %Compute MMSE, P-MMSE, and P-RZF precoding
        V_MMSE_nfl = p*((p*(Hhatallj_active_nfl*Hhatallj_active_nfl')+p*C_tot_blk_nfl+eye(La*N))\Hhatallj_active_nfl(:,k));
        V_P_MMSE_nfl = p*((p*(Hhatallj_active_nfl(:,servedUEs)*Hhatallj_active_nfl(:,servedUEs)')+p*C_tot_blk_partial_nfl+eye(La*N))\Hhatallj_active_nfl(:,k));
        V_P_RZF_nfl = p*((p*(Hhatallj_active_nfl(:,servedUEs)*Hhatallj_active_nfl(:,servedUEs)')+eye(La*N))\Hhatallj_active_nfl(:,k));
        
        %Compute scaling factor by Monte Carlo methods
        scaling_MMSE_nfl(k) = scaling_MMSE_nfl(k) + sum(abs(V_MMSE_nfl).^2,1)/nbrOfRealizations;
        scaling_P_MMSE_nfl(k) = scaling_P_MMSE_nfl(k) + sum(abs(V_P_MMSE_nfl).^2,1)/nbrOfRealizations;
        scaling_P_RZF_nfl(k) = scaling_P_RZF_nfl(k) + sum(abs(V_P_RZF_nfl).^2,1)/nbrOfRealizations;
        
        %Go through all the serving APs
        for l=1:La
            
            %Extract the portions of the centralized precoding vectors
            V_MMSE2 = V_MMSE((l-1)*N+1:l*N,:);
            V_P_MMSE2 = V_P_MMSE((l-1)*N+1:l*N,:);
            V_P_RZF2 = V_P_RZF((l-1)*N+1:l*N,:);
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Theorem 6.1
            
            portionScaling_MMSE(servingAPs(l),k) = portionScaling_MMSE(servingAPs(l),k) ...
                + sum(abs(V_MMSE2).^2,1)/nbrOfRealizations;
            
            portionScaling_PMMSE(servingAPs(l),k) = portionScaling_PMMSE(servingAPs(l),k) ...
                + sum(abs(V_P_MMSE2).^2,1)/nbrOfRealizations;
            
            portionScaling_PRZF(servingAPs(l),k) = portionScaling_PRZF(servingAPs(l),k) ...
                + sum(abs(V_P_RZF2).^2,1)/nbrOfRealizations;
            
            
            %_nfl
            %Extract the portions of the centralized precoding vectors
            V_MMSE2_nfl = V_MMSE_nfl((l-1)*N+1:l*N,:);
            V_P_MMSE2_nfl = V_P_MMSE_nfl((l-1)*N+1:l*N,:);
            V_P_RZF2_nfl = V_P_RZF_nfl((l-1)*N+1:l*N,:);
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Theorem 6.1
            
            portionScaling_MMSE_nfl(servingAPs(l),k) = portionScaling_MMSE_nfl(servingAPs(l),k) ...
                + sum(abs(V_MMSE2_nfl).^2,1)/nbrOfRealizations;
            
            portionScaling_PMMSE_nfl(servingAPs(l),k) = portionScaling_PMMSE_nfl(servingAPs(l),k) ...
                + sum(abs(V_P_MMSE2_nfl).^2,1)/nbrOfRealizations;
            
            portionScaling_PRZF_nfl(servingAPs(l),k) = portionScaling_PRZF_nfl(servingAPs(l),k) ...
                + sum(abs(V_P_RZF2_nfl).^2,1)/nbrOfRealizations;
        end
    end
    
end

%Normalize the norm squares of the portions for the normalized centralized precoders
portionScaling_MMSE = portionScaling_MMSE./repmat(scaling_MMSE.',[L 1]);

portionScaling_PMMSE = portionScaling_PMMSE./repmat(scaling_P_MMSE.',[L 1]);

portionScaling_PRZF = portionScaling_PRZF./repmat(scaling_P_RZF.',[L 1]);

%_nfl
%Normalize the norm squares of the portions for the normalized centralized precoders
portionScaling_MMSE_nfl = portionScaling_MMSE_nfl./repmat(scaling_MMSE_nfl.',[L 1]);

portionScaling_PMMSE_nfl = portionScaling_PMMSE_nfl./repmat(scaling_P_MMSE_nfl.',[L 1]);

portionScaling_PRZF_nfl = portionScaling_PRZF_nfl./repmat(scaling_P_RZF_nfl.',[L 1]);

%% The parameters for the scalable centralized downlink power allocation in
%(7.43)
upsilon = -0.5;
kappa = 0.5;

%Compute the power allocation coefficients for centralized precoding
%according to (7.43)

rho_MMSE = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_MMSE,upsilon,kappa);

rho_PMMSE = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PMMSE,upsilon,kappa);

rho_PRZF = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PRZF,upsilon,kappa);

%_nfl
rho_MMSE_nfl = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_MMSE_nfl,upsilon,kappa);

rho_PMMSE_nfl = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PMMSE_nfl,upsilon,kappa);

rho_PRZF_nfl = functionCentralizedPowerAllocation(K,gainOverNoisedB,D,rho_tot,portionScaling_PRZF_nfl,upsilon,kappa);

%% Desired signal and interference calculations
w_pmmse=zeros(L*N,K,nbrOfRealizations);
w_przf=zeros(L*N,K,nbrOfRealizations);
w_mmse=zeros(L*N,K,nbrOfRealizations);
antenna_AP_map=reshape(1:L*N,[N L]);
for n = 1:nbrOfRealizations
    
    %Matrix to store Monte-Carlo results for this realization
    interf_MR_n = zeros(K,K);
    interf_MMSE_n = zeros(K,K);
    interf_P_MMSE_n = zeros(K,K);
    interf_P_RZF_n = zeros(K,K);
    % interf_MMSE_n_nfl = zeros(K,K);
    interf_P_MMSE_n_nfl = zeros(K,K);
    interf_P_RZF_n_nfl = zeros(K,K);
    % interf_L_MMSE_n = zeros(K,K);
    interf_LP_MMSE_n = zeros(K,K);
    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract channel estimates from all UEs to AP l
        Hhatallj = reshape(Hhat_nfl((l-1)*N+1:l*N,n,:),[N K]);
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute sum of estimation error covariance matrices of the UEs
        %served by AP l
        Cserved = sum(C_nfl(:,:,l,servedUEs),4);
        
        %Compute MR combining
        V_MR = Hhatallj(:,servedUEs);
        
        % Compute L-MMSE combining
        % V_L_MMSE = p*((p*(Hhatallj*Hhatallj')+p*sum(C_nfl(:,:,l,:),4)+eyeN)\V_MR);
        
        %Compute LP-MMSE combining
        V_LP_MMSE = p*((p*(V_MR*V_MR')+p*Cserved+eyeN)\V_MR);
        
        
        %Go through all UEs served by the AP
        for ind = 1:length(servedUEs)
            
            %Extract UE index
            k = servedUEs(ind);
            
            %Normalize MR precoding
            w = V_MR(:,ind)*sqrt(rho_dist(l,k)/scaling_MR(l,k));
            
            %Compute gain of the signal from UE that arrives at other UEs
            interUserGains_MR(:,k,n) = interUserGains_MR(:,k,n) + Hallj'*w;
            
            signal_MR1(k) = signal_MR1(k) + (Hallj(:,k)'*w)/nbrOfRealizations;%alpha(quantization) is added latter
            interf_MR_n(:,k) = interf_MR_n(:,k) + (Hallj'*w);%alpha(quantization) is added latter
            % %Normalize L-MMSE precoding
            % w = V_L_MMSE(:,ind)*sqrt(rho_dist(l,k)/scaling_L_MMSE(l,k));
            
            % %Compute gain of the signal from UE that arrives at other UEs
            % interUserGains_L_MMSE(:,k,n) = interUserGains_L_MMSE(:,k,n) + Hallj'*w;
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Corollary 6.3
            %signal_L_MMSE(k) = signal_L_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;%alpha(quantization) is added latter
            %interf_L_MMSE_n(:,k) = interf_L_MMSE_n(:,k) + (Hallj'*w);%alpha(quantization) is added latter
            
            %Normalize LP-MMSE precoding
            w = V_LP_MMSE(:,ind)*sqrt(rho_dist(l,k)/scaling_LP_MMSE(l,k));
            
            %Compute gain of the signal from UE that arrives at other UEs
            interUserGains_LP_MMSE(:,k,n) = interUserGains_LP_MMSE(:,k,n) + Hallj'*w;
           
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Corollary 6.3
            signal_LP_MMSE(k) = signal_LP_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;%alpha(quantization) is added latter
            interf_LP_MMSE_n(:,k) = interf_LP_MMSE_n(:,k) + (Hallj'*w);%alpha(quantization) is added latter
            
            
        end
        
    end  
    
    %Consider the centralized schemes
    
    
    %Go through all UEs
    for k = 1:K
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        servingAntennas=antenna_AP_map(:,servingAPs);
        servingAntennas=servingAntennas(:);
        La = length(servingAPs);
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
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
        %Compute P-MMSE precoding
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        %Apply power allocation
        w_pmmse(servingAntennas,k,n) = w*sqrt(rho_PMMSE(k)/scaling_P_MMSE(k));
        %Compute P-RZF combining
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));
        %Apply power allocation
        w_przf(servingAntennas,k,n) = w*sqrt(rho_PRZF(k)/scaling_P_RZF(k));
        
        %%Compute MMSE combining
        %w = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        %%Apply power allocation
        %w_mmse(servingAntennas,k,n) = w*sqrt(rho_MMSE(k)/scaling_MMSE(k));

        %Compute realizations of the terms inside the expectations for P-MMSE
        signal_P_MMSE(k) = signal_P_MMSE(k) + (Hallj_active(:,k)'*w_pmmse(servingAntennas,k,n))/nbrOfRealizations;
        interf_P_MMSE_n(:,k) = interf_P_MMSE_n(:,k) + (Hallj_active'*w_pmmse(servingAntennas,k,n));
        %Compute realizations of the terms inside the expectations for P-RZF
        signal_P_RZF(k) = signal_P_RZF(k) + (Hallj_active(:,k)'*w_przf(servingAntennas,k,n))/nbrOfRealizations;
        interf_P_RZF_n(:,k) = interf_P_RZF_n(:,k) + (Hallj_active'*w_przf(servingAntennas,k,n));
        
        %%Compute realizations of the terms inside the expectations for MMSE
        %signal_MMSE(k) = signal_MMSE(k) + (Hallj_active(:,k)'*w_mmse(servingAntennas,k,n))/nbrOfRealizations;
        %interf_MMSE_n(:,k) = interf_MMSE_n(:,k) + (Hallj_active'*w_mmse(servingAntennas,k,n));
        
        %_nfl
        Hhatallj_active_nfl = zeros(N*La,K);
        C_tot_blk_nfl = zeros(N*La,N*La);
        C_tot_blk_partial_nfl = zeros(N*La,N*La);
        for l = 1:La
            Hhatallj_active_nfl((l-1)*N+1:l*N,:) = reshape(Hhat_nfl((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk_nfl((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C_nfl(:,:,servingAPs(l),:),4);
            C_tot_blk_partial_nfl((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C_nfl(:,:,servingAPs(l),servedUEs),4);
        end
        %Compute P-MMSE precoding
        w = p*((p*(Hhatallj_active_nfl(:,servedUEs)*Hhatallj_active_nfl(:,servedUEs)')+p*C_tot_blk_partial_nfl+eye(La*N))\Hhatallj_active_nfl(:,k));
        %Apply power allocation
        w = w*sqrt(rho_PMMSE_nfl(k)/scaling_P_MMSE_nfl(k));
        interUserGains_P_MMSE(:,k,n) = interUserGains_P_MMSE(:,k,n) + Hallj_active'*w;
        signal_P_MMSE_nfl(k) = signal_P_MMSE_nfl(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_MMSE_n_nfl(:,k) = interf_P_MMSE_n_nfl(:,k) + (Hallj_active'*w);
        %Compute P-RZF combining
        w = p*((p*(Hhatallj_active_nfl(:,servedUEs)*Hhatallj_active_nfl(:,servedUEs)')+eye(La*N))\Hhatallj_active_nfl(:,k));
        %Apply power allocation
        w = w*sqrt(rho_PRZF_nfl(k)/scaling_P_RZF_nfl(k));
        interUserGains_P_RZF(:,k,n) = interUserGains_P_RZF(:,k,n) + Hallj_active'*w;
        signal_P_RZF_nfl(k) = signal_P_RZF_nfl(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        interf_P_RZF_n_nfl(:,k) = interf_P_RZF_n_nfl(:,k) + (Hallj_active'*w);
        
        %%Compute MMSE combining
        % w = p*((p*(Hhatallj_active_nfl*Hhatallj_active_nfl')+p*C_tot_blk_nfl+eye(La*N))\Hhatallj_active_nfl(:,k));
        %Apply power allocation
        
        %w = w*sqrt(rho_MMSE_nfl(k)/scaling_MMSE_nfl(k));
        % interUserGains_MMSE(:,k,n) = interUserGains_MMSE(:,k,n) + Hallj_active'*w;
        %signal_MMSE_nfl(k) = signal_MMSE_nfl(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        %interf_MMSE_n_nfl(:,k) = interf_MMSE_n_nfl(:,k) + (Hallj_active'*w);
    end
    
    %%Compute interference power in one realization
    %interf_MMSE = interf_MMSE + sum(abs(interf_MMSE_n).^2,2)/nbrOfRealizations;
    interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;
    interf_P_RZF = interf_P_RZF + sum(abs(interf_P_RZF_n).^2,2)/nbrOfRealizations;
    
    %interf_MMSE_nfl = interf_MMSE_nfl + sum(abs(interf_MMSE_n_nfl).^2,2)/nbrOfRealizations;
    interf_P_MMSE_nfl = interf_P_MMSE_nfl + sum(abs(interf_P_MMSE_n_nfl).^2,2)/nbrOfRealizations;
    interf_P_RZF_nfl = interf_P_RZF_nfl + sum(abs(interf_P_RZF_n_nfl).^2,2)/nbrOfRealizations;
    
    interf_MR1 = interf_MR1 + sum(abs(interf_MR_n).^2,2)/nbrOfRealizations;
    %interf_L_MMSE = interf_L_MMSE + sum(abs(interf_L_MMSE_n).^2,2)/nbrOfRealizations;
    interf_LP_MMSE = interf_LP_MMSE + sum(abs(interf_LP_MMSE_n).^2,2)/nbrOfRealizations;
end

%% fiding the centralized precoder statistics
%q_w_mmse=zeros(L*N,K);
q_w_przf=zeros(L*N,K);
q_w_pmmse=zeros(L*N,K);
for k=1:K
    %q_w_mmse(:,k)=mean(abs(squeeze(w_mmse(:,k,:))).^2,2);
    q_w_przf(:,k)=mean(abs(squeeze(w_przf(:,k,:))).^2,2);
    q_w_pmmse(:,k)=mean(abs(squeeze(w_pmmse(:,k,:))).^2,2);
end
%rho_AP_mmse=sum(reshape(sum(q_w_mmse,2),[N length(sum(q_w_mmse,2))/N]))';
rho_AP_pmmse=sum(reshape(sum(q_w_pmmse,2),[N length(sum(q_w_pmmse,2))/N]))';
rho_AP_przf=sum(reshape(sum(q_w_przf,2),[N length(sum(q_w_przf,2))/N]))';
%q_w_mmse=sum(q_w_mmse,2);
q_w_przf=sum(q_w_przf,2);
q_w_pmmse=sum(q_w_pmmse,2);
D_expanded=AF_repmat_bellow_row(D,length(H(:,1))/length(D(:,1)));
interf_P_RZF_quant=transpose(sum(q_w_przf.*mean(abs(D_expanded.*permute(H,[1 3 2])).^2,3)));
interf_P_MMSE_quant=transpose(sum(q_w_pmmse.*mean(abs(D_expanded.*permute(H,[1 3 2])).^2,3)));

% interf_MMSE_quant=transpose(sum(q_w_mmse.*mean(abs(D_expanded.*permute(H,[1 3 2])).^2,3)));

%
signal_P_MMSE=(abs(signal_P_MMSE).^2);
signal_P_RZF=(abs(signal_P_RZF).^2);

% signal_MMSE_nfl=(abs(signal_MMSE_nfl).^2);
% signal_P_MMSE_nfl=(abs(signal_P_MMSE_nfl).^2);
% signal_P_RZF_nfl=(abs(signal_P_RZF_nfl).^2);

signal_MR1=(abs(signal_MR1).^2);
signal_LP_MMSE=(abs(signal_LP_MMSE).^2);


% signal_L_MMSE=(abs(signal_L_MMSE).^2);
% signal_MMSE=(abs(signal_MMSE).^2);


%% Compute the SEs without fronthaul limitation
SE_MR_nfl = prelogFactor*real(log2(1+ signal_MR1 ./ (interf_MR1 - signal_MR1 + 1)));

%Compute SE  in Corollary 6.3 with LP-MMSE
SE_LP_MMSE_nfl = prelogFactor*real(log2(1+signal_LP_MMSE ./ (interf_LP_MMSE - signal_LP_MMSE + 1)));

%Compute SE in Theorem 6.1 with P-MMSE
SE_P_MMSE_nfl = prelogFactor*real(log2(1+signal_P_MMSE ./ (interf_P_MMSE - signal_P_MMSE + 1)));

%Compute SE in Theorem 6.1 with P-RZF
SE_P_RZF_nfl = prelogFactor*real(log2(1+signal_P_RZF ./ (interf_P_RZF - signal_P_RZF + 1)));

% % %Compute SE in Theorem 6.1 with MMSE
% % SE_MMSE_nfl = prelogFactor*real(log2(1+signal_MMSE_nfl ./ (interf_MMSE_nfl - signal_MMSE_nfl + 1)));
% % 
% % %Compute SE in Corollary 6.3 with L-MMSE
% % SE_L_MMSE_nfl = prelogFactor*real(log2(1+signal_L_MMSE ./ (interf_L_MMSE - signal_L_MMSE + 1)));

%% Compute the SEs and obtain fronthaul requirements
%Compute SE in Corollary 6.3 with MR  using the closed-form expressions in Corollary 6.4
% SE_MR = prelogFactor*real(log2(1+ (abs(signal_MR).^2)./ (interf_MR + sum(cont_MR,2) - (abs(signal_MR).^2) + 1)));
[SE_MR, bitMR]=AF_SE_and_FrontBitReq_dist(K, signal_MR1, interf_MR1, prelogFactor,...
    acceptable_deg_data, bits_downlink_dist);

%Compute SE  in Corollary 6.3 with LP-MMSE
[SE_LP_MMSE, bitLP_MMSE]=AF_SE_and_FrontBitReq_dist(K, signal_LP_MMSE, interf_LP_MMSE,...
    prelogFactor, acceptable_deg_data, bits_downlink_dist);

%Compute SE in Theorem 6.1 with P-MMSE
[SE_P_MMSE, bit_P_MMSE]=AF_SE_and_FrontBitReq_cent(signal_P_MMSE, interf_P_MMSE,...
    interf_P_MMSE_quant, prelogFactor, acceptable_deg_data, bits_downlink_cent);

%Compute SE in Theorem 6.1 with P-RZF
[SE_P_RZF, bit_P_RZF]=AF_SE_and_FrontBitReq_cent(signal_P_RZF, interf_P_RZF,...
    interf_P_RZF_quant, prelogFactor, acceptable_deg_data, bits_downlink_cent);

Rate_MR=2*sum((((tau_c-tau_p)*BW/tau_c)*bitMR).*D');
Rate_LP_MMSE=2*sum((((tau_c-tau_p)*BW/tau_c)*bitLP_MMSE).*D');
Rate_P_MMSE=2*BW*bit_P_MMSE*N*(tau_c-tau_p)/tau_c+2*BW*bits_per_sample_up*N*(tau_p)/tau_c;
Rate_P_RZF=2*BW*bit_P_RZF*N*(tau_c-tau_p)/tau_c+2*BW*bits_per_sample_up*N*(tau_p)/tau_c;


% %Compute SE in Theorem 6.1 with MMSE
% [SE_MMSE, bit_MMSE]=AF_SE_and_FrontBitReq_cent(signal_MMSE, interf_MMSE,...
%     interf_MMSE_quant, prelogFactor, acceptable_deg_data, bits_downlink_cent);
% 
% %Compute SE in Corollary 6.3 with L-MMSE
% [SE_L_MMSE, bitL_MMSE]=AF_SE_and_FrontBitReq_dist(K, signal_L_MMSE, interf_L_MMSE,...
%     prelogFactor, acceptable_deg_data, bits_downlink_dist);
% 
% Rate_L_MMSE=2*sum((((tau_c-tau_p)*BW/tau_c)*bitL_MMSE).*D');
% Rate_MMSE=2*BW*bit_MMSE*N*(tau_c-tau_p)/tau_c+2*BW*bits_per_sample_up*N*(tau_p)/tau_c;


end