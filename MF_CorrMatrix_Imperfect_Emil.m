% This function computes the impferct channel coraviance of the channels
% between UEs and APs. It reproduces the method proposed in the paper
% "Massive MIMO with Imperfect Channel Covariance Information"

% INPUT:
% nbrOfRealizations = the number of channel realizations
% nbrOfAPs          = the number of APs
% nbrOfUEs          = the number of UEs
% N                 = the number of antenna elements per AP
% p                 = the uplink transmit power per UE
% tau_p             = the number of orthogonal pilots
% H                 = the real channel, comprising LOS and NLOS components
%                     Dim: nbrOfAPs*N x nbrOfRealizations x nbrOfUEs
% noiseVariance_dBm = the noise variance in dBm

% OUTPUT:
% gainOverNoisedB   = the imperfect channel gain normalized by noise
%                     variance of UE k regarding all APs normalized by the
%                     noise variance
%                     Dim: nbrOfAPs x 1
% R_imp             = imperfect spatial correlation matrix
%                     Dim: N x N x nbrOfAPs x nbrOfUEs

% PsiInv_imp        = imperfect spatial correlation matrix of UEs sharing
%                     the same pilot 
%                     Dim: N x N x nbrOfAPss x nbrOfUEs

% Note: In our simulations we consider N = 1

function [gainOverNoise_dB, R_imp] = MF_CorrMatrix_Imperfect_Emil(nbrOfRealizations, nbrOfAPs, nbrOfUEs, N, p, tau_p, H, noiseVariance_dBm)

% -------------------------------------------------------------------------
% PREPARE TO STORE RESULTS
% -------------------------------------------------------------------------

% Additional pilots
tau_p_aux_R = min(tau_p, nbrOfUEs);

% For the estimation of Q
pilotIndex_add_R = zeros(nbrOfUEs,1);
nbrOfRealizations_yp_n_R = nbrOfRealizations;
nbrOfRealizations_yp_k   = ceil(nbrOfRealizations_yp_n_R/2);

yp_n_R = zeros(N,nbrOfAPs,tau_p_aux_R,nbrOfRealizations_yp_n_R);
yp_n_TIMES_yp_n_H_R = zeros(N,N,nbrOfAPs,tau_p_aux_R,nbrOfRealizations_yp_n_R);
sum_yp_n_TIMES_yp_n_H_R = zeros(N,N,nbrOfAPs,tau_p_aux_R);

% For the large-scale fading estimation
gainOverNoise_lin = zeros(nbrOfAPs,nbrOfUEs);

% For the estimation of R
yp_n_k_TIMES_yp_n_k_H = zeros(N,N,nbrOfAPs,nbrOfUEs,nbrOfRealizations_yp_k);
sum_yp_n_k_TIMES_yp_n_k_H = zeros(N,N,nbrOfAPs,nbrOfUEs);

yp_n_k_minus = zeros(N,nbrOfAPs,nbrOfUEs,nbrOfRealizations_yp_k);
yp_n_k = zeros(N,nbrOfAPs,nbrOfUEs,nbrOfRealizations_yp_k);

noiseVariance_mW = 10^(noiseVariance_dBm/10);
Np = sqrt(0.5)*sqrt(noiseVariance_mW)*(randn(N,nbrOfRealizations_yp_n_R,nbrOfAPs,tau_p_aux_R) + 1i*randn(N,nbrOfRealizations_yp_n_R,nbrOfAPs,tau_p_aux_R));

% -------------------------------------------------------------------------
% COVARIANCE ESTIMATION "Q" (TOTAL)
% -------------------------------------------------------------------------
% Pilot assignment
pilotIndex_add_R(1:tau_p_aux_R) = 1:tau_p_aux_R;
if nbrOfUEs > tau_p_aux_R
    pilotIndex_add_R(tau_p_aux_R+1:nbrOfUEs) = randi(tau_p_aux_R,[1 nbrOfUEs-tau_p_aux_R]);
end

for realization = 1:nbrOfRealizations_yp_n_R
    for indexAP = 1:nbrOfAPs
        for t = 1:tau_p_aux_R % Go through all pilots
            yp_n_R(:,indexAP,t,realization) = sqrt(p)*sqrt(tau_p_aux_R)*sum(H((indexAP-1)*N+1:indexAP*N,realization,t==pilotIndex_add_R),3) + Np(:,realization, indexAP,t);
            yp_n_R(:,indexAP,t,realization) = yp_n_R(:,indexAP,t,realization)/(sqrt(p)*sqrt(tau_p_aux_R));
            yp_n_TIMES_yp_n_H_R(:,:,indexAP,t,realization) = yp_n_R(:,indexAP,t,realization)*transpose(conj(yp_n_R(:,indexAP,t,realization)));
            sum_yp_n_TIMES_yp_n_H_R(:,:,indexAP,t) = sum_yp_n_TIMES_yp_n_H_R(:,:,indexAP,t) + yp_n_TIMES_yp_n_H_R(:,:,indexAP,t,realization)/nbrOfRealizations_yp_n_R;
        end
    end
end

PsiInv_imp = sum_yp_n_TIMES_yp_n_H_R;

% -------------------------------------------------------------------------
% COVARIANCE ESTIMATION "R" (INDIVIDUAL)
% -------------------------------------------------------------------------
for indexAP = 1:nbrOfAPs
    for t = 1:tau_p_aux_R % Go through all pilots
        servedUEs_id = find(pilotIndex_add_R == t);
        if length(servedUEs_id) == 1
            sum_yp_n_k_TIMES_yp_n_k_H(:,:,indexAP,servedUEs_id) = PsiInv_imp(:,:,indexAP,t);
        else
            for indexUE = 1:length(servedUEs_id)
                servedUEs_id_aux = setdiff(servedUEs_id,servedUEs_id(indexUE));
                for realization = 1:nbrOfRealizations_yp_k
                    yp_n_k_minus(:,indexAP,servedUEs_id(indexUE),realization) = sqrt(p)*sqrt(tau_p_aux_R)*sum(H((indexAP-1)*N+1:indexAP*N,realization,servedUEs_id_aux),3) + Np(:,realization,indexAP,t);
                    yp_n_k_minus(:,indexAP,servedUEs_id(indexUE),realization) = yp_n_k_minus(:,indexAP,servedUEs_id(indexUE),realization)/(sqrt(p)*sqrt(tau_p_aux_R)); % Eq (4)
                    yp_n_k(:,indexAP,servedUEs_id(indexUE),realization) = yp_n_R(:,indexAP,t,realization)-yp_n_k_minus(:,indexAP,servedUEs_id(indexUE),realization);
                    yp_n_k_TIMES_yp_n_k_H(:,:,indexAP,servedUEs_id(indexUE),realization) = yp_n_k(:,indexAP,servedUEs_id(indexUE),realization)*transpose(conj(yp_n_k(:,indexAP,servedUEs_id(indexUE),realization)));
                    sum_yp_n_k_TIMES_yp_n_k_H(:,:,indexAP,servedUEs_id(indexUE)) = sum_yp_n_k_TIMES_yp_n_k_H(:,:,indexAP,servedUEs_id(indexUE)) + yp_n_k_TIMES_yp_n_k_H(:,:,indexAP,servedUEs_id(indexUE),realization)/nbrOfRealizations_yp_k;
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% IMPERFECT INDIVIDUAL COVARIANCE MATRIX
% -------------------------------------------------------------------------
R_imp = sum_yp_n_k_TIMES_yp_n_k_H;

clear H Np yp_n yp_n_TIMES_yp_n_H sum_yp_n_TIMES_yp_n_H servedUEs_id servedUEs_id_aux yp_n_k yp_n_k_TIMES_yp_n_k_H;

% For N> 1
% % Regularization presented in Eq. (10) of the reference paper
% neta = 0.25;
% mi = 0.25;
% PsiInv_imp = neta*PsiInv_imp_sample + (1-neta)*PsiInv_imp_diagonal;
% R_imp = mi*R_imp_sample + (1-mi)*R_imp_diagonal;
% A = 0;

% -------------------------------------------------------------------------
% IMPERFECT LARGE SCALE FADING, GAIN OVER NOISE
% -------------------------------------------------------------------------
% Compute the imperfect gain over noise, i.e., large scale fading

for indexAP = 1:nbrOfAPs
    for indexUE = 1:nbrOfUEs
        gainOverNoise_lin(indexAP,indexUE) = (1/N)*abs(trace(R_imp(:,:,indexAP,indexUE))) + eps;
    end
end

gainOverNoise_dB = pow2db(gainOverNoise_lin);
