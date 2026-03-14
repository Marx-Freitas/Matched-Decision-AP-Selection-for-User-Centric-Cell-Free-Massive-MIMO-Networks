% This Matlab function was developed to generate simulation results to:
% Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
% "Foundations of User-Centric Cell-Free Massive MIMO",
% Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
% pp 162-472. DOI: 10.1561/2000000109
%
% This is version 1.0 (Last edited: 2021-01-31)
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.

% THIS CODE WAS ALSO EDITED BY Marx M. M. Freitas, Daynara D. Souza, and
% André L. P. Fernandes in some parts

function [Hhat,C] = functionChannelEstimatesWithRician(R, Hmean, HmeanFase, H, nbrOfRealizations,...
    L, K , N, tau_p, pilotIndex, p, estimationMethod, covariance_matrix, R_imp)
% Generate the channel realizations and estimates of these channels for all
% UEs in the entire network. The channels are assumed to be correlated
% Although the function name contains the term Rician, the small-scale
% fading in this code can be simulated as both Rician and Rayleigh fading

% INPUT:
% R                  = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                     the spatial correlation matrix between AP l and UE k,
%                     normalized by noise variance
% Hmean             = the LOS component scaled by the Rician factor
%                     Dim: nbrOfAPs*N x nbrOfUEs
% HmeanFase         = the LOS component Hmean scaled by phase shifts
% nbrOfRealizations  = Number of channel realizations
% L                  = Number of APs
% K                  = Number of UEs in the network
% N                  = Number of antennas per AP
% tau_p              = Number of orthogonal pilots
% pilotIndex         = Vector containing the pilot assigned to each UE
% p                  = Uplink transmit power per UE (same for everyone)
% estimationMethod  = the channel estimation method
% covariance_matrix = the covariance matrix estimation method, i.e.,
%                     perfect (closed form expression) or imperfect
%                     (without a closed-form expression)
%                     Dim: N x N x nbrOfAPs x tau_p
% R_imp             = imperfect spatial correlation matrix
%                     Dim: N x N x nbrOfAPs x nbrOfUEs

% OUTPUT:
% Hhat         = Matrix with dimension L*N x nbrOfRealizations x K where
%               (:,n,k) is the estimated collective channel to UE k in
%               channel realization n.
% C            = Matrix with dimension N x N x L x K where (:,:,l,k) is the
%               spatial correlation matrix of the channel estimation error
%               between AP l and UE k, normalized by noise variance

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));

%Prepare to store results
Hhat = zeros(L*N,nbrOfRealizations,K);

C = zeros(N,N,L,K);

% Estimating the correlation matrices of UEs sharing the same pilot
if strcmp(covariance_matrix, 'imperfect')
    Q_imp = Psi_imp_estimation(N, L, tau_p, nbrOfRealizations, p, H, Np, pilotIndex);
end

%% Perform channel estimation
switch estimationMethod

    % =================================================================== %
    % ======================= Phase Aware MMSE ========================== %
    % =================================================================== %

    % Phase-aware MMSE
    case 'phaseAwareMMSE'

        %Go through all APs
        for l = 1:L

            %Go through all pilots
            for t = 1:tau_p

                %Compute processed pilot signal for all UEs that use pilot t
                %according to (4.4) with an additional scaling factor \sqrt{tau_p}
                yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);
                % Please, not that yp is not divided by tau_p. Instead, Hhat is divided by tau_p.
                % Thus, one can note sqrt(p)*tau_p instead of sqrt(p)*sqrt(tau_p) in yp

                yp_bar = sqrt(p)*tau_p*sum(HmeanFase((l-1)*N+1:l*N,:,t==pilotIndex),3);

                % MF MODIFICATION
                %---------------------------------------------------------%
                % -------------- COVARIANCE MATRICES, USAGE ------------- %
                %---------------------------------------------------------%
                %Compute the matrix in (4.6) that is inverted in the MMSE estimator % in (4.5)
                PsiInv = (p*tau_p*sum(R(:,:,l,t==pilotIndex),4) + eyeN);

                if strcmp(covariance_matrix, 'imperfect')
                    % Do not use a closed-form expression
                    R = R_imp;
                    PsiInv = Q_imp(:,:,l,t);

                end
                %---------------------------------------------------------%

                %Go through all UEs that use pilot t
                for k = find(t==pilotIndex)'

                    %Compute the MMSE estimate
                    RPsi = R(:,:,l,k) / PsiInv;
                    Hhat((l-1)*N+1:l*N,:,k) = HmeanFase((l-1)*N+1:l*N,:,k) + sqrt(p)*RPsi*(yp - yp_bar);

                end
            end
        end


        % =============================================================== %
        % ====================== LINEAR MMSE ============================ %
        % =============================================================== %
        % Perform LMMSE channel estimation, phase completely unknown
    case 'LMMSE'

        %Go through all APs
        for l = 1:L

            %Go through all pilots
            for t = 1:tau_p

                %Compute processed pilot signal for all UEs that use pilot t
                %according to (4.4) with an additional scaling factor \sqrt{tau_p}
                yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);
                % Please, not that yp is not divided by tau_p. Instead, Hhat is divided by tau_p.
                % Thus, one can note sqrt(p)*tau_p instead of sqrt(p)*sqrt(tau_p) in yp

                HmeanPower = zeros(N,N,K);
                for j=1:K
                    HmeanPower(:,:,j) = Hmean((l-1)*N+1:l*N,j)*Hmean((l-1)*N+1:l*N,j)';
                end

                % Marx M.M. Freitas MODIFICATION
                %---------------------------------------------------------%
                % -------------- COVARIANCE MATRICES, USAGE ------------- %
                %---------------------------------------------------------%
                % Compute the matrix in (4.6) that is inverted in the MMSE estimator in (4.5)
                PsiInv = (p*tau_p*(sum(R(:,:,l,t==pilotIndex),4) + sum(HmeanPower(:,:,t==pilotIndex),3)) + eyeN);

                if strcmp(covariance_matrix, 'imperfect')
                    % Do not use a closed-form expression
                    PsiInv = Q_imp(:,:,l,t);
                end
                %---------------------------------------------------------%


                % Go through all UEs that use pilot t
                for k = find(t==pilotIndex)'

                    % Compute the MMSE estimate FOR PERFECT R
                    RPsi = (R(:,:,l,k) + HmeanPower(:,:,k))/ PsiInv;

                    if strcmp(covariance_matrix, 'imperfect')
                        % Compute the MMSE estimate FOR IMPERFECT R
                        R = R_imp;
                        RPsi = R(:,:,l,k)/ PsiInv; % OK

                    end

                    Hhat((l-1)*N+1:l*N,:,k) = sqrt(p)*RPsi*yp;

                end
            end
        end

        % =============================================================== %
        % ========================== LEAST SQUARED ====================== %
        % =============================================================== %
    case 'LS'

        %Go through all APs
        for l = 1:L

            %Go through all pilots
            for t = 1:tau_p

                %Compute processed pilot signal for all UEs that use pilot t
                %according to (4.4) with an additional scaling factor \sqrt{tau_p}
                yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);


                %Go through all UEs that use pilot t
                for k = find(t==pilotIndex)'

                    %Compute the LS estimate
                    Hhat((l-1)*N+1:l*N,:,k) = yp/(sqrt(p)*tau_p);
                end

            end

        end

end

%% Covariance matrix of estimation error
for l=1:L
    for k=1:K
        for n=1:nbrOfRealizations
            Htil=H((l-1)*N+1:l*N,n,k)-Hhat((l-1)*N+1:l*N,n,k);
            C(:,:,l,k)=C(:,:,l,k) + (Htil*Htil')/nbrOfRealizations;
        end
    end
end

end

function Q_imp = Psi_imp_estimation(N, nbrOfAPs, tau_p, nbrOfRealizations, p, H, Np, pilotIndex)

% -------------------------------------------------------------------------
% COVARIANCE ESTIMATION "Q" (TOTAL)
% -------------------------------------------------------------------------
nbrOfRealizations_yp_n_R = nbrOfRealizations;
yp_n_R = zeros(N, nbrOfAPs, tau_p, nbrOfRealizations_yp_n_R);
yp_n_TIMES_yp_n_H_R = zeros(N, N, nbrOfAPs, tau_p, nbrOfRealizations_yp_n_R);
sum_yp_n_TIMES_yp_n_H_R = zeros(N, N, nbrOfAPs, tau_p);

for realization = 1:nbrOfRealizations_yp_n_R
    for indexAP = 1:nbrOfAPs
        for t = 1:tau_p % Go through all pilots
            yp_n_R(:, indexAP, t, realization) = sqrt(p)*sqrt(tau_p)*sum(H((indexAP-1)*N+1:indexAP*N, realization, t==pilotIndex), 3) + Np(:, realization, indexAP, t);
            yp_n_TIMES_yp_n_H_R(:, :, indexAP, t, realization) = yp_n_R(:,indexAP, t, realization)*transpose(conj(yp_n_R(:, indexAP, t, realization)));
            sum_yp_n_TIMES_yp_n_H_R(:, :, indexAP, t) = sum_yp_n_TIMES_yp_n_H_R(: , :, indexAP, t) + yp_n_TIMES_yp_n_H_R(:,:, indexAP, t, realization)/nbrOfRealizations_yp_n_R;
        end
    end
end

Q_imp = sum_yp_n_TIMES_yp_n_H_R;

end