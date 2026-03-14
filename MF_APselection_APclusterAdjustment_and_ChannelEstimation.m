%=========================================================================%
% =========================== AP SELECTION ============================== %
%=========================================================================%
[D, pilotIndex, masterAPs] = MF_APselectionSchemes(APselectionMethod, ...
    nbrOfAPs, nbrOfUEs, tau_p, U_max, gainOverNoise_dB);

%=========================================================================%
% ================= AP CLUSTER ADJUSTMENT (FINE TUNING) ================= %
%=========================================================================%
% Fine-tuning based on power allocation
if strcmp(fineTunning,'Algorithm_3') 
    Gamma_kl = 0.98; % Γ% of the cumulative sum in (24)
    D = MF_fineTuning_Power(D, PowerDL_mW, nbrOfAPs, nbrOfUEs, ...
        gainOverNoise_dB, masterAPs, Gamma_kl);

% Fine-tuning based on spectral efficiency (SE)
elseif strcmp(fineTunning,'Algorithm_4')
    varepsilon = 0.02; % threshold for SE losses in (26)
    D = MF_fineTuning_SE(D, R, noiseVariance_dBm, PowerDL_mW, U_max, ...
        nbrOfAPs, nbrOfUEs, varepsilon, tau_p, tau_c, masterAPs);

% Fine-tuning based on energy efficiency (EE)
elseif strcmp(fineTunning,'Algorithm_5')
    zeta = 10; % Ranged from 0 to 10 in the paper
    zeta = zeta/2;
    D = MF_fineTuning_EE(D, R, noiseVariance_dBm, PowerDL_mW, nbrOfAPs, ...
        nbrOfUEs, tau_p, tau_c, masterAPs, gainOverNoise_dB, zeta);

else
    % No_FineTuning
end

%=========================================================================%
% ========================== CHANNEL ESTIMATION ========================= %
%=========================================================================%
% Set channel estimation method
channelestimationMethod = 'LMMSE';% 'phaseAwareMMSE';'LMMSE';'LS';

% Generate channel realizations with estimates and estimation error correlation matrices
[Hhat,C] = functionChannelEstimatesWithRician(R, Hmean ,HmeanFase, H, ...
    nbrOfRealizations, nbrOfAPs, nbrOfUEs, N, tau_p, pilotIndex, ...
    poweUL_mW, channelestimationMethod, covariance_matrix, R_imp);