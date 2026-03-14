function D = MF_FineTunning (fineTunning, D, R, gainOverNoise_dB, noiseVariance_dBm, PowerDL_mW, U_max, nbrOfAPs, nbrOfUEs, masterAPs,tau_p,tau_c,N,BW_Hz,ClusterType,EE_factor)
if strcmp (fineTunning, 'centSEbased')
    epsilon = 0.02;
    D = MF_fineTuning_SE(D,R,noiseVariance_dBm,PowerDL_mW,U_max,nbrOfAPs,nbrOfUEs,epsilon,tau_p,tau_c, masterAPs);

elseif strcmp (fineTunning, 'distPwBased')
    D = MF_fineTuning_Power(D,PowerDL_mW,nbrOfAPs,nbrOfUEs,gainOverNoise_dB,masterAPs);

% elseif strcmp (fineTunning, 'centEEbased')
%     epsilon = 10e-3; Gamma_prime = 0.9;
%     D = MF_fineTuning_EE_2(D,R,noiseVariance_dBm,PowerDL_mW,U_max,nbrOfAPs,nbrOfUEs,epsilon,tau_p,tau_c, masterAPs,Gamma_prime,N,BW_Hz);

elseif strcmp (fineTunning, 'centEEbased')
    epsilon = 10e-3; Gamma_prime = 0.9;
    D = MF_fineTuning_EE_2(D,R,noiseVariance_dBm,PowerDL_mW,U_max,nbrOfAPs,nbrOfUEs,epsilon,tau_p,tau_c, masterAPs,...
        Gamma_prime,N,BW_Hz,gainOverNoise_dB,ClusterType,EE_factor);

else
    
end