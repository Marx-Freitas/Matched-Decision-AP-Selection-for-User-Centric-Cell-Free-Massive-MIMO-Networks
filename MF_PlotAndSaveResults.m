%% Plot simulation results
figure; hold on; box on;
plot(sort(SE_P_MMSE_tot(:)),linspace(0,1,nbrOfUEs*nbrOfSetups),'r-.','LineWidth',2);
plot(sort(SE_LP_MMSE_tot(:)),linspace(0,1,nbrOfUEs*nbrOfSetups),'b--','LineWidth',2);
plot(sort(SE_MR_tot(:)),linspace(0,1,nbrOfUEs*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_P_RZF_tot(:)),linspace(0,1,nbrOfUEs*nbrOfSetups),'y-','LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'P-MMSE (Scalable)','LP-MMSE (Scalable)','MR (Scalable)','P-RZF (Scalable)'},'Interpreter','Latex','Location','SouthEast');
xlim([0 14]);

clear R p nSetup Hhat H gainOverNoisedB C B D pilotIndex...
    UEpositions APpositions APPositionsXYZ_m APpositionsWrapped

save results.mat