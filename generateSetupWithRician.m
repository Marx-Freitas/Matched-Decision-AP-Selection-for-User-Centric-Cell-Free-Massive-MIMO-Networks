function [gainOverNoisedB,R,Hmean,HmeanFase,H,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetupWithRician(L,K,N,tau_p,nbrOfRealizations,nbrOfSetups,seed,ASD_varphi,ASD_theta)
%This function generates realizations of the simulation setup described in
%Section 5.3.
%
%INPUT:
%L               = Number of APs per setup
%K               = Number of UEs in the network
%N               = Number of antennas per AP
%tau_p           = Number of orthogonal pilots
%nbrOfSetups     = Number of setups with random UE and AP locations
%seed            = Seed number of pseudorandom number generator
%ASD_varphi      = Angular standard deviation in the local scattering model 
%                  for the azimuth angle (in radians)
%ASD_theta       = Angular standard deviation in the local scattering model
%                  for the elevation angle (in radians)
%
%OUTPUT:
%gainOverNoisedB = Matrix with dimension L x K x nbrOfSetups where
%                  element (l,k,n) is the channel gain (normalized by the
%                  noise variance) between AP l and UE k in setup n
%R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                  where (:,:,l,k,n) is the spatial correlation matrix
%                  between AP l and UE k in setup n, normalized by noise
%pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                  pilots assigned to the UEs
%D               = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for cell-free setup
%D_small         = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for small-cell setup
%APpositions     = Vector of length L with the AP locations, where the real
%                  part is the horizontal position and the imaginary part
%                  is the vertical position
%UEpositions     = Vector of length K with UE positions, measured in the
%                  same way as APpositions
%distances       = Matrix with same dimension as gainOverNoisedB containing
%                  the distances in meter between APs and UEs
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


%% Define simulation setup

%Set the seed number if it is specified other than zero
if (nargin>5)&&(seed>0)
    rng(seed)
end

%Size of the coverage area (as a square with wrap-around)
squareLength = 1000; %meter

%Communication bandwidth (Hz)
B = 20e6;

%Parameters for 3GPP path loss model
model = 'Umi'; % 'Umi' 'InH-Mixed' 'InH-Open' 'InF-SL' 'InF-DL'

fc_GHz = 3.6;

h_AP_m = 11.65;

h_UE_m = 1.65;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters for the model in (5.42)
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Height difference between an AP and a UE (in meters)
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Prepare to save results
gainOverNoisedB = zeros(L,K,nbrOfSetups);
R = zeros(N,N,L,K,nbrOfSetups);
distances = zeros(L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);
D = zeros(L,K,nbrOfSetups);
D_small = zeros(L,K,nbrOfSetups);

masterAPs = zeros(K,1); %the indices of master AP of each UE k 

LosVerifier = zeros(L,K);
PrLOS = zeros(L,K);
ricianFactor = zeros(L,K);

Hnlos = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));
Hmean = zeros(L*N,K);
HmeanFase = zeros(L*N,nbrOfRealizations,K);
H = zeros(L*N,nbrOfRealizations,K);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Random AP locations with uniform distribution
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;
    
    %Prepare to compute UE locations
    UEpositions = zeros(K,1);
    
    
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    
    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowAPrealizations = zeros(K,L);
    
    
    
    %Add UEs
    for k = 1:K
        
        %Generate a random UE location in the area
        UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
        
        %Compute distances assuming that the APs are 10 m above the UEs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEposition,size(APpositionsWrapped))),[],2);
        distances(:,k,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
        
        %If this is not the first UE
        if k-1>0
            
            %Compute distances from the new prospective UE to all other UEs
            shortestDistances = zeros(k-1,1);
            
            for i = 1:k-1
                shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
            end
            
            %Compute conditional mean and standard deviation necessary to
            %obtain the new shadow fading realizations, when the previous
            %UEs' shadow fading realization have already been generated.
            %This computation is based on Theorem 10.2 in "Fundamentals of
            %Statistical Signal Processing: Estimation Theory" by S. Kay
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
            meanvalues = term1*shadowAPrealizations(1:k-1,:);
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
            
        else %If this is the first UE
            
            %Add the UE and begin to store shadow fading correlation values
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];
            
        end
        
        %Generate the shadow fading realizations
        shadowing = meanvalues + stdvalue*randn(1,L);
        
        %Compute the channel gain divided by noise power
        [gainOverNoise_dB,PrLOS(:,k),LosVerifier(:,k)] = AF_LargeScaleGainOverNoise(model,fc_GHz,noiseVariancedBm,distances(:,k,n),h_AP_m,h_UE_m);
        gainOverNoisedB(:,k,n) = gainOverNoise_dB + shadowing';
        
        %Rician factor
        ricianFactor(:,k) = PrLOS(:,k)./(1-PrLOS(:,k));
        
        %Update shadowing correlation matrix and store realizations
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k,:) = shadowing;
        
        %Store the UE position
        UEpositions(k) = UEposition;
        
        
        %Determine the master AP for UE k by looking for AP with best
        %channel condition
        [~,master] = max(gainOverNoisedB(:,k,n));
        D(master,k,n) = 1;
        masterAPs(k) = master;
        
        %Assign orthogonal pilots to the first tau_p UEs according to
        %Algorithm 4.1
        if k <= tau_p
            
            pilotIndex(k,n) = k;
            
        else %Assign pilot for remaining UEs
            
            %Compute received power to the master AP from each pilot
            %according to Algorithm 4.1
            pilotinterference = zeros(tau_p,1);
            
            for t = 1:tau_p
                
                pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1,n)==t,n)));
                
            end
            
            %Find the pilot with the least receiver power according to
            %Algorithm 4.1
            [~,bestpilot] = min(pilotinterference);
            pilotIndex(k,n) = bestpilot;
            
        end
        
        
        
        %Go through all APs
        for l = 1:L
            
            %Compute nominal angle between UE k and AP l
            angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); %azimuth angle
            angletoUE_theta = asin(distanceVertical/distances(l,k,n));  %elevation angle
            %Generate spatial correlation matrix using the local
            %scattering model in (2.18) and Gaussian angular distribution
            %by scaling the normalized matrices with the channel gain
            if nargin>7
                R(:,:,l,k,n) = functionRlocalscatteringWithRician(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            else
                R(:,:,l,k,n) = eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
            end
            
            
            %% Generate Rician channel realizations
       
            if LosVerifier(l,k)==0
                
                H((l-1)*N+1:l*N,:,k) = sqrt(0.5*db2pow(gainOverNoisedB(l,k)))*sqrtm(R(:,:,l,k,n))*Hnlos((l-1)*N+1:l*N,:,k);
                
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k))*R(:,:,l,k,n);
            else
                if ricianFactor(l,k) == Inf
                    factorNLoS = 0;
                    factorLoS = 1;
                else
                    factorNLoS = 1/(ricianFactor(l,k)+1);
                    factorLoS = ricianFactor(l,k)/(ricianFactor(l,k)+1);
                end
                          
                %Apply correlation to the Rician channel realizations
                % LOS channel using 3D ULA response
                dn = linspace(1,N,N)';
                Hmean((l-1)*N+1:l*N,k) = (sqrt(factorLoS*db2pow(gainOverNoisedB(l,k)))*exp(-1j*pi*sin(angletoUE_varphi)*cos(angletoUE_theta).*(dn-1)));
                
                HmeanFase((l-1)*N+1:l*N,:,k) = Hmean((l-1)*N+1:l*N,k)*exp(1j*rand(1,nbrOfRealizations)*2*pi);                
                
                H((l-1)*N+1:l*N,:,k) = HmeanFase((l-1)*N+1:l*N,:,k) ...
                    + sqrt(0.5*factorNLoS*db2pow(gainOverNoisedB(l,k)))*sqrtm(R(:,:,l,k,n))*Hnlos((l-1)*N+1:l*N,:,k);
                
                R(:,:,l,k,n) = factorNLoS*db2pow(gainOverNoisedB(l,k))*R(:,:,l,k,n);
            end            
                           
        end
        
    end
    
    
    %Each AP serves the UE with the strongest channel condition on each of
    %the pilots in the cell-free setup
    for l = 1:L
        
        for t = 1:tau_p
            
            pilotUEs = find(t==pilotIndex(:,n));
            [~,UEindex] = max(gainOverNoisedB(l,pilotUEs,n));
            D(l,pilotUEs(UEindex),n) = 1;
           
        end
        
    end
    
    %Determine the AP serving each UE in the small-cell setup according to
    %(5.47) by considering only the APs from the set M_k for UE k, i.e.,
    %where D(:,k,n) is one.
    for k=1:K
        
        tempmat = -inf*ones(L,1);
        tempmat(D(:,k,n)==1,1) = gainOverNoisedB(D(:,k,n)==1,k,n);
        [~,servingAP] = max(tempmat);
        D_small(servingAP,k,n) = 1;
        
    end
    
    
end

