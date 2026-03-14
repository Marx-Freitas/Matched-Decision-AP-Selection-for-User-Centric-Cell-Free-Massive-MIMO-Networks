
%Determine the master AP for UE k by looking for AP with best
%channel condition
[~,master] = max(gainOverNoisedB(:,indexUE,nSetup));
D(master,indexUE,nSetup) = 1;
mastersOfUEs(indexUE,1) = master;

%Assign orthogonal pilots to the first tau_p UEs
if indexUE <= tau_p
    
    pilotIndex(indexUE,nSetup) = indexUE;
    
else %Assign pilot for remaining UEs
    
    %Compute received power from to the master AP from each pilot
    pilotinterference = zeros(tau_p,1);
    
    for t = 1:tau_p
        
        pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:indexUE-1,nSetup)==t,nSetup)));
        
    end
    
    %Find the pilot with the least receiver power
    [~,bestpilot] = min(pilotinterference);
    pilotIndex(indexUE,nSetup) = bestpilot;
    
end