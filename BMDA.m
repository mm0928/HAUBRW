function [Rt]=BMDA(Wrr,Wdd,A,alpha)
% L1=2;L2=1;
% L1=2;
L2=1;
L1=1;
% alpha=0.5;
%bi-random walk on the heterogeneous network to calculate association scores for each miRNA-disease pair.

%Wrr: adjacency matrix of the miRNA similarity network
%Wdd: adjacency matrix of the disease similarity network
%A: adjacency matrix of the miRNA-disease association network
normWrr = normFun1(Wrr);
normWdd = normFun1(Wdd);
R0 = main_HeatS_ProbS1(A,0.9);
R0=R0/sum(R0(:));
Rt=R0;

% bi-random walk on the heterogeneous network
for t=1:max(L1,L2)
    
    ftl = 0;
    ftr = 0;
    
    %random walk on the miRNA similarity network
    if(t<=L1)
        nRtleft =(1- alpha) * normWrr * Rt + alpha*R0;
        ftl = 1;
    end
    %random walk on the disease similarity network
    if(t<=L2)
        nRtright = (1-alpha) *  Rt * normWdd + alpha*R0;
        ftr = 1;
    end
    
    %Rt: predictive association scores between each miRNA-disease pair
    Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
    
    
    
end


end

