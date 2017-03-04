function [ W] = wronskian(solution)
    %SOLUTION 
%     solution = [phiL,dphiL,phiR,dphiR];

    phiL = solution{1}(:,1);
    dphiL = solution{2}(:,1);
    phiR = solution{1}(:,2);
    dphiR = solution{2}(:,2);
    dphiL = (dphiL(1:end-1)+dphiL(2:end))/2;
    dphiR = (dphiR(1:end-1)+dphiR(2:end))/2;
    
    W = phiL.*dphiR - dphiL.*phiR;
    
end

 