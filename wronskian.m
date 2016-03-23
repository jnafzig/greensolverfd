function [ W] = wronskian(solution)
    %SOLUTION 
%     solution = [phiL,dphiL,phiR,dphiR];

    phiL = solution(:,1);
    dphiL = solution(:,2);
    phiR = solution(:,3);
    dphiR = solution(:,4);
    
    W = phiL.*dphiR - dphiL.*phiR;
    
end

 