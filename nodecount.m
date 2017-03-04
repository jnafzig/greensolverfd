function [ nodes] = nodecount(solution)
    %SOLUTION 
%     solution = [phiL,dphiL,phiR,dphiR];

    phiL = solution{1}(:,1);
    dphiL = solution{2}(:,1);
    phiR = solution{1}(:,2);
    dphiR = solution{2}(:,2);
    
    
    nodes = sum(abs(diff(heaviside(phiL))))+...
        heaviside(dphiR(end)./phiR(end)-dphiL(end)./phiL(end));
     
end

 