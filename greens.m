function [ G] = greens(solution)
    
    phiL = solution(:,1);
    dphiL = solution(:,2);
    phiR = solution(:,3);
    dphiR = solution(:,4);
    
    W = phiL.*dphiR - dphiL.*phiR;
    
    
    A = triu((phiL./W)*transpose(phiR));
    G = (-1i/pi)*(A + transpose(A) - diag(diag(A)));
    
end
