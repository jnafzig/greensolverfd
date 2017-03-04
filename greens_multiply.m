function [ prod] = greens_multiply(solution,b)
    
    phiL = solution{1}(:,1);
    dphiL = solution{2}(:,1);
    phiR = solution{1}(:,2);
    dphiR = solution{2}(:,2);
    dphiL = (dphiL(1:end-1)+dphiL(2:end))/2;
    dphiR = (dphiR(1:end-1)+dphiR(2:end))/2;
    
    W = phiL.*dphiR - dphiL.*phiR;
    
    csumL = cumsum(bsxfun(@times,phiL,b));
    csumR = flipud(cumsum(flipud(bsxfun(@times,phiR,b))));
    diagLR = bsxfun(@times,(phiL.*phiR),b);

    prod = (-1i/pi)*bsxfun(@rdivide,(bsxfun(@times,phiR,csumL) ...
        + bsxfun(@times,phiL,csumR)...
        - diagLR),W);    

%     prod = diag(1./W)*(-1i/pi)*(diag(phiR)*cumsum(diag(phiL)*b) ...
%         + diag(phiL)*flipud(cumsum(flipud(diag(phiR)*b)))...
%         - diag(phiL.*phiR)*b);
%     A = triu((phiL./W)*transpose(phiR));
%     G = (-1i/pi)*(A + transpose(A) - diag(diag(A)));
    
end
