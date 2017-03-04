function [ chi] = response(solution)
    
    phiL = solution{1}(:,1);
    dphiL = solution{2}(:,1);
    phiR = solution{1}(:,2);
    dphiR = solution{2}(:,2);
    dphiL = (dphiL(1:end-1)+dphiL(2:end))/2;
    dphiR = (dphiR(1:end-1)+dphiR(2:end))/2;
    
    
    W = phiL.*dphiR - dphiL.*phiR;
    
%     A = triu((phiL.^2./W.^2)*transpose(phiR.^2));
%     chi2 = -2i*(A + transpose(A) - diag(diag(A)))/pi;
    chi = resp_mex((-2i/pi)*(phiL./W).^2,complex(real(phiR.^2),imag(phiR.^2)));
end
