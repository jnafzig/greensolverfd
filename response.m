function [ chi] = response(solution)
    
    phiL = solution(:,1);
    dphiL = solution(:,2);
    phiR = solution(:,3);
    dphiR = solution(:,4);
    
    W = phiL.*dphiR - dphiL.*phiR;
    
%     A = triu((phiL.^2./W.^2)*transpose(phiR.^2));
%     chi2 = -2i*(A + transpose(A) - diag(diag(A)))/pi;
    chi = resp_mex((-2i/pi)*(phiL./W).^2,complex(real(phiR.^2),imag(phiR.^2)));
end
