function [ combined] = density_response(solution)
    
% Combining density and response integrals apparently doesn't speed things 
% up, but this is the neccesary function to do that.
    
    phiL = solution{1}(:,1);
    dphiL = solution{2}(:,1);
    phiR = solution{1}(:,2);
    dphiR = solution{2}(:,2);
    dphiL = (dphiL(1:end-1)+dphiL(2:end))/2;
    dphiR = (dphiR(1:end-1)+dphiR(2:end))/2;
    
    W = phiL.*dphiR - dphiL.*phiR;
    
    dn = 1/(pi*1i)*phiL.*phiR./W;
    dchi = resp_mex((-2i/pi)*(phiL./W).^2,complex(real(phiR.^2),imag(phiR.^2)));
    
    combined = [dn,dchi];
end
