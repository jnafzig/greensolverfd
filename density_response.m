function [ combined] = density_response(solution)
    
% Combining density and response integrals apparently doesn't speed things 
% up, but this is the neccesary function to do that.
    
    phiL = solution(:,1);
    dphiL = solution(:,2);
    phiR = solution(:,3);
    dphiR = solution(:,4);
    
    W = phiL.*dphiR - dphiL.*phiR;
    
    dn = 1/(pi*1i)*phiL.*phiR./W;
    dchi = resp_mex((-2i/pi)*(phiL./W).^2,complex(real(phiR.^2),imag(phiR.^2)));
    
    combined = [dn,dchi];
end
