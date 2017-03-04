function [ solver_fh ] = density_fh(Nelem, dx)
    
    solver_fh = @shootsolve;

    function [n,dndX] = shootsolve(X,vL,vR)

        E = X(end);
        
        kR = acos(1-(E-vR)*dx^2)/dx; 
        dkRdE = dx/sqrt(1-(1-dx^2*(E-vR))^2);
        if imag(kR)<0; kR = -kR; end
        r = exp(1i*kR*dx);
        drdkR = 1i*dx*r;

        intR = dx*X(Nelem)^2*r^2/(1-r^2);
        dintRdr = dx*X(Nelem)^2*(2*r^3/(1-r^2)^2+2*r/(1-r^2));
        dintRdXend = 2*dx*X(Nelem)*r^2/(1-r^2);

        kL = acos(1-(E-vL)*dx^2)/dx;
        dkLdE = dx/sqrt(1-(1-dx^2*(E-vL))^2);
        if imag(kL)<0; kL = -kL; end 
        r = exp(1i*kL*dx);
        drdkL = 1i*dx*r;

        intL = dx*X(1)^2*r^2/(1-r^2);
        dintLdr = dx*X(1)^2*(2*r^3/(1-r^2)^2+2*r/(1-r^2));
        dintLdX1 = 2*dx*X(1)*r^2/(1-r^2);

        int = (sum(X(1:Nelem).^2)*dx + intL + intR);
        dintdE = (dintRdr*drdkR*dkRdE + dintLdr*drdkL*dkLdE);
        dintdX = 2*dx*X(1:Nelem);
        
        n = X(1:Nelem).^2/int;
        dndint = -X(1:Nelem).^2/int^2;

        dndE = dndint*dintdE;
        dndX1 = dndint*dintLdX1;
        dndXend = dndint*dintRdXend;
        dndphi = dndint*dintdX';

        dndX = [sparse(1:Nelem,1:Nelem,...
            2*X(1:Nelem),Nelem,2*Nelem+1,Nelem),dndE];
        dndX(:,1) = dndX(:,1) + dndX1;
        dndX(:,Nelem) = dndX(:,Nelem) + dndXend;
        dndX(:,1:Nelem) = dndX(:,1:Nelem) + dndphi; 
        
    end
        
end
        

