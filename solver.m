function [ solver_fh ] = solver(Nelem, dx)
    %SOLVER 
    
    loc = -1/2:1/2;

    coeff = fd_coeff(loc,0,dx);
    D0 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1));
    MID = spdiags(ones([Nelem,1]), 1,Nelem-1,Nelem+1);
    
    coeff = fd_coeff(loc,1,dx);
    D1 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem-1,Nelem));       
    D1dphi = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1)); 
        
    bcL = [[1,zeros(1,2*Nelem)];[zeros(1,Nelem),1,zeros(1,Nelem)]];
    bcR = [[zeros(1,Nelem-1),1,zeros(1,Nelem+1)];[zeros(1,2*Nelem),1]];

    matL = [bcL;[D1,-MID];[zeros(Nelem),D1dphi]];
    matR = [bcR;[D1,-MID];[zeros(Nelem),D1dphi]];

    
    matVE =  spdiags(ones(Nelem,1),-Nelem-1,2*Nelem+1,2*Nelem+1);
       
    [ival,jval] = find(matVE);
    
    solver_fh = @(E,v,vL,vR) shoot(E,v,vL,vR);

    function [solution] = shoot(E,v,vL,vR)
        
        matVE = sparse(ival,jval,2*(E-v),2*Nelem+1,2*Nelem+1,Nelem);

        kL = sqrt(2*(E-vL));
        if imag(kL)<0; kL = -kL; end
        
        matchL = matL+matVE;
        rhsL = [1;-1i*kL;zeros(2*Nelem-1,1)];
        solL = matchL\rhsL;
        phiL = solL(1:Nelem);
        dphiL = D0*solL((1:Nelem+1)+Nelem);
        
        kR = sqrt(2*(E-vR)); 
        if imag(kR)<0; kR = -kR; end
        
        matchR = matR+matVE;
        rhsR = [1;1i*kR;zeros(2*Nelem-1,1)]; 
        solR = matchR\rhsR;
        phiR = solR(1:Nelem);
        dphiR = D0*solR((1:Nelem+1)+Nelem);
        
        solution = [phiL,dphiL,phiR,dphiR];

    end
end

