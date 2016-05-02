function [ solver_fh ] = eigsolver(Nelem, dx)
    %SOLVER provides solver for boundstates via a quadratic eigenvalue
    %solver recast as a linear eigenvalue problem.
    
    loc = -1/2:1/2;

    coeff = fd_coeff(loc,0,dx);
    D0 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1));
    MID = spdiags(ones([Nelem,1]), 1,Nelem-1,Nelem+1);

    coeff = fd_coeff(loc,1,dx);
    D1 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem-1,Nelem));       
    D1dphi = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1)); 

    bc0 = [zeros(2,Nelem),[[1,zeros(1,Nelem)];[zeros(1,Nelem),1]]];
    bc1 = [[[1,zeros(1,Nelem-1)];[zeros(1,Nelem-1),-1]],zeros(2,Nelem+1)];

    A0 = [bc0;[D1,-MID];[zeros(Nelem),D1dphi]];
    A1 = [bc1;zeros(2*Nelem-1,2*Nelem+1)];
    A2 =  spdiags(ones(Nelem,1),-Nelem-1,2*Nelem+1,2*Nelem+1);
    B = [[A1,A2];[speye(2*Nelem+1),zeros(2*Nelem+1)]];
    [ival,jval] = find(A2);

    solver_fh = @(N,v) eigsolve(N,v);

    function [Evals] = eigsolve(N,v)
        matV = sparse(ival,jval,-2*v,2*Nelem+1,2*Nelem+1,Nelem);

        A = blkdiag(A0+matV,speye(2*Nelem+1));
        
        kvals = eigs(A,B,N,sqrt(abs(2*min(v))));
        Evals = -kvals.^2/2;
    end
    
end
