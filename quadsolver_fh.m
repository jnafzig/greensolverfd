function [ solver_fh ] = quadsolver_fh(Nelem, dx)
    %SOLVER provides solver for boundstates via a quadratic eigenvalue
    %solver recast as a linear eigenvalue problem.
    %
    % Returns a handle to a function which can provide density and response
    %

    % Schrodinger equation is eig problem is H*phi = E*phi,
    % We consider a vector x = [phi,dphi] then our first eq
    % requires that D1*phi - dphi = 0;
    
    % but boundary conditions depend on k = sqrt(E):
    % dphi(1) = k*phi(1) and dphi(end) = k*phi(end)
    %
    % So we have a quadratic eigenvalue problem. 
    % We convert to a linear eigenvalue problem:
    %

    loc = -1/2:1/2;
    
    MID = spdiags(ones([Nelem,1]), 1,Nelem-1,Nelem+1);

    % D1 is first order finite difference
    coeff = fd_coeff(loc,1,dx);
    D1 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem-1,Nelem));       
    D1dphi = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1)); 

    % Boundary conditions: 
    % dphi(1) = k*phi(1) and dphi(end) = -k*phi(end)
    bc_lhs = [sparse(2,Nelem),[[1,sparse(1,Nelem)];...
                             [sparse(1,Nelem),1]]];
    bc_rhs = [[[1,sparse(1,Nelem-1)];...
              [sparse(1,Nelem-1),-1]],sparse(2,Nelem+1)];

    % LHS operator without potential
    A0 = [bc_lhs;          ...
         [D1,-MID];     ...
         [sparse(Nelem,Nelem),D1dphi]];
    
    % RHS operator
    A1 = [bc_rhs;sparse(2*Nelem-1,2*Nelem+1)];
    A2 =  spdiags(ones(Nelem,1),-Nelem-1,2*Nelem+1,2*Nelem+1);
    B = [[A1             ,A2              ];...
        [speye(2*Nelem+1),sparse(2*Nelem+1,2*Nelem+1)]];
    
    [ival,jval] = find(A2);

    % return function handle
    solver_fh = @densitysolve;

    function [n,response] = densitysolve(N,v,~,~)
        if N == 0
            n = zeros(Nelem,1);
            response = zeros(Nelem);
            return;
        end
        Nbs = ceil(N); % number of bound states to calculate
        Nfs = floor(N); % number of completely full states
        
        matV = sparse(ival,jval,-2*v,2*Nelem+1,2*Nelem+1,Nelem);

        A = blkdiag(A0+matV,speye(2*Nelem+1));
        
        [Evecs,kvals] = eigs(A,B,Nbs,sqrt(abs(2*min(v))));

        % normalization
        int = (Evecs(1,:).^2+Evecs(Nelem,:).^2)/(2*kvals);
        C = sum(Evecs(1:Nelem,:).^2)*dx + int;
        Evecs = Evecs*diag(C.^-(1/2));
        
        Lvecs =  Evecs(1:Nelem,:);
        Lvecs(1,:) = Lvecs(1,:) + Lvecs(1,:)/(2*kvals*dx);
        Lvecs(Nelem,:) = Lvecs(Nelem,:) + Lvecs(Nelem,:)/(2*kvals*dx);
        
        n = sum(Evecs(1:Nelem,1:Nfs).^2,2) ...
            + (N-Nfs)*Evecs(1:Nelem,end).^2;
        
        if nargout>1
            response = zeros(Nelem);
            for i = 1:Nbs
                Evec = Evecs(:,i);
                Lvec = Lvecs(:,i);
                kval = kvals(i,i);
                lhs = [[A-kval*B,B*Evec];...
                    [transpose(Lvec(1:Nelem)),...
                    zeros(1,3*Nelem+2),...
                    -(Evec(1)^2+Evec(Nelem)^2)/kval^2]];
                rhs = sparse(ival,jval,2*Evec(1:Nelem),4*Nelem+3,Nelem,Nelem);
                dndvec1 = sparse(1:Nelem,1:Nelem,2*Evec(1:Nelem),Nelem,4*Nelem+3,Nelem);
                
                dvecdv = (lhs\rhs);
                
                if i<=Nfs
                    scale = 1;
                else
                    scale = N-Nfs;
                end
                response = response + scale*(dndvec1*dvecdv ...
                    - 2*dx*Evec(1:Nelem).^2*(Evec(1:Nelem)'*dvecdv(1:Nelem,:)) ...
                    - Evec(1)/kval*Evec(1:Nelem).^2*dvecdv(1,:) ...
                    - Evec(Nelem)/kval*Evec(1:Nelem).^2*dvecdv(Nelem,:) ...
                    - Evec(1:Nelem).^2*(Evec(1)^2+Evec(Nelem)^2)/kval^2/2*dvecdv(end,:));
            end
            
        end
    end
        
end

