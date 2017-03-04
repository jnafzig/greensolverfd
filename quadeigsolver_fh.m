function [ solver_fh ] = quadeigsolver_fh(Nelem, dx)
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
    
    % We will formulate the problem as follows:
    %
    %  A2*k^2 + A1*k + A0 = 0
    %
    % and then solve the linearization:
    %
    %      /        \   / \         /       \   / \
    %     | A0   A1  | | x |       | 0   -A2 | | x |
    %     |          | |   | = k * |         | |   |
    %     | 0     I  | | y |       | I    0  | | y |
    %      \        /   \ /         \       /   \ /
    %
    % the y has been constructed so that y = k*x 
    %
    % lhs operator will be A and rhs operator will be B
    
    loc = -1/2:1/2;    
    MID = spdiags(ones([Nelem,1]), 1,Nelem-1,Nelem+1);

    % D1 is first order finite difference
    coeff = fd_coeff(loc,1,dx);
    D1 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem-1,Nelem));       
    D1dphi = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1)); 

    % Linearized Boundary conditions: 
    % dphi(1) = k*phi(1) and dphi(end) = -k*phi(end)
    bc_dphi = [zeros(2,Nelem),[[1,zeros(1,Nelem)];...
                             [zeros(1,Nelem),1]]];
    bc_phi = [[[1,zeros(1,Nelem-1)];...
              [zeros(1,Nelem-1),-1]],zeros(2,Nelem+1)];

  % k dependent portion of boundry condition
    bck = [[1,zeros(1,Nelem-1),zeros(1,Nelem+1)];...
            [zeros(1,Nelem-1),1,zeros(1,Nelem+1)]];
    [bck_i,bck_j] = find(bck);

          
    % the part of A0 that doesn't depend on the potential:
    A0a = [bc_dphi;          ...
         [D1,-MID];     ...
         [sparse(Nelem,Nelem),D1dphi]];
    
    % RHS operator
    A1 = [bc_phi;sparse(2*Nelem-1,2*Nelem+1)];
    A2 =  spdiags(ones(Nelem,1),-Nelem-1,2*Nelem+1,2*Nelem+1);
    B = [[A1             ,A2              ];...
        [speye(2*Nelem+1),sparse(2*Nelem+1,2*Nelem+1)]];
    
    [ival,jval] = find(A2);

    % return function handle
    solver_fh = @eigsolve;

    function [Evals,Evecs] = eigsolve(N,v,~,~)
        
        % the part of A0 that does depend on the potential
        A0b = sparse(ival,jval,-2*v,2*Nelem+1,2*Nelem+1,Nelem);

        A = blkdiag(A0a+A0b,speye(2*Nelem+1));
        
        [Evecs,kvals] = eigs(A,B,ceil(N),sqrt(abs(2*min(v))));
        Evecs = Evecs(1:Nelem,:);
        Evals = -diag(kvals).^2/2;
        vL = 0;
        vR = 0;
        
        function mat = T(E)
%             kL = acos(1-(E-vL)*dx^2)/dx;
%             if imag(kL)<0; kL = -kL; end  % ensure we choose the sqrt with 
%                                           % negative imaginary part
%             kR = acos(1-(E-vR)*dx^2)/dx;
%             if imag(kR)<0; kR = -kR; end
            
            matBC = sparse(bck_i,bck_j,...
                [dx*(E-vL)-1i*sqrt((E-vL)*2-dx^2*(E-vL)^2),...
                -dx*(E-vR)+1i*sqrt((E-vR)*2-dx^2*(E-vR)^2)],...
                2*Nelem+1,2*Nelem+1,2);

            A0b = sparse(ival,jval,-2*(v-E),2*Nelem+1,2*Nelem+1,Nelem);

            mat = A0a+A0b+matBC;

        end
    end
        
end

