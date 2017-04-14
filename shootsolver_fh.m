function [ solver_fh ] = shootsolver_fh(Nelem, dx)
    %SOLVER provides solver for boundstates via a shooting method
    %
    % Returns a handle to a function which can provide a requested numer of
    % bound state eigenvalues
    %

    MID = spdiags(ones([Nelem,1]), 1,Nelem-1,Nelem+1);

    % D1 is first order finite difference
    loc = -1/2:1/2;
    coeff = fd_coeff(loc,1,dx);
    D1 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem-1,Nelem));       
    D1dphi = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1)); 

    % Boundary conditions: 
    % dphi(1) = (1-exp(1i*kL*dx))/dx*phi(1)
    % and dphi(end) = -(1-exp(1i*kR*dx))/dx*phi(end)
    
    % k independent portion of boundary condition
    bc = [[zeros(1,Nelem),1,zeros(1,Nelem)];...
          [zeros(1,Nelem),zeros(1,Nelem),1]];

    % k dependent portion of boundry condition
    bck = [[1,zeros(1,Nelem-1),zeros(1,Nelem+1)];...
            [zeros(1,Nelem-1),1,zeros(1,Nelem+1)]];
    [bck_i,bck_j] = find(bck);
        
    % k, v and E independent portion of lhs operator
    A0a = [bc;          ...
         [D1,-MID];     ...
         [zeros(Nelem),D1dphi]];
    
    % RHS operator
    A1 =  spdiags(2*ones(Nelem,1),-Nelem-1,2*Nelem+1,2*Nelem+1);
    [ival,jval] = find(A1);
    
    shoot = shoot_fh(Nelem,dx);
    
    solver_fh = @shootsolve;

    function [n,response] = shootsolve(N,v,vL,vR)
        if nargin < 3
            vL = 0;
            vR = 0;
        end
        % calculate scale factors
        if (N~=0)
            f = ones(ceil(N),1);
            f(end) = 1+mod(N,-1);
        else
            f = 0;
        end
        
        n = zeros(Nelem,1);
        response = zeros(Nelem);
        
        maxN = nodecount(shoot(min(vL,vR),v,vL,vR));
        if N>maxN
            warning('Max N is %i for this case',maxN);
            return;
        end
            
        
        Emin = min(v);
        for i = 1:ceil(N)
            
            % use bisection to find a quick upperbound for state with i nodes
            nd = @(E) nodecount(shoot(E,v,vL,vR))-i;
            Emax = fzero(nd,[Emin,min(vL,vR)]);
 
            % find the boundstate between Emin and Emax 
            wron = @(E) mean(wronskian(shoot(E,v,vL,vR)));
            [E,wronval] = fzero(wron,[Emin,Emax]);
          
            % updated Emin
            Emin = Emax;
            
            % get corresponding eigenvector
            [solution] = shoot(E,v,vL,vR);
            phi = solution{1}(:,1);
            dphi = solution{2}(:,1);
            
            X = [phi;dphi;E];
            
            % if wronval is not close enough to zero than do a step of
            % inverse iteration
            if abs(wronval) > 1e-8
                % get derivative of normfactor with respect to X
                [~,dCdX] = normfactor(X);

                % construct our matrix and the derivative of that matrix with
                % respect to E
                [mat,dmatdE] = A(X(end));

                lhs = [[mat,dmatdE*X(1:end-1)];dCdX];
                rhs = [mat*X(1:end-1);0];
                dX = lhs\rhs;
                
                X = X - dX;
            end
            
            % normalize.  (normalization depends on energy)
            C = normfactor(X);
            X = [X(1:end-1)/C^(1/2);X(end)];
            
            % add orbital contributions weighted by scale factor
            ni = X(1:Nelem).^2;
            n = n + f(i)*ni;
            
            if nargout > 1
                % get derivative of normfactor with respect to X
                [~,dCdX] = normfactor(X);

                % construct our matrix and the derivative of that matrix with
                % respect to E
                [mat,dmatdE] = A(X(end));

                lhs = [[mat,dmatdE*X(1:end-1)];dCdX];
                rhs = sparse(ival,jval,2*X(1:Nelem),2*Nelem+2,Nelem,Nelem);

                dXdv = lhs\rhs;              

                response = response + f(i)*2*bsxfun(@times,X(1:Nelem),dXdv(1:Nelem,:));
            end
        end
        
        function [mat,dmatdE] = A(E)

            sqL = 1i*sqrt((E-vL)*2-dx^2*(E-vL)^2);
            sqR = 1i*sqrt((E-vR)*2-dx^2*(E-vR)^2);
            matBC = sparse(bck_i,bck_j,...
                [-dx*(E-vL)+sqL,...
                  dx*(E-vR)-sqR],...
                2*Nelem+1,2*Nelem+1,2);

            A0b = sparse(ival,jval,-2*v,2*Nelem+1,2*Nelem+1,Nelem);
            mat = A0a+A0b+matBC+A1*E;

            if nargout == 1; return; end;
            
            dmatBC = sparse(bck_i,bck_j,...
                [-dx-(1 - dx^2*(E-vL))./sqL,...
                dx+(1 - dx^2*(E-vR))./sqR],...
                2*Nelem+1,2*Nelem+1,2);

            dmatdE = A1 + dmatBC;
        end
        
        function [int,dintdX] = normfactor(X)

            E = X(end);

            kR = acos(1-(E-vR)*dx^2)/dx; 
            dkRdE = dx/sqrt(1-(1-dx^2*(E-vR))^2);
            if imag(kR)<0; 
                kR = -kR;
                dkRdE = -dkRdE;
            end
            r = exp(1i*kR*dx);
            drdkR = 1i*dx*r;

            intR = dx*X(Nelem)^2*r^2/(1-r^2);
            dintRdr = dx*X(Nelem)^2*(2*r^3/(1-r^2)^2+2*r/(1-r^2));
            dintRdXnelem = 2*dx*X(Nelem)*r^2/(1-r^2);

            kL = acos(1-(E-vL)*dx^2)/dx;
            dkLdE = dx/sqrt(1-(1-dx^2*(E-vL))^2);
            if imag(kL)<0;
                kL = -kL;
                dkLdE = -dkLdE;
            end 
            l = exp(1i*kL*dx);
            dldkL = 1i*dx*l;

            intL = dx*X(1)^2*l^2/(1-l^2);
            dintLdl = dx*X(1)^2*(2*l^3/(1-l^2)^2+2*l/(1-l^2));
            dintLdX1 = 2*dx*X(1)*l^2/(1-l^2);

            int = (sum(X(1:Nelem).^2)*dx + intL + intR);
            dintdE = (dintRdr*drdkR*dkRdE + dintLdl*dldkL*dkLdE);
            dintdphi = 2*dx*transpose(X(1:Nelem));
            dintdphi(1) = dintdphi(1) + dintLdX1;
            dintdphi(end) = dintdphi(end) + dintRdXnelem;
            
            dintdX = [dintdphi,zeros(1,Nelem+1),dintdE];
            
        end
        
    end
        
end

