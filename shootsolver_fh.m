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
        f = ones(ceil(N),1);
        f(end) = 1+mod(N,-1);
        
        Emin = min(v);
        Ediv = Emin;
        Evecs = zeros(Nelem,ceil(N));
        Evals = zeros(ceil(N),1);
        n = zeros(Nelem,1);
        response = zeros(Nelem);
        for i = 1:ceil(N)

            nd = @(E) nodecount(shoot(E,v,vL,vR))-i;
            Ediv = fzero(nd,[Ediv,min(vL,vR)]);
 
            wron = @(E) mean(wronskian(shoot(E,v,vL,vR)));
            
            E = fzero(wron,[Emin,Ediv]);
          
            Emin = Ediv;
            
            [solution] = shoot(E,v,vL,vR);
%              
            phi = solution{1}(:,1);
            dphi = solution{2}(:,1);

%             phi = solution{1}(:,1);
%             dphi = solution{2}(:,1);            
%             dphi = dphi/(sum(phi.^2)^(1/2)*dx);
%             phi = phi/(sum(phi.^2)^(1/2)*dx);
            
            C = normfactor([phi;dphi;E]);
            X = [[phi;dphi]/C^(1/2);E];
            
            [C,dCdX] = normfactor(X);
            
%             dmatBC = sparse(bck_i,bck_j,...
%                 [dx-1i*(2 - 2*dx^2*(X(end)-vL))./sqrt((X(end)-vL)*2 -dx^2*(X(end)-vL)^2),...
%                 -dx+1i*(2 - 2*dx^2*(X(end)-vR))./sqrt((X(end)-vR)*2 -dx^2*(X(end)-vR)^2)],...
%                 2*Nelem+1,2*Nelem+1,2);

            [mat,dmatdE] = A(X(end));
            lhs = [[mat,dmatdE*X(1:end-1)];dCdX];
                
%             res = [A(X(end))*X(1:end-1);0]; 
%             X = X - lhs\res;
%             
%             dmatBC = sparse(bck_i,bck_j,...
%                 [dx-1i*(2 - 2*dx^2*(X(end)-vL))./sqrt((X(end)-vL)*2 -dx^2*(X(end)-vL)^2),...
%                 -dx+1i*(2 - 2*dx^2*(X(end)-vR))./sqrt((X(end)-vR)*2 -dx^2*(X(end)-vR)^2)],...
%                 2*Nelem+1,2*Nelem+1,2);
% 
%             lhs = [[A(X(end)),(A1+dmatBC)*X(1:end-1)];[X(1:Nelem)',zeros(1,Nelem+2)]];

            rhs = sparse(ival,jval,2*X(1:Nelem),2*Nelem+2,Nelem,Nelem);

            dXdv = lhs\rhs; 
%             
            ni = X(1:Nelem).^2/C;
%             dndX = diag(2*X(1:Nelem));
%             dndC = -X(1:Nelem).^2/C^2;
            
    
%             numdndX = numfuncderiv(@(density,X));
            n = n + f(i)*ni;
            response = response + f(i)*2*bsxfun(@times,X(1:Nelem),dXdv(1:Nelem,:));
        end
        
        function [mat,dmatdE] = A(E)
            sqL = 1i*sqrt((E-vL)*2-dx^2*(E-vL)^2);
            sqR = 1i*sqrt((E-vR)*2-dx^2*(E-vR)^2);
            matBC = sparse(bck_i,bck_j,...
                [-dx*(E-vL)+sqL,...
                  dx*(E-vR)-sqR],...
                2*Nelem+1,2*Nelem+1,2);

            dmatBC = sparse(bck_i,bck_j,...
                [-dx-(1 - dx^2*(E-vL))./sqL,...
                dx+(1 - dx^2*(E-vR))./sqR],...
                2*Nelem+1,2*Nelem+1,2);

            
            A0b = sparse(ival,jval,-2*v,2*Nelem+1,2*Nelem+1,Nelem);

            mat = A0a+A0b+matBC+A1*E;

            dmatdE = A1 + dmatBC;
        end
        
        function [int,dintdX] = normfactor(X)

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
            dintdphi = 2*dx*transpose(X(1:Nelem));
            dintdphi(1) = dintdphi(1) + dintLdX1;
            dintdphi(end) = dintdphi(end) + dintRdXend;
            
            dintdX = [dintdphi,zeros(1,Nelem+1),dintdE];
            

        end
        
    end
        
end

