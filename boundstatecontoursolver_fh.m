function [ solver_fh ] = boundstatecontoursolver_fh(Nelem, dx, RelTol,AbsTol)
    % Provides solver for boundstates via contour integral method 
    % described in: "An integral method for solving nonlinear eigenvalue 
    % problems" Wolf-Jürgen Beyn

    % Schrodinger equation is linear eig problem, but
    % boundary conditions are nonlinear
    
    % We consider a vector x = [phi,dphi] then our first matrix block
    % requires that 
    %
    % D1*phi - dphi = 0;
    %
    % and the second requires that
    %
    % -1/2*D1*dphi + (v-E)*phi = 0  
    
    % but boundary conditions depend on 
    % k ->  E-v = 1/dx^2 *(1-cos(k*dx)):
    % dphi(1) = 1/dx*(1-exp(-1i*k*dx))*phi(1) and 
    % dphi(end) = 1/dx*(exp(1i*k*dx)-1)*phi(end)
    %
    if nargin <3
        RelTol = eps;
        AbsTol = eps;
    end

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
    bc = [[zeros(1,Nelem),-1,zeros(1,Nelem)];...
          [zeros(1,Nelem),zeros(1,Nelem),-1]];

    % k dependent portion of boundry condition
    bck = [[1,zeros(1,Nelem-1),zeros(1,Nelem+1)];...
            [zeros(1,Nelem-1),1,zeros(1,Nelem+1)]];
    [bck_i,bck_j] = find(bck);
        
    % k, v and E independent portion of lhs operator
    A0 = [bc;          ...
         [D1,-MID];     ...
         [zeros(Nelem),D1dphi]];
    
    % RHS operator
    A2 =  spdiags(ones(Nelem,1),-Nelem-1,2*Nelem+1,2*Nelem+1);
    [ival,jval] = find(A2);

    % return function handle
    solver_fh = @eigsolve;

    function [Evals,Evecs] = eigsolve(N,v,vL,vR)

        b = rand(Nelem,N)-1/2;
        minE = min(v)-1;
        maxE = min(vL,vR)-.25;

        E0 = (maxE+minE)/2; 
        R = (maxE-minE)/2; 
        E = @(theta) E0 + R*exp(1i*theta);
        dEdt = @(theta) 1i*R*exp(1i*theta); 

        shoot = shoot_fh(Nelem,dx);
        greenmult = @(E) greens_multiply(shoot(E,v,vL,vR),b);

        T0 = integral(@(theta) ((greenmult(E(theta)))*dEdt(theta)),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
        T0 = T0 + conj(T0);
%         T0 = 2*real(T0);
        
        T1 = integral(@(theta) E(theta)*(greenmult(E(theta)))*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
        T1 = T1 + conj(T1);
        
        [V,Sig,W] = svd(T0);
        B = (V(:,1:N)'*T1*W)/Sig(1:N,1:N);

        [s,Evals] = eig(B);

        Evecs = V(:,1:N)*s;
%         
%         function mat = T(E)
%             
% %             kL = acos(1-(E-vL)*dx^2)/dx;
% %             if imag(kL)<0; kL = -kL; end  % ensure we choose the sqrt with 
% %                                           % negative imaginary part
% %             kR = acos(1-(E-vR)*dx^2)/dx;
% %             if imag(kR)<0; kR = -kR; end
%             
%             matBC = sparse(bck_i,bck_j,...
%                 [dx*(E-vL)-1i*sqrt((E-vL)*2-dx^2*(E-vL)^2),...
%                 -dx*(E-vR)+1i*sqrt((E-vR)*2-dx^2*(E-vR)^2)],...
%                 2*Nelem+1,2*Nelem+1,2);
% 
%             matV = sparse(ival,jval,-2*(v-E),2*Nelem+1,2*Nelem+1,Nelem);
% 
%             A = A0+matV+matBC;
%             
%             mat = A; 
% 
%         end
    end
    
    
end
