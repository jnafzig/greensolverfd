function [ solver_fh ] = contoureigsolver_fh(Nelem, dx, RelTol,AbsTol)
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

    
    shoot = shoot_fh(Nelem,dx);
    
    % return function handle
    solver_fh = @eigsolve;

    function [Evals,Evecs] = eigsolve(N,v,vL,vR)
        if nargin < 3
            vL = 0;
            vR = 0;
        end
        b = rand(Nelem,N)-1/2;
        minE = min(v)-1;
        maxE = min(vL,vR)-.25;

        % Define the contour
        E0 = (maxE+minE)/2; 
        R = (maxE-minE)/2; 
        E = @(theta) E0 + R*exp(1i*theta);
        dEdt = @(theta) 1i*R*exp(1i*theta); 

        greenmult = @(E) greens_multiply(shoot(E,v,vL,vR),b);

        T0 = integral(@(theta) ((greenmult(E(theta)))*dEdt(theta)),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
        T0 = T0 + conj(T0);
        
        T1 = integral(@(theta) E(theta)*(greenmult(E(theta)))*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
        T1 = T1 + conj(T1);
        
        [V,Sig,W] = svd(T0);
        B = (V(:,1:N)'*T1*W)/Sig(1:N,1:N);

        [s,Evals] = eig(B);

        Evals = sort(diag(Evals));
        Evecs = V(:,1:N)*s;

    end
    
    
end
