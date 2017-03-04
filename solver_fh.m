function [ solver_fh ] = solver_fh(Nelem, dx, Tol)
    if nargin<3
        RelTol = eps;
        AbsTol = eps;
    else
        RelTol = Tol;
        AbsTol = Tol;
    end
    
    shoot = shoot_fh(Nelem,dx);
    dens = @(E,v,vL,vR) ldos(shoot(E,v,vL,vR));
    resp = @(E,v,vL,vR) response(shoot(E,v,vL,vR));
    
    % return function handle
    solver_fh = @densitysolve;

    function [n,response] = densitysolve(mu,v,vL,vR)
        
        E0 = min(v)-1;

        R = (E0+mu)/2; 
        A = mu-R; 
        E = @(theta) R + A*exp(1i*theta);
        dEdt = @(theta) 1i*A*exp(1i*theta); 
 
        n = integral(@(theta) dens(E(theta),v,vL,vR)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
        n = n+conj(n);
        
        if nargout>1
            response = dx*integral(@(theta) resp(E(theta),v,vL,vR)*dEdt(theta),0,pi,...
                    'ArrayValued',true,...
                    'RelTol',1e-1,...
                    'AbsTol',1e-1);
            response = response+conj(response);
        end
    end
        
end

