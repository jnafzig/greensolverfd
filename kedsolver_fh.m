function [ kedsolver_fh ] = kedsolver_fh(Nelem, dx, Tol)
    % provides a handle to a function which can return density and kinetic
    % energy density.  
    if nargin<3
        RelTol = eps;
        AbsTol = eps;
    else
        RelTol = Tol;
        AbsTol = Tol;
    end
    
    shoot = shoot_fh(Nelem,dx);
    dens = @(E,v,vL,vR) ldos(shoot(E,v,vL,vR));
    
    % return function handle
    kedsolver_fh = @kedsolve;

    function [n,ts] = kedsolve(mu,v,vL,vR)
        
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
        
        e = integral(@(theta) E(theta)*dens(E(theta),v,vL,vR)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
        e = e+conj(e);
        
        ts = e - v.*n;
        
    end
        
end

