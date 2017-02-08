function  v  = invert( solver, n0, mu, v0, vL, vR, TolFun )
    if nargin < 7
        TolFun = 1e-6;
    end
    
    %INVERT tries to find v such that v -> n0 at a chemical potential mu

    vsR = 0; % currently the inverted potential is assumed to be zero 
    vsL = 0; % outside the gridded region.  For now vsR and vsL can be 
             % adjusted manually to play around with that assumption.
    
%     RelTol = eps;
%     AbsTol = eps;
    
%     Nelem = numel(n0);
%     shoot = shoot_fh(Nelem,dx);

%     dens = @(E,v,vL,vR) ldos(shoot(E,v,vL,vR));
%     resp = @(E,v,vL,vR) response(shoot(E,v,vL,vR));

% Combining density and response integrals apparently doesn't speed things 
% up, but this is part of the code to do that: 
%     dens_resp = @(E,v,vL,vR) density_response(shoot(E,v,vL,vR));
    
%     E0 = min(v0)-1;
% 
%     R = (E0+mu)/2;
%     A = mu-R;
%     E = @(theta) R + A*exp(1i*theta);
%     dEdt = @(theta) 1i*A*exp(1i*theta);

    options = optimoptions('fsolve',...
        'Jacobian','on',...
        'TolFun',TolFun,...
        'Display','iter');
    v = fsolve(@eqn,v0,options);

%
% This is the code to do a manual newton-raphson search instead of fsolve
%
%     v = v0;
%     fprintf(' iter     res_ncon        res_Ncon\n');
%     fprintf(' -----------------------------------------\n');
%     
%     maxiter = 20;
%     tol = 1e-15;
%     optimality = 1;
%     iter = 0; 
%     while iter<=maxiter && optimality>tol;
%         [err,grad] = eqn(v);  
%         optimality = max(abs(err));
%         dv = - grad\err;
%         v = v + dv; 
%         fprintf('   %i    %e    %e\n',iter,optimality,sum(err));
%         iter = iter+1;   
%     end
    
    function [err,grad] = eqn(v)
        [n,grad] = solver(mu,v,vL+vsL,vR+vsR);
        err = n-n0;
    end
end

