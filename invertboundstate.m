function  v  = invertboundstate( bssolver, n0, N, v0, TolFun )
    if nargin < 5
        TolFun = 1e-6;
    end
    
    %INVERTBOUNDSTATE
    % tries to find v such that v -> n0 for fixed particle number

    % bssolver must be function handle for a boundstatesolver
    % n0 is target density
    % N is number of bound states
    % v0 is initial guess for the potential
    
    options = optimoptions('fsolve',...
        'Jacobian','on',...
        'TolFun',TolFun,...
        'Display','iter');
    v = fsolve(@eqn,v0,options);
%    
%     fprintf(' iter     res_ncon        res_Ncon\n');
%     fprintf(' -----------------------------------------\n');
%     
%     v = v0;
%     maxiter = 20;
%     tol = 1e-13;
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
        
        [n,grad] = bssolver(N,v);
        err = n-n0;
        
    end
end

