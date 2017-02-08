function  v  = invertvpmixed( solver, bssolver, n0, mu, v1, vL1, vR1, N2, v2, TolFun)
    if nargin < 10
        TolFun = 1e-6;
    end
    
    %INVERTVPMIXED tries to find find vp for cases where on fragment is
    % fixed by chemical potential and the other is fixed by occupation
    % number.
    
    vpR = 0; % currently the partition potential is assumed to be zero 
    vpL = 0; % outside the gridded region.  For now vpR and vpL can be 
             % adjusted manually to play around with that assumption.
    
    RelTol = eps;
    AbsTol = eps;

    v0 = zeros(size(n0));

    options = optimoptions('fsolve',...
        'Jacobian','on',...
        'TolFun',TolFun,...
        'Display','iter');
    v = fsolve(@eqn,v0,options);
    
%    
%     fprintf(' iter     res_ncon        res_Ncon\n');
%     fprintf(' -----------------------------------------\n');
%     
%     maxiter = 20;
%     tol = 1e-13;
%     optimality = 1; 
%     iter = 0; 
%     while iter<=maxiter && optimality>tol;
%         [err,grad] = eqn(v);   
%         n1;n2;
%         optimality = max(abs(err));
%         dv = - grad\err;
%         if abs(dv(end))>1
%             dv(end) = sign(dv(end));
%         end
%         if abs(dv(1))>1
%             dv(end) = sign(dv(1));
%         end
%         
%         v = v + dv; 
%         
%         fprintf('   %i    %e    %e\n',iter,optimality,sum(err));
%         iter = iter+1;   
% end
    
    
    function [err,grad] = eqn(vp)

        [n1,grad1] = solver(mu,v1+vp,vL1+vpL,vR1+vpR);
        [n2,grad2] = bssolver(N2,v2+vp);
        
        err = n1+n2-n0;
        
        grad = grad1+grad2;
        
    end
end

