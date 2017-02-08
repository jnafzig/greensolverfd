function  v  = invertvp( solver, n0, mu, v1, vL1, vR1, v2, vL2, vR2, TolFun)
    if nargin < 10
        TolFun = 1e-6;
    end
    
    %INVERT tries to find v such that v -> n0 at a chemical potential mu

    vpR = 0; % currently the partition potential is assumed to be zero 
    vpL = 0; % outside the gridded region.  For now vpR and vpL can be 
             % adjusted manually to play around with that assumption.

    v = zeros(size(n0));

    options = optimoptions('fsolve',...
        'Jacobian','on',...
        'TolFun',TolFun,...
        'Display','iter');
    v = fsolve(@eqn,v,options);
    
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
        [n2,grad2] = solver(mu,v2+vp,vL2+vpL,vR2+vpR);
        
        err = n1+n2-n0;
        grad = grad1+grad2;
        
    end
end

