function  vp  = invertvpboundstate( bssolver, n0, N1, v1, N2, v2, TolFun)
    if nargin < 7
        TolFun = 1e-6;
    end
    
    %INVERT tries to find v such that v -> n0 at a chemical potential mu

    % currently the partition potential is assumed to be zero 
    % outside the gridded region.  
    
    % bssolver must be function handle for a boundstatesolver
    % n0 is target density
    % N1 is number of bound states in fragment 1
    % v1 is initial guess for the potential
    % N2 is number of bound states in fragment 2
    % v2 is initial guess for the potential
    
    v0 = zeros(size(n0));
    
    options = optimoptions('fsolve',...
        'Jacobian','on',...
        'TolFun',TolFun,...
        'Display','iter');
    vp = fsolve(@eqn,v0,options);
        
    function [err,grad] = eqn(vp)

        [n1,grad1] = bssolver(N1,v1+vp);
        [n2,grad2] = bssolver(N2,v2+vp);
        
        err = n1+n2-n0;
        grad = grad1+grad2;
        
    end
end

