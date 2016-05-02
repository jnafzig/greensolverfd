function  v  = invert( dx, n0, mu, v0, vL, vR )
    %INVERT tries to find v such that v -> n0 at a chemical potential mu

    vsR = 0; % currently the inverted potential is assumed to be zero 
    vsL = 0; % outside the gridded region.  For now vsR and vsL can be 
             % adjusted manually to play around with that assumption.
    
    RelTol = eps;
    AbsTol = eps;
    
    Nelem = numel(n0);
    shoot = solver(Nelem,dx);

    dens = @(E,v,vL,vR) ldos(shoot(E,v,vL,vR));
    resp = @(E,v,vL,vR) response(shoot(E,v,vL,vR));

    E0 = min(v0)-1;

    R = (E0+mu)/2;
    A = mu-R;
    E = @(theta) R + A*exp(1i*theta);
    dEdt = @(theta) 1i*A*exp(1i*theta);

    v = v0;
    
    fprintf(' iter     res_ncon        res_Ncon\n');
    fprintf(' -----------------------------------------\n');
    
    maxiter = 20;
    tol = 1e-15;
    optimality = 1;
    iter = 0; 
    while iter<=maxiter && optimality>tol;
        [err,grad] = eqn(v);  
        optimality = max(abs(err));
        dv = - grad\err;
        v = v + dv; 
        fprintf('   %i    %e    %e\n',iter,optimality,sum(err));
        iter = iter+1;   
    end
    
    function [err,grad] = eqn(v)

        n = integral(@(theta) dens(E(theta),v,vL+vsL,vR+vsR)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
        n = n+conj(n);

        grad = dx*integral(@(theta) resp(E(theta),v,vL+vsL,vR+vsR)*dEdt(theta),0,pi,...
                    'ArrayValued',true,...
                    'RelTol',1e-1,...
                    'AbsTol',1e-1);
        grad = grad+conj(grad);

        err = n-n0;
        
    end
end

