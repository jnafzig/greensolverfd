function testInvert1( testCase )

    RelTol = eps;
    AbsTol = eps;

    A = 0;
    B = pi;
    Nelem = 100;
    x = linspace(A,B,Nelem)';
    dx = x(2)-x(1);

    v = -ones(Nelem,1);
    vdiff = 1e-2*cos(x); 
    vL = -1;
    vR = 0;

    shoot = solver_fh(Nelem,dx);

    dens = @(E,v) ldos(shoot(E,v,vL,vR));
    resp = @(E,v) response(shoot(E,v,vL,vR));

    E0 = min(v)-1;
    mu = -.25;
    kf = sqrt(-2*mu);

    R = (E0+mu)/2;
    A = mu-R;
    E = @(theta) R + A*exp(1i*theta);
    dEdt = @(theta) 1i*A*exp(1i*theta);

    tic;
    n = integral(@(theta) dens(E(theta),v)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    n = n+conj(n);
    toc;

    tic;
    vinv = invert(dx,n,mu,v+vdiff,vL,vR,AbsTol);
    toc;

    tic;
    ninv = integral(@(theta) dens(E(theta),vinv)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    ninv = ninv+conj(ninv);
    toc;
    
    Check = ninv-n;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);
    
    Check = vinv-v;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);

end
