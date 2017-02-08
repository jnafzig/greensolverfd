function testInvertvp1( testCase )


    RelTol = eps;
    AbsTol = eps;

    A = -10;
    B = 10;
    Nelem = 100;
    x = linspace(A,B,Nelem)';
    dx = x(2)-x(1);

    R = 3;
    v1 = -cosh(x+R/2).^-2;
    v2 = -cosh(x-R/2).^-2;
    vL1 = 0;
    vR1 = 0;
    vL2 = 0;
    vR2 = 0;

    shoot = solver_fh(Nelem,dx);

    dens = @(E,v,vL,vR) ldos(shoot(E,v,vL,vR));
    resp = @(E,v,vL,vR) response(shoot(E,v,vL,vR));

    E0 = min(v1+v2)-1;
    mu = -.25;
    kf = sqrt(-2*mu);

    R = (E0+mu)/2;
    A = mu-R;
    E = @(theta) R + A*exp(1i*theta);
    dEdt = @(theta) 1i*A*exp(1i*theta);

    tic;
    nm = integral(@(theta) dens(E(theta),v1+v2,vL1+vL2,vR1+vR2)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    nm = nm+conj(nm);
    toc;

    tic;
    vp = invertvp(dx,nm,mu,v1,vL1,vR1,v2,vL2,vR2,AbsTol);
    toc;

    tic;
    n1 = integral(@(theta) dens(E(theta),v1+vp,vL1,vR1)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    n1 = n1+conj(n1);
    toc;

    tic;
    n2 = integral(@(theta) dens(E(theta),v2+vp,vL2,vR2)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    n2 = n2+conj(n2);
    toc;
    nf = n1+n2;
    
    
    Check = nf-nm;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);


end
