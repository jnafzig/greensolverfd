function testInvertvp2( testCase )

    RelTol = eps;
    AbsTol = eps;

    A = -6;
    B = 6;
    Nelem = 100;
    x = linspace(A,B,Nelem)';
    dx = x(2)-x(1);

    vL1 = -1;
    vR1 = 0;
    R = 2;
    v1 = zeros(Nelem,1);
    v1(x<0) = vL1;

    vL2 = 0;
    vR2 = 0;
    v2 = -cosh(x-R).^-2;

    solver = solver_fh(Nelem,dx);

    mu = -.25;
    
    nm = solver(mu,v1+v2,vL1+vL2,vR1+vR2);

    vp = invertvp(solver,nm,mu,v1,vL1,vR1,v2,vL2,vR2,AbsTol);

    n1 = solver(mu,v1+vp,vL1,vR1);
    n2 = solver(mu,v2+vp,vL2,vR2);
    nf = n1+n2;
    
    Check = nf-nm;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);

end