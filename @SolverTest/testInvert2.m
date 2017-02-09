function testInvert2( testCase )

    RelTol = eps;
    AbsTol = eps;

    A = -6;
    B = 6;
    Nelem = 100;
    x = linspace(A,B,Nelem)';
    dx = x(2)-x(1);

    R = 2;
    v1 = -cosh(x+R/2).^-2;
    v2 = -cosh(x-R/2).^-2;
    vL = 0;
    vR = 0;

    solver = solver_fh(Nelem,dx);

    mu = -.25;
    
    n1 = solver(mu,v1,vL,vR);
    n2 = solver(mu,v2,vL,vR);
    nf = n1+n2;

    vinv = invert({solver},nf,mu,v1+v2,vL,vR,AbsTol);
    ninv = solver(mu,vinv+v1+v2,vL,vR);
    
    Check = ninv-nf;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);

end

