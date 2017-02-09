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

    solver = solver_fh(Nelem,dx);

    mu = -.25;
    nm = solver(mu,v1+v2,vL1+vL2,vR1+vR2);

    vp = invert({solver,solver},nm,[mu,mu],[v1,v2],[vL1,vL2],[vR1,vR2],AbsTol);
    
    n1 = solver(mu,v1+vp,vL1,vR1);
    n2 = solver(mu,v2+vp,vL2,vR2);
    
    nf = n1+n2;
    
    Check = nf-nm;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);

end
