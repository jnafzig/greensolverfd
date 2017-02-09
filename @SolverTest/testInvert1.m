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

    solver = solver_fh(Nelem,dx);

    mu = -.25;
    n = solver(mu,v,vL,vR);

    vinv = invert({solver},n,mu,v+vdiff,vL,vR,AbsTol);
    ninv = solver(mu,vinv+v+vdiff,vL,vR);
    
    Check = ninv-n;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);
    
    Check = vinv+vdiff;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);

end
