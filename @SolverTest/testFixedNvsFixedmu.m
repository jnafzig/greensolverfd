function testFixedNvsFixedmu( testCase )
    %TESTRESPONSE test whether fixed mu calculation matches density of
    %fixed N calculation
    
    Tol = eps;
    
    A = -7;
    B = 7;
    Nelem = 100;
    x = linspace(A,B,Nelem)';
    dx = x(2)-x(1);

    v = -10*cosh(x).^-2;
    vL = 1;
    vR = -.1;

    solver = solver_fh(Nelem,dx,Tol);
    bssolver = shootsolver_fh(Nelem,dx);

    N = 4;
    tic;
    nb = bssolver(N,v,vL,vR);
    toc;

    mu = -.25;
    tic;
    n = solver(mu,v,vL,vR);
    toc;
    
    % The density should be the same whether calculated using bound states
    % or fixed chemical potential
    Check = n-nb;

    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',2*Nelem*eps);

end

