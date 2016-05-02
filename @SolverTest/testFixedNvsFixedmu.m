function testResponse( testCase )
    %TESTRESPONSE test whether fixed mu calculation matches density of
    %fixed N calculation
    
    RelTol = eps;
    AbsTol = eps;

    A = -7;
    B = 7;
    Nelem = 100;
    x = linspace(A,B,Nelem)';
    dx = x(2)-x(1);

    v = -10*cosh(x).^-2;
    vL = 0;
    vR = 0;

    shoot = solver(Nelem,dx);
    boundstates = eigsolver(Nelem,dx);
    dens = @(E) ldos(shoot(E,v,vL,vR));

    N = 4;
    tic;
    Evals = boundstates(N,v);
    toc;

    nb = zeros(Nelem,1);
    for i = 1:N
        nb = nb + boundstatedensity(Evals(i),shoot(Evals(i),v,vL,vR), dx);
    end

    E0 = min(v)-1;
    mu = -.25;

    R = (E0+mu)/2;
    A = mu-R;
    E = @(theta) R + A*exp(1i*theta);
    dEdt = @(theta) 1i*A*exp(1i*theta);

    tic;
    n = integral(@(theta) dens(E(theta))*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    n = n+conj(n);
    toc;
    
    % The density should be the same whether calculated using bound states
    % or fixed chemical potential
    Check = n-nb;

    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',2*Nelem*eps);

     
end

