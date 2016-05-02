function testInvert3( testCase )

    RelTol = eps;
    AbsTol = eps;

    A = -6;
    B = 6;
    Nelem = 100;
    x = linspace(A,B,Nelem)';
    dx = x(2)-x(1);

    vL = -1;
    vR = 0;
    R = 2;
    v1 = zeros(Nelem,1);
    v1(x<0) = vL;
    v2 = -cosh(x-R).^-2;


    shoot = solver(Nelem,dx);

    dens = @(E,v,vL) ldos(shoot(E,v,vL,vR));
    resp = @(E,v,vL) response(shoot(E,v,vL,vR));

    E0 = min(v1+v2)-1;
    mu = -.25;
    kf = sqrt(-2*mu);

    R = (E0+mu)/2;
    A = mu-R;
    E = @(theta) R + A*exp(1i*theta);
    dEdt = @(theta) 1i*A*exp(1i*theta);

    tic;
    n1 = integral(@(theta) dens(E(theta),v1,vL)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    n1 = n1+conj(n1);
    toc;
    tic;
    n2 = integral(@(theta) dens(E(theta),v2,0)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    n2 = n2+conj(n2);
    toc;
    nf = n1+n2;

    tic;
    vinv = invert(dx,nf,mu,v1+v2,vL,vR);
    toc;

    tic;
    ninv = integral(@(theta) dens(E(theta),vinv,vL)*dEdt(theta),0,pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);
    ninv = ninv+conj(ninv);
    toc;
    
    Check = ninv-nf;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);


end