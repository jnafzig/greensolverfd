function testPoschlTeller( testCase )
    %TESTPOSCHLTELLER test using the poschlteller potential
    

    Nelem = 400;

    R = -5;
    B = 5;
    x = linspace(R,B,Nelem)';
    dx = x(2)-x(1);

    v = -cosh(x).^-2;
    vL = 0;
    vR = 0;

    shoot = shoot_fh(Nelem,dx);
    wron = @(E) mean(wronskian(shoot(E,v,vL,vR)));

    E0 = fzero(wron,-.5);

    Check = E0+1/2;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx^2);
    
    solution = shoot(E0,v,vL,vR);
    
    phiL = solution{1}(:,1);
    dphiL = solution{2}(:,1);
    phiR = solution{1}(:,2);
    dphiR = solution{2}(:,2);
    dphiL = (dphiL(1:end-1)+dphiL(2:end))/2;
    dphiR = (dphiR(1:end-1)+dphiR(2:end))/2;
    
    W = phiL.*dphiR - dphiL.*phiR;

    phi = (1-tanh(x).^2).^(1/2);

    Check = (phi/phi(1)-phiL)./(phi/phi(1));
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx);
    
    Check = W-mean(W);
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);
    
    Check = phiL-phiR;
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e-10);
    
    dens = @(E) ldos(shoot(E,v,vL,vR));
    
    R = 1/4;
    E = @(theta) E0 + R*exp(1i*theta);
    dEdt = @(theta) 1i*R*exp(1i*theta);

    RelTol = 1e-6;
    AbsTol = 1e-6;

    n = integral(@(theta) dens(E(theta))*dEdt(theta),0,2*pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);

    Check = imag(n);
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',eps*Nelem);
    
    Check = real(n-(1-tanh(x).^2)/2);
    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx^2);
    
    E = @(theta) -1 + R*exp(1i*theta);

    % Integral should be zeros since it includes no poles:
    Check = integral(@(theta) dens(E(theta))*dEdt(theta),0,2*pi,...
                'ArrayValued',true,...
                'RelTol',RelTol,...
                'AbsTol',AbsTol);

    testCase.verifyEqual(max(abs(Check)),0,'AbsTol',eps*Nelem);
    
end

