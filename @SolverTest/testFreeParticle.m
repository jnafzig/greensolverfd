function testFreeParticle( testCase )
    %TESTPOSCHLTELLER test using the poschlteller potential
    

    % testing set: [E, dv, Nelem]
    testparameters = [[-1/2,        100]; ...
                      [1/2,         300]; ...
                      [-1/2 + 1i/2, 201]; ...
                      [-1/2 - 1i/2, 200]; ...
                      [1/2 - 1i/2,  101]; ...
                      [1/2 + 1i/2,  100]; ...
                      [1i/2,        301]; ...
                      [-1i/2,       200]; ...
                      [1/4+1i/2,    201]; ...
                      [1/4-1i/2,    501]];
    Ntest = size(testparameters,1);
    
    for i = 1:Ntest

        E = testparameters(i,1);
        Nelem = testparameters(i,2);

        A = 0;
        B = 1;
        x = linspace(A,B,Nelem)';
        dx = x(2)-x(1);

        v = zeros(size(x));
        vL = 0;
        vR = 0;

        shoot = shoot_fh(Nelem,dx);

        solution = shoot(E,v,vL,vR);

        wron = wronskian(solution);

        phiL = solution(:,1);
        dphiL = solution(:,2);
        phiR = solution(:,3);
        dphiR = solution(:,4);

        W = phiL.*dphiR - dphiL.*phiR;

        k = sqrt(2*E);
        if imag(k)<0; k = -k; end

        Check = dphiL./phiL+1i*k;
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx);

        ikxLcheck = -1i*k*x+1i*k*x(1);
        Check = real(log(phiL)-ikxLcheck)./abs(ikxLcheck);
        Check(1) = 0;
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx);
        Check = unwrap(imag(log(phiL)-ikxLcheck))./abs(ikxLcheck);
        Check(1) = 0;
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx);

        Check = dphiR./phiR-1i*k;
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx);

        ikxRcheck = 1i*k*x-1i*k*x(end);
        Check = real(log(phiR)-ikxRcheck)./abs(ikxRcheck);
        Check(end) = 0;
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx);
        Check = max(abs(unwrap(imag(log(phiR)-ikxRcheck))./abs(ikxRcheck)));
        Check(end) = 0;
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dx);

        Check = max(abs((wron-mean(wron))/mean(wron)));
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',eps*Nelem);

    end
    
end

