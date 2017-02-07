function testResponse( testCase )
    %TESTRESPONSE test whether response correctly predicts density
    % response to small change in potential 
    

    % testing set: [E, dv, Nelem]
    testparameters = [[-.25,      1e-3,  10]; ...
                      [-.25,      1e-6,  30]; ...
                      [-.25+1i/2, 1e-4,  20]; ...
                      [-.5+1i/2,  1e-10, 20]; ...
                      [.25,       1e-11, 100]; ...
                      [.25,       1e-3,  10]; ...
                      [.25,       1e-6,  30]; ...
                      [.25+1i/2,  1e-4,  20]; ...
                      [.25+1i/2,  1e-10, 20]; ...
                      [.25+1i/2,  1e-9,  50]];
    Ntest = size(testparameters,1);
    
    for i = 1:Ntest

        E = testparameters(i,1);
        dv = testparameters(i,2);
        Nelem = testparameters(i,3);
                
        A = -2;
        B = 2;
        x = linspace(A,B,Nelem)';
        dx = x(2)-x(1);

        v = -cosh(x).^-2;
        vdiff = dv * rand(size(v));
        v1 = v - vdiff;
        v2 = v + vdiff;
        vL = 0;
        vR = 0;

        shoot = solver(Nelem,dx);

        chi = response(shoot(E,v,vL,vR));

        n1 = ldos(shoot(E,v1,vL,vR));
        n2 = ldos(shoot(E,v2,vL,vR));

        % The change in density should be equal to response * change in
        % potential
        Check = (n2-n1)/2 - chi*vdiff*dx;
         
        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',dv*dx^2);
    end
     
end

