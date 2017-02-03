function testBoundStateResponse( testCase )
    %TESTRESPONSE test whether response correctly predicts density
    % response to small change in potential 
    

    % testing set: [E, dv, Nelem]
    testparameters = [[1e-3,  10]; ...
                      [1e-6,  30]; ...
                      [1e-4,  20]; ...
                      [1e-10, 20]; ...
                      [1e-11, 100]];
    Ntest = size(testparameters,1);
    
    for i = 1:Ntest

        dv = testparameters(i,1);
        Nelem = testparameters(i,2);
                
        A = -5;
        B = 5;
        x = linspace(A,B,Nelem)';
        dx = x(2)-x(1);

        v = -cosh(x).^-2;
        vdiff = dv * x.*cosh(x).^-2;
        v1 = v - vdiff;
        v2 = v + vdiff;

        bssolver = boundstatesolver(Nelem,dx);

        N = 1;
        tic;
        [~,chi] = bssolver(N,v);
        n1 = bssolver(N,v1);
        n2 = bssolver(N,v2);
        toc;
        
        % The change in density should be equal to response * change in
        % potential
        Check = (n2-n1)/2 - chi*vdiff*dx;
         
        testCase.verifyEqual(max(abs(Check)),0,...
            'AbsTol',max(2*dv*dx^2,2e-10));
    end
     
end

