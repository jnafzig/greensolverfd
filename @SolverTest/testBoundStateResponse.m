function testBoundStateResponse( testCase )
    %TESTRESPONSE test whether response correctly predicts density
    % response to small change in potential 
    

    Ntest = 10;
    
    for i = 1:Ntest

        dv = 1e-6;
        Nelem = 50;
                
        A = -2;
        B = 2;
        x = linspace(A,B,Nelem)';
        dx = x(2)-x(1);

        v = -cosh(x).^-2;
        vdiff = dv * rand(size(v));
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
            'AbsTol',1e-14);
    end
     
end

