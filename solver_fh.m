function [ solver_fh ] = solver_fh(Nelem, dx)
    %SOLVER_FH provides function handle for solver function
    %
    % Solver method shoots solutions from both left and right boundaries. 
    % the output is a matrix composed of four column vectors: phi_L, dphi_L
    % and phi_R, dphi_R.  The phi vectors are solutions to the schrodinger
    % equation which satisfy either the left or right boundary condition
    % and the dphi are the corresponding spatial derivatives.
    %
    % While solving each solution is stored in a single column solution 
    % vector sol = [phi;dphi].  Where phi has Nelem elements and dphi has
    % Nelem + 1 elements, so that each element of phi spatially fits in 
    % between two dphi elements.
    %
    
    loc = -1/2:1/2;
    
    % Average operator for dphi
    coeff = fd_coeff(loc,0,dx);
    D0 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1));
    
    MID = spdiags(ones([Nelem,1]), 1,Nelem-1,Nelem+1);
    
    % D1 is forward difference for phi and D1dphi is forward difference for
    % dphi
    coeff = fd_coeff(loc,1,dx);
    D1 = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem-1,Nelem));       
    D1dphi = (spdiags(repmat(coeff',[Nelem,1]), 0:1,Nelem,Nelem+1)); 
        
    % left boundary conditions specify leftmost phi and dphi values 
    bcL = [[1,zeros(1,2*Nelem)];...
           [zeros(1,Nelem),1,zeros(1,Nelem)]];

    % right boundary conditions specify rightmost phi and dphi values
    bcR = [[zeros(1,Nelem-1),1,zeros(1,Nelem+1)];...
           [zeros(1,2*Nelem),1]];

    % lhs operators for shooting method:
    % top two rows specify boundary conditions: bcL and bcR
    % next Nelem-1 rows require that dphi is the spatial derivative of phi:
    % [D1,-MID] and the last Nelem rows specify shrodinger equation.  the 
    % potential and energy of the solution are added in later using the 
    % matVE matrix.
    matL = [bcL;...
           [D1,-MID];...
           [zeros(Nelem),D1dphi]];
    matR = [bcR;...
           [D1,-MID];...
           [zeros(Nelem),D1dphi]];

    
    matVE =  spdiags(ones(Nelem,1),-Nelem-1,2*Nelem+1,2*Nelem+1);
       
    [ival,jval] = find(matVE);
    
    solver_fh = @(E,v,vL,vR) shoot(E,v,vL,vR);

    function [solution] = shoot(E,v,vL,vR)
        
        matVE = sparse(ival,jval,2*(E-v),2*Nelem+1,2*Nelem+1,Nelem);

        kL = sqrt(2*(E-vL));
        if imag(kL)<0; kL = -kL; end  % ensure we choose the sqrt with 
                                      % negative imaginary part
        
        % first two elements of rhs are specifying boundary conditions for
        % phi and dphi.
        rhsL = [1;-1i*kL;zeros(2*Nelem-1,1)];
        lhsL = matL+matVE;
        solL = lhsL\rhsL;
        phiL = solL(1:Nelem);
        dphiL = D0*solL((1:Nelem+1)+Nelem);
        
        % Same but for right boundary condition.
        kR = sqrt(2*(E-vR)); 
        if imag(kR)<0; kR = -kR; end
        
        lhsR = matR+matVE;
        rhsR = [1;1i*kR;zeros(2*Nelem-1,1)]; 
        solR = lhsR\rhsR;
        phiR = solR(1:Nelem);
        dphiR = D0*solR((1:Nelem+1)+Nelem);
        
        solution = [phiL,dphiL,phiR,dphiR];

    end
end
