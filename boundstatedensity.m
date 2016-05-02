function [ n ] = boundstatedensity( E, solution, dx )
    %BOUNDSTATEDENSITY calculate orbital density with energy E.
  
    phiL = solution(:,1);
    phiR = solution(:,3);

    alpha = 2*sqrt(-2*E);
    int = (phiL(1).*phiR(1)+phiL(end).*phiR(end))/alpha;
    C = sum(phiL.*phiR)*dx + int;
    n = phiL.*phiR./C;  
    
end

