function  v  = partialinvert( solver_array, n0, N_mu_array, vi, vLi, vRi,msk,TolFun,v0)
    %PARTIALINVERT( solver_array, n0, N_mu_array, vi, vLi, vRi,TolFun)
    % Attempts to invert the density in the locations indicated by msk
    % solver_array -> Cell array of function handles to bound state or
    % continuum solvers.
    % n0 -> Target density
    % N_mu_array -> Array of occupation numbers, N, and chemical potentials
    % mu, of the same size as solver_array.  Values will be interpreted as
    % N or mu depending on function_handle in solver_array.
    % msk -> logical array of the same size as n0
    
    if nargin < 8
        TolFun = 1e-6;
    end
    
    numfragments = numel(solver_array);
    nelem = numel(n0);
    if nargin < 6
        vRi = zeros(size(N_mu_array));
    end
    if nargin < 5
        vLi = zeros(size(N_mu_array));
    end
    
    % Check sizes
    assert(all(size(n0) == [size(vi,1),1]), 'n0 size does not match vi');
    assert(numel(N_mu_array) == numfragments, 'N_mu_array is the wrong size');
    assert(size(vi,2) == numfragments, 'vi is wrong size');
    assert(numel(vLi) == numfragments, 'vLi is wrong size');
    assert(numel(vRi) == numfragments, 'vRi is wrong size');
    
    vpR = 0; % currently the partition potential is assumed to be zero 
    vpL = 0; % outside the gridded region.  For now vpR and vpL can be 
             % adjusted manually to play around with that assumption.
    
    if nargin < 9
        v0 = zeros(size(n0(msk)));
    else
        v0 = v0(msk);
    end

    options = optimoptions('fsolve',...
        'Jacobian','on',...
        'TolFun',TolFun,...
        'Display','iter');
    v = fsolve(@eqn,v0,options);
    
    function [err,grad] = eqn(vp)

        vpfill = zeros(nelem,1);
        vpfill(msk) = vp;
        n = zeros(nelem,1);
        grad = zeros(nelem);
        for i = 1:numfragments;
            [ni,gradi] = solver_array{i}...
                (N_mu_array(i),vi(:,i)+vpfill,vLi(i)+vpL,vRi(i)+vpR);
            n = n + ni;
            grad = grad + gradi;
        end
        err = n(msk)-n0(msk);
        grad = grad(msk,msk); 
        
    end
end

