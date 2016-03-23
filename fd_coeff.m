function stencil = fd_coeff(range,order,h)
    % calculates finite difference coefficients for points in range for
    % approximating derivative given by order.
    
    if (nargin<3); h = 1; end
    
    N = numel(range);
    
    if (order+1>N); 
        error('Order is too high for given range'); 
    end
    
    % Powers for polynomials: one polynomial per point starting at the zeroeth
    power = (0:(N-1))';  

    % Form grid
    [exponents,xvalues] = ndgrid(power,range);
    A = xvalues.^exponents;

    % Analytical derivatives at zero.
    b = zeros(N,1);
    b(order+1) = factorial(order);

    % Calculate weights which will be exact for given orders of polynomials
    if license('test','symbolic_toolbox')
        A = sym(A);
        b = sym(b);
        h = sym(h);
        order = sym(order);
        stencil = double((A\b)/h^order);
    else
        stencil = (A\b)/h^order;
    end
 
end
