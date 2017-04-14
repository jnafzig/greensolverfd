function [ logtfit_fcn ] = logtfit_fcn( ~ )
    % Return a fit function for lsqcurvefit
    % in particular it returns the log of a skewed student's t-distribution
    
    student = makedist('tlocationscale');
    
    % Distribution is centered on 0 and shifted manually
    student.mu = 0;
        
    logtfit_fcn = @logtfit;
    
    function logtfit = logtfit(x,xdata)
        % x contains an array of the parameters: x -> [mu,sigma,nu,skew]
        % logtfit will return the log of the skewed student's
        % t-distribution with parameters contained in x, at the values
        % contained in xdata.
        
        student.sigma = x(2);
        student.nu = x(3);
        
        logtfit = log(2*student.pdf(xdata-x(1)).*student.cdf(x(4)*(xdata-x(1))));
    end
end

