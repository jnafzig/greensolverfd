function [ logtfit_fcn ] = logtfit_fcn( ~ )
    %LOSS Summary of this function goes here
    %   Detailed explanation goes here
    
%     
%     lb = -1;
%     ub = -0.42;
% 
%     student = truncate(makedist('tlocationscale'),lb,ub);
    student = makedist('tlocationscale');
    student.mu = 0;
        
    
    logtfit_fcn = @logtfit;
    
    function logtfit = logtfit(x,xdata)
        student.sigma = x(2);
        student.nu = x(3);
        
        logtfit = log(2*student.pdf(xdata-x(1)).*student.cdf(x(4)*(xdata-x(1))));
    end
end

