clear;

% Use of Aitken's delta squared process to improve eigenvalue estimate

ColOrd = get(gca,'ColorOrder');

A = -15;
B = 15;

NelemValues = 6:13;
NumVals = numel(NelemValues);

x = cell(NumVals,1);
Evecs = cell(NumVals,1);
Evals0 = zeros(NumVals,1);
Eval1 = zeros(NumVals,1);
dx = zeros(NumVals,1);

for i = 1:NumVals;
    Nelem = 2^NelemValues(i);
    x{i} = linspace(A,B,Nelem)';
    dx(i) = x{i}(2)-x{i}(1);

    lambda = 1;
    v = -lambda*(lambda+1)/2*cosh(x{i}).^-2;
    vL = 0;
    vR = 0;

    solver = shooteigsolver_fh(Nelem,dx(i));

    N = 1;
    mu = -.25;

    [Evals0(i),Evecs{i}] = solver(N,v);
    
end
hold off;


dvals0 = diff(Evals0);
d2vals0 = diff(Evals0,2);
Evals1 = (Evals0(1:end-2)-dvals0(1:end-1).^2./d2vals0);


loglog(2.^NelemValues,(abs(Evals0+.5)),'+-',...
    (2.^NelemValues(3:end)),(abs(Evals1+.5)),'+-');

xlabel('number of grid points');
ylabel('error');
