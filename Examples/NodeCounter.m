clear;

RelTol = eps;
AbsTol = eps;

A = -10;
B = 10;
Nelem = 100;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

lambda = 3;
R = 6;
v = -lambda*(lambda+1)/2*cosh(x-R/2).^-2 ...
    -lambda*(lambda+1)/2*cosh(x+R/2).^-2;
vL = 0;
vR = 0;

shoot = shoot_fh(Nelem,dx);
bssolver = shooteigsolver_fh(Nelem,dx);

N = 7;
[Evals,Evecs] = bssolver(N,v,vL,vR);

numE = 1000;
Espace = linspace(min(v),0,numE);
nodes = zeros(numE,1);

for i = 1:numE
    nodes(i) = nodecount(shoot(Espace(i),v,vL,vR));
end

plot(Espace,nodes);


for i = 1:N
    line([Evals(i),Evals(i)],ylim,'Color',[.5,.5,.5],'LineStyle','--');
end

xlabel('Energy');
ylabel('Number of States');