clear;

RelTol = eps;
AbsTol = eps;

A = -4;
B = 4;
Nelem = 100;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

v = -10*cosh(x).^-2;
vL = 0;
vR = 0;

solver = solver_fh(Nelem,dx);
bssolver = shootsolver_fh(Nelem,dx);
shootbs = shooteigsolver_fh(Nelem,dx);
contourbs = contoureigsolver_fh(Nelem,dx);
quadbs = quadeigsolver_fh(Nelem,dx);


N = 4;
mu = -.25;

% These routines do not use consistent normalization yet
[Evals,Evecs] = shootbs(N,v);
[Evalsc,Evecsc] = contourbs(N,v);
[Evalsq,Evecsq] = quadbs(N,v);

nb = bssolver(N,v);
n = solver(mu,v,vL,vR);

subplot(2,2,1);
plot(x,[n,nb]);
xlim([min(x),max(x)]);

title('densities from fixed N and fixed mu');

subplot(2,2,3);
plot(x,v,'r');
for i = 1:N
    line(xlim,[Evals(i),Evals(i)],'Color',[.5,.5,.5],'LineStyle','--');
end

ylim([min(v)-1/2,max(v)+1/2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,2);

plot(x,[n-nb]);
xlim([min(x),max(x)]);

title('density difference between fixed N and fixed mu');

subplot(2,2,4);

E0 = min(v)-1;
R = (E0+mu)/2;
A = mu-R; 
E = @(theta) R + A*exp(1i*theta);
plot(E(linspace(0,2*pi,200)),'Color',[.5,.5,.5],'LineStyle','--');
hold on;
plot(Evals,imag(Evals),'x','LineStyle', 'none')
hold off;

title('complex contour and poles in the complex energy plane');