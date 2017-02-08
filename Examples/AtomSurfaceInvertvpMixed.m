
RelTol = eps;
AbsTol = eps;

A = -15;
B = 10;
Nelem = 400;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

vL1 = -1;
vR1 = 0;
R = 2;
v1 = zeros(Nelem,1);
v1(x<0) = vL1;

N2 = 1;
v2 = -cosh(x-R).^-2;

solver = solver_fh(Nelem,dx);
bssolver = boundstatesolver_fh(Nelem,dx);

mu = -.25;
nm = solver(mu,v1+v2,vL1+vL2,vR1+vR2);

vp = invertvpmixed(solver,bssolver,nm,mu,v1,vL1,vR1,N2,v2,eps);

n1 = solver(mu,v1+vp,vL1,vR1);
n2 = bssolver(N2,v2+vp);

nf = n1+n2;

subplot(2,2,1);
plot(x,[n1,n2,nm,nf]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v1,v2,vp]);

ylim([min(v1+v2)-.2,max(vp+.2)]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,vp,'r');

ylim([min(vp)-.03,max(vp)+.03]);
xlim([min(x),max(x)]);

title('vp');

subplot(2,2,2);

plot(x,nf-nm);
xlim([min(x),max(x)]);

title('density error');