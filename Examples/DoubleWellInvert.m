tic;

RelTol = eps;
AbsTol = eps;

A = -6;
B = 6;
Nelem = 100;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

R = 3;
v1 = -cosh(x+R/2).^-2;
v2 = -cosh(x-R/2).^-2;
vL = 0;
vR = 0;

solver = solver_fh(Nelem,dx);

mu = -.25;

n1 = solver(mu,v1,0,0);
n2 = solver(mu,v2,0,0);
nf = n1+n2;

vinv = invert(solver,nf,mu,v1+v2,vL,vR,eps);
ninv = solver(mu,vinv,0,0);

subplot(2,2,1);
plot(x,[n1,n2,ninv,nf]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v1,v2,vinv,vinv-v1-v2]);

ylim([min(vinv)-.2,max(vinv)+.2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,[vinv-v1-v2]);

ylim([min(vinv-v1-v2)-.03,max(vinv-v1-v2)+.03]);
xlim([min(x),max(x)]);

title('vinv-v1-v2');

subplot(2,2,2);

plot(x,[ninv-nf]);
xlim([min(x),max(x)]);

title('density error');

toc;