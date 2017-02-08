
RelTol = eps;
AbsTol = eps;

A = -7;
B = 10;
Nelem = 200;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

vL = -1;
vR = 0;
R = 2;
v1 = zeros(Nelem,1);
v1(x<0) = vL;
v2 = -cosh(x-R).^-2;

solver = solver_fh(Nelem,dx);

mu = -.25;

n1 = solver(mu,v1,vL,vR);
n2 = solver(mu,v2,0,vR);
nf = n1+n2;

vinv = invert(solver,nf,mu,v1+v2,vL,vR,eps);
ninv = solver(mu,vinv,vL,vR);

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

