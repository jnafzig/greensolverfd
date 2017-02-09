RelTol = eps;
AbsTol = eps;

A = 0;
B = pi;
Nelem = 100;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

v = -ones(Nelem,1);
vdiff = cos(x)/6;
vL = 0;
vR = 0;

solver = solver_fh(Nelem,dx);

mu = -.25;
n = solver(mu,v+vdiff,vL,vR);
n0 = solver(mu,v,vL,vR);

vinv = invert({solver},n,mu,v,vL,vR,eps);
ninv = solver(mu,vinv+v,vL,vR);

subplot(2,2,1);
plot(x,[n,ninv,n0]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v+vdiff,vinv+v,v]);

ylim([min(vinv+v)-.2,max(vinv+v)+.2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,vinv-(v+vdiff));

xlim([min(x),max(x)]);

title('vinv-(v+vdiff)');

subplot(2,2,2);

plot(x,ninv-n);
xlim([min(x),max(x)]);

title('density error');

