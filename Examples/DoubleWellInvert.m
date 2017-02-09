clear;

A = -6;
B = 6;
Nelem = 100;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

R = 3;
v1 = -cosh(x+R/2).^-2;
v2 = -cosh(x-R/2).^-2;
vi = [v1,v2];
vLi = [0,0];
vRi = [0,0];

solver = solver_fh(Nelem,dx);

mu = -.25;

n1 = solver(mu,v1,0,0);
n2 = solver(mu,v2,0,0);
nf = n1+n2;

vinv = invert({solver},nf,mu,sum(vi,2),sum(vLi),sum(vRi),eps);
ninv = solver(mu,vinv+v1+v2,0,0);

subplot(2,2,1);
plot(x,[n1,n2,ninv,nf]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v1,v2,vinv]);

ylim([min(vinv+v1+v2)-.2,max(vinv+v1+v2)+.2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,vinv);

ylim([min(vinv)-.03,max(vinv)+.03]);
xlim([min(x),max(x)]);

title('vinv');

subplot(2,2,2);

plot(x,ninv-nf);
xlim([min(x),max(x)]);

title('density error');
