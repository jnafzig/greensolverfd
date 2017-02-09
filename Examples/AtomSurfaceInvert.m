clear;

A = -7;
B = 10;
Nelem = 200;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

vLi = [-1,0];
vRi = [0,0];
R = 2;
v1 = zeros(Nelem,1);
v1(x<0) = vLi(1);
v2 = -cosh(x-R).^-2;
vi = [v1,v2];

solver = solver_fh(Nelem,dx);

mu = -.25;

n1 = solver(mu,v1,vLi(1),vRi(1));
n2 = solver(mu,v2,vLi(2),vRi(2));
nf = n1+n2;

vinv = invert({solver},nf,mu,v1+v2,sum(vLi),sum(vRi),eps);
ninv = solver(mu,vinv+v1+v2,sum(vLi),sum(vRi));

subplot(2,2,1);
plot(x,[n1,n2,ninv,nf]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v1,v2,vinv,vinv+v1+v2]);

ylim([min(vinv+v1+v2)-.2,max(vinv+v1+v2)+.2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,[vinv]);

ylim([min(vinv)-.03,max(vinv)+.03]);
xlim([min(x),max(x)]);

title('vinv-v1-v2');

subplot(2,2,2);

plot(x,[ninv-nf]);
xlim([min(x),max(x)]);

title('density error');

