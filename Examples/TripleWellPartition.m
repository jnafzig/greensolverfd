clear;

A = -10;
B = 10;
Nelem = 200;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

Nf = 3;
R = 3.5;
v1 = -cosh(x+R).^-2;
v2 = -cosh(x).^-2;
v3 = -cosh(x-R).^-2;
vLi = zeros(1,3);
vRi = zeros(1,3);
vi = [v1,v2,v3];
v = sum(vi,2);

solver = solver_fh(Nelem,dx);
bssolver = shootsolver_fh(Nelem,dx);

mu = -.25;
nm = solver(mu,v,sum(vLi),sum(vRi));
N = 1;

vp = invert({solver,bssolver,solver},nm,[mu,N,mu],vi,vLi,vRi,eps);

n1 = solver(mu,v1+vp,vLi(1),vRi(1));
n2 = bssolver(N,v2+vp);
n3 = solver(mu,v3+vp,vLi(1),vRi(1));
nf = n1+n2+n3;

subplot(2,2,1);
plot(x,[n1,n2,n3,nm,nf]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v1,v2,v3,vp]);

ylim([min(v1+v2)-.2,max(v1+v2)+.2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,vp,'r');

ylim([min(vp)-.003,max(vp)+.003]);
xlim([min(x),max(x)]);

title('vp');

subplot(2,2,2);

plot(x,nf-nm);
xlim([min(x),max(x)]);

title('density error');

toc;