clear;

A = -15;
B = 10;
Nelem = 400;
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
bssolver = quadsolver_fh(Nelem,dx);

mu = -.25;
nm = solver(mu,v1+v2,sum(vLi),sum(vRi));
N = 1;

vp = invert({solver,bssolver},nm,[mu,N],vi,vLi,vRi,eps);

n1 = solver(mu,v1+vp,vLi(1),vRi(1));
n2 = bssolver(N,v2+vp);

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