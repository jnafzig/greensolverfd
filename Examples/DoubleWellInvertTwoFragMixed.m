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
bssolver = boundstatesolver_fh(Nelem,dx);

mu = -.4;
nm = solver(mu,sum(vi,2),sum(vLi),sum(vRi));
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

ylim([min(v1+v2)-.2,max(v1+v2)+.2]);
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
