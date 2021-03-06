clear;

A = -12;
B = 12;
Nelem = 600;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

R = 3;
v1 = -cosh(x+R/2).^-2;
v2 = -cosh(x-R/2).^-2;
vi = [v1,v2];

bssolver = shootsolver_fh(Nelem,dx);

N1 = .5;
N2 = .5;
Ni = [N1,N2];
nm = bssolver(sum(Ni),sum(vi,2));

vp = invert({bssolver,bssolver},nm,Ni,vi,[0,0],[0,0],eps);

n1 = bssolver(N1,v1+vp);
n2 = bssolver(N2,v2+vp);
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

toc;