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

bssolver = boundstatesolver(Nelem,dx);

N1 = 1;
tic;
n1 = bssolver(N1,v1);
toc;

N2 = 1;
tic;
n2 = bssolver(N2,v2);
toc;

nf = n1+n2;
Nf = N1+N2;

tic;
vinv = invertboundstate(bssolver,nf,Nf,v1+v2,eps);
toc;

tic;
ninv = bssolver(Nf,vinv);
toc;

subplot(2,2,1);
plot(x,[n1,n2,ninv,nf]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v1,v2,vinv,vinv-v1-v2]);

ylim([min(v1)-.2,max(vinv)+.2]);
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
