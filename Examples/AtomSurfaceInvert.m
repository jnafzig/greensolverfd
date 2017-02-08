
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

shoot = solver_fh(Nelem,dx);

dens = @(E,v,vL) ldos(shoot(E,v,vL,vR));
resp = @(E,v,vL) response(shoot(E,v,vL,vR));

E0 = min(v1+v2)-1;
mu = -.25;
kf = sqrt(-2*mu);

R = (E0+mu)/2;
A = mu-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);

tic;
n1 = integral(@(theta) dens(E(theta),v1,vL)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n1 = n1+conj(n1);
toc;
tic;
n2 = integral(@(theta) dens(E(theta),v2,0)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n2 = n2+conj(n2);
toc;
nf = n1+n2;

tic;
vinv = invert(dx,nf,mu,v1+v2,vL,vR,eps);
toc;

tic;
ninv = integral(@(theta) dens(E(theta),vinv,vL)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
ninv = ninv+conj(ninv);
toc;
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

