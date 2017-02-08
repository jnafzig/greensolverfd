tic;

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
vL = 0;
vR = 0;

shoot = solver_fh(Nelem,dx);

dens = @(E,v) ldos(shoot(E,v,vL,vR));
resp = @(E,v) response(shoot(E,v,vL,vR));

E0 = min(v1+v2)-1;
mu = -.25;
kf = sqrt(-2*mu);

R = (E0+mu)/2;
A = mu-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);

n1 = integral(@(theta) dens(E(theta),v1)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n1 = n1+conj(n1);

n2 = integral(@(theta) dens(E(theta),v2)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n2 = n2+conj(n2);
nf = n1+n2;

vinv = invert(dx,nf,mu,v1+v2,vL,vR,eps);

ninv = integral(@(theta) dens(E(theta),vinv)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
ninv = ninv+conj(ninv);

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

toc;