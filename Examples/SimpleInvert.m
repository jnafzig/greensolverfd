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

shoot = solver_fh(Nelem,dx);

dens = @(E,v) ldos(shoot(E,v,vL,vR));
resp = @(E,v) response(shoot(E,v,vL,vR));

E0 = min(v+vdiff)-1;
mu = -.25;

R = (E0+mu)/2;
A = mu-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);

tic;
n = integral(@(theta) dens(E(theta),v+vdiff)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n = n+conj(n);
toc;

tic;
n0 = integral(@(theta) dens(E(theta),v)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n0 = n0+conj(n0);
toc;

tic;
vinv = invert(dx,n,mu,v,vL,vR,eps);
toc;

tic;
ninv = integral(@(theta) dens(E(theta),vinv)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
ninv = ninv+conj(ninv);
toc;

subplot(2,2,1);
plot(x,[n,ninv,n0]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v+vdiff,vinv,v]);

ylim([min(vinv)-.2,max(vinv)+.2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,[vinv-(v+vdiff)]);

xlim([min(x),max(x)]);

title('vinv-(v+vdiff)');

subplot(2,2,2);

plot(x,ninv-n);
xlim([min(x),max(x)]);

title('density error');

