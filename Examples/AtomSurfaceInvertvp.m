

RelTol = eps;
AbsTol = eps;

A = -6;
B = 10;
Nelem = 200;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);


vL1 = -1;
vR1 = 0;
R = 2;
v1 = zeros(Nelem,1);
v1(x<0) = vL1;

vL2 = 0;
vR2 = 0;
v2 = -cosh(x-R).^-2;

shoot = solver(Nelem,dx);

dens = @(E,v,vL,vR) density(shoot(E,v,vL,vR));
resp = @(E,v,vL,vR) response(shoot(E,v,vL,vR));

E0 = min(v1+v2)-1;
mu = -.25;
kf = sqrt(-2*mu);

R = (E0+mu)/2;
A = mu-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);

tic;
nm = integral(@(theta) dens(E(theta),v1+v2,vL1+vL2,vR1+vR2)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
nm = nm+conj(nm);
toc;

tic;
vp = invertvp(dx,nm,mu,v1,vL1,vR1,v2,vL2,vR2);
toc;

tic;
n1 = integral(@(theta) dens(E(theta),v1+vp,vL1,vR1)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n1 = n1+conj(n1);
toc;

tic;
n2 = integral(@(theta) dens(E(theta),v2+vp,vL2,vR2)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n2 = n2+conj(n2);
toc;
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