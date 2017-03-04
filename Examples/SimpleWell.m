
clear;

RelTol = 1e-10;
AbsTol = 1e-10;

A = -5;
B = 5;
Nelem = 100;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

v = -cosh(x).^-2;
vL = 0;
vR = 0;

shoot = shoot_fh(Nelem,dx);

dens = @(E) ldos(shoot(E,v,vL,vR));
green = @(E) greens(shoot(E,v,vL,vR));
resp = @(E) response(shoot(E,v,vL,vR));

E0 = min(v)-1;
mu = -.5;
kf = sqrt(-2*mu);

R = (E0+mu)/2;
A = mu-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);

% Ek = @(ki) (kf+ki).^2/2;
% dEdk = @(ki) (kf+ki);


tic;
n = integral(@(theta) dens(E(theta))*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n = n+conj(n);
toc;


tic;
chi = integral(@(theta) resp(E(theta))*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
chi = chi+conj(chi);
toc;

tic;
densmat = integral(@(theta) green(E(theta))*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
densmat = densmat+conj(densmat);
toc;

subplot(2,2,1);
plot(x,n);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,v,'r');

ylim([min(v)-1/2,max(v)+1/2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
contourf(x,x,chi);
colorbar;

title('response');

subplot(2,2,2);
contourf(x,x,densmat); 
colorbar;

title('density matrix');

