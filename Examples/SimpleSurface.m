
clear;

RelTol = 1e-6;
AbsTol = 1e-6;

A = -10;
B = 5/2;
Nelem = 200;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

v = -cosh(x).^-2;

vL = -1;
vR = 0;
R = 5;
v(x<-R) = vL;

shoot = solver(Nelem,dx);

dens = @(E) ldos(shoot(E,v,vL,vR));
green = @(E) greens(shoot(E,v,vL,vR));
resp = @(E) response(shoot(E,v,vL,vR));

E0 = min(v)-1;
mu = -.25;

R = (E0+mu)/2;
A = mu-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);

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
