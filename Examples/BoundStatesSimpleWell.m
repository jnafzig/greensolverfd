clear;

RelTol = eps;
AbsTol = eps;

A = -7;
B = 7;
Nelem = 100;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

v = -10*cosh(x).^-2;
vL = 0;
vR = 0;

shoot = solver(Nelem,dx);
boundstates = eigsolver(Nelem,dx);
dens = @(E) ldos(shoot(E,v,vL,vR));

N = 4;
tic;
Evals = boundstates(N,v);
toc;

nb = zeros(Nelem,1);
for i = 1:N
    nb = nb + boundstatedensity(Evals(i),shoot(Evals(i),v,vL,vR), dx);
end

E0 = min(v)-1;
mu = -.25;
kf = sqrt(-2*mu);

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

subplot(2,2,1);
plot(x,[n,nb]);
xlim([min(x),max(x)]);

title('densities from fixed N and fixed mu');

subplot(2,2,3);
plot(x,v,'r');
for i = 1:N
    line(xlim,[Evals(i),Evals(i)],'Color',[.5,.5,.5],'LineStyle','--');
end

ylim([min(v)-1/2,max(v)+1/2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,2);

plot(x,[n-nb]);
xlim([min(x),max(x)]);

title('density difference between fixed N and fixed mu');

subplot(2,2,4);


plot(E(linspace(0,2*pi,200)),'Color',[.5,.5,.5],'LineStyle','--');
hold on;
plot(Evals,imag(Evals),'x','LineStyle', 'none')
hold off;

title('complex contour and poles in the complex energy plane');