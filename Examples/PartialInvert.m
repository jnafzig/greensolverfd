clear;

A = -20;
B = 20;
Nelem = 500;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

Nfreq = 75;
vfft = zeros(size(x));
vfft((1:Nfreq)+1) = Nelem*rand(Nfreq,1).*exp(2i*pi*rand(Nfreq,1))/Nfreq;
v = ifft(vfft,'symmetric');

vL = 0;
vR = 0;

solver = solver_fh(Nelem,dx);
kedsolver = kedsolver_fh(Nelem,dx);

mu = 1;
msk = abs(x)<10;

[n,ts] = kedsolver(mu,v,0,0);

vfft = zeros(size(x));

npotentials = 1;

valt = zeros(Nelem,npotentials);
nalt = zeros(Nelem,npotentials);
tsalt = zeros(Nelem,npotentials);
for i = 1:npotentials

    Nfreq = randi(100);
    vfft = zeros(size(x));
    vfft((2:Nfreq+1)) = Nelem*rand(Nfreq,1).*exp(2i*pi*rand(Nfreq,1))/Nfreq;
    vdiff = ifft(vfft,'symmetric');

    vL = v(1)+vdiff(1);
    vR = v(end)+vdiff(end);
    
    vinv = partialinvert({solver},n,mu,v+vdiff,vL,vR,msk,eps);
    valt(:,i) = vinv+v+vdiff;

    [nalt(:,i),tsalt(:,i)] = kedsolver(mu,valt(:,i),vR,vR);

end


subplot(2,2,1);
plot(x,[n,nalt]);
xlim([min(x),max(x)]);

title('density');

subplot(2,2,3);
plot(x,[v,valt]);

ylim([min(min(valt))-.2,max(max(valt))+.2]);
xlim([min(x),max(x)]);

title('potential');

subplot(2,2,4);
plot(x,[ts,tsalt]);

ylim([min(ts)-.03,max(ts)+.03]);
xlim([min(x),max(x)]);

title('ts');
