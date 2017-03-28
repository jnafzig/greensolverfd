clear;

A = -15;
B = 15;
Nelem = 400;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

vLi = [-1,0];
vRi = [0,0];
R = 5;
va = zeros(Nelem,1);
va(x<0) = vLi(1);
vb = -cosh(x-R).^-2;
vi = [va,vb];

solver = solver_fh(Nelem,dx);
bssolver = shootsolver_fh(Nelem,dx);
evsolver = shooteigsolver_fh(Nelem,dx);

nummu = 1;
% muspace = linspace(-.6,-.4,nummu);
muspace = -.5;

Nrec = zeros(nummu,1);
mub = zeros(nummu,1);
vprec = zeros(Nelem,nummu);
narec = zeros(Nelem,nummu);
nbrec = zeros(Nelem,nummu);

for i = 1:nummu

    mu = muspace(i);
    nm = solver(mu,va+vb,sum(vLi),sum(vRi));

    n10 = solver(mu,va,vLi(1),vRi(1));
    n20 = bssolver(1,vb);
    N = sum((nm-n10).*n20)/sum(n20.^2);

    N = max(min(N,1),0);
    n20 = bssolver(N,vb);
    vp = zeros(size(nm));
    dvp = zeros(size(nm));

    dN = 0;
    diffmu = 1;
    iter = 0;
    maxiter = 10;
    while (abs(diffmu)>eps && iter <= maxiter)

        vp = invert({solver,bssolver},nm,[mu,N],vi,vLi,vRi,eps,vp+dvp);

        [na,chia] = solver(mu,va+vp,vLi(1),vRi(1));
        [nb,chib] = bssolver(N,vb+vp);
        [eval,evec] = evsolver(N,vb+vp);
        fb = evec.^2;
        nf = na+nb;

        dvp = -(chib+chia)\fb;
        diffmu = mu-eval;

        deval_dN = (sum(dvp.*fb)*dx);
        dN = diffmu/deval_dN;

        dvp = dvp*dN; 

        N = N + dN;
        N = max(min(N,1),0);

        iter = iter + 1;
    end
    
    Nrec(i) = N;
    mub(i) = eval;
    vprec(:,i) = vp;
    
    narec(:,i) = na;
    nbrec(:,i) = nb;
    
end

% [n1,chi1] = solver(mu,v1+vp1,vLi(1),vRi(1));
% [n2,chi2] = bssolver(N,v2+vp1);

% 
% 
% subplot(2,2,1);
% plot(x,[na,nb,nm,nf]);
% xlim([min(x),max(x)]);
% 
% title('density');
% 
% subplot(2,2,3);
% plot(x,[va,vb,vp]);
% 
% ylim([min(va+vb)-.2,max(vp+.2)]);
% xlim([min(x),max(x)]);
% 
% title('potential');
% 
% subplot(2,2,4);
% plot(x,vp,'r');
% 
% ylim([min(vp)-.03,max(vp)+.03]);
% xlim([min(x),max(x)]);
% 
% title('vp');
% 
% subplot(2,2,2);
% 
% plot(x,nf-nm);
% xlim([min(x),max(x)]);
% 
% title('density error');