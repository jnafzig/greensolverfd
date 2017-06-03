clear;

% mu = -2.019;
% mu = -2.010;
% mu = -2.003;
mu = -.5;

L = 20; % Length of potential well
s = 5; % steepness parameter
V0 = 3.5; % Depth of potential well

R = 3;
Z = 3;

start = 10;
dx = 0.1; %spatial step size
xmin = -L-R;
xmax = start;

A = xmin;
B = xmax;
x = (xmin:dx:xmax)'; %initialize spatial grid
Nelem = max(size(x)); %Number of grid points

%Initialize potential of fragment 1
%va=-V0./(1+exp(-s*(x+R+L)))+V0./(1+exp(-s*(x+R)));
va = -V0+V0./(1+exp(-s*(x+R)));

%initialize potential of fragment 2
vb = -Z*cosh(x).^-2;

%Create fragment potentials & left, right BCs on fragments
vi = [va,vb];
vLi = [-V0,0];
vRi = [0,0];

solver = solver_fh(Nelem,dx);
shoot = shoot_fh(Nelem,dx);
bssolver = shootsolver_fh(Nelem,dx);
evsolver = shooteigsolver_fh(Nelem,dx);

nm = solver(mu,va+vb,sum(vLi),sum(vRi));

n10 = solver(mu,va,vLi(1),vRi(1));
MaxEval = nodecount(shoot(0,vb,0,0));
[evals,evecs] = evsolver(MaxEval,vb);
diffs = evals - mu;
[~,index]=min(abs(diffs));
if index == 1
    ncore = zeros(Nelem,1);
else
    ncore = sum(evecs(:,1:index-1).^2,2);
end

nfuk = evecs(:,index).^2;
Natom = sum((nm-n10-ncore).*nfuk)/sum(nfuk.^2);
Natom = (index-1)+max(min(Natom,1),0);

n20 = bssolver(Natom,vb);
vp = zeros(size(nm));
vp_pred = zeros(size(nm));
dvp = zeros(size(nm));

msk1 = nm<1e-10;
msk2 = nm>1e-15;

dN = 0;
eval = 0;
diffmu = 1;
iter = 0;
maxiter = 10;
tol = eps;
done = false;
while (~done && iter <= maxiter)
    vp_old = vp;
    eval_old = eval;

    vp_pred(msk1) = vp_pred(find(~msk1,1,'last'));
%     vp = invert({solver,bssolver},nm,[mu,Natom],vi,vLi,vRi,eps,vp_pred);

    vp_part = partialinvert({solver,bssolver},nm,[mu,Natom],vi,vLi,vRi,msk2,eps,vp_pred);
    vp(msk2) = vp_part;
    vp(~msk2) = vp_pred(find(~msk1,1,'last'));
    
    [na,chia] = solver(mu,va+vp,vLi(1),vRi(1));
    [nb,chib] = bssolver(Natom,vb+vp);
    [eval,evec] = evsolver(Natom,vb+vp);
    fb = evec(:,end).^2;
    nf = na+nb;
 
    % differential for how vp changes for small changes in Natom with
    % total density fixed
    dvp_dNatom = -(chib+real(chia))\fb; 
    % corresponding change in atom's chemical potential
    dmuatom_dN = (sum(dvp_dNatom.*fb)*dx);

    diffmu = mu-eval(ceil(Natom));

    if (abs(diffmu) < tol)
        done = true;
    end

    dN = diffmu/dmuatom_dN;

    if mod(Natom,1)==0
        moddiff = 0;
    else
        moddiff=1-mod(Natom,1);
    end
    dN = min(max(-Natom,dN),moddiff);  % bound dN so N stays between integer numbers

    % if N is at the edge of the interval and doesn't go anywhere we are
    % done
    if ((mod(Natom,1)==0) && dN == 0)
        done = true;
    else
        Natom = Natom + dN;
    end

    dvp = dvp_dNatom*dN;
    vp_pred = vp + dvp;
    eval_pred = eval + dmuatom_dN*dN;

    iter = iter + 1;
end


