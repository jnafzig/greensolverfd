clear;

A = -15;
B = 8;
Nelem = 150;
x = linspace(A,B,Nelem)';
dx = x(2)-x(1);

vLi = [-1,0];
vRi = [0,0];
R = 5;
vmetal = zeros(Nelem,1);
vmetal(x<0) = vLi(1);
scale = 4;
vatom = -(1+scale)*cosh((x-R)*(1+2*scale)).^-2;
vi = [vmetal,vatom];

% grab the solver functions
solver = solver_fh(Nelem,dx);
shoot = shoot_fh(Nelem,dx);
bssolver = shootsolver_fh(Nelem,dx);
evsolver = shooteigsolver_fh(Nelem,dx);

% scan through a range of mu values
nummu = 10;
muspace = linspace(-.55,-.45,nummu+1);
muspace = muspace(2:end);
dmu = muspace(2)-muspace(1);

Nrec = zeros(nummu,1);
Hrec = zeros(nummu,1);
dNatom_dmu = zeros(nummu,1);
dNmetal_dmu = zeros(nummu,1);
dNtot_dmu = zeros(nummu,1);
mu_atom = zeros(nummu,1);
vprec = zeros(Nelem,nummu);
narec = zeros(Nelem,nummu);
nbrec = zeros(Nelem,nummu);

for i = 1:nummu

    mu = muspace(i);
    nm = solver(mu,vmetal+vatom,sum(vLi),sum(vRi));

    n10 = solver(mu,vmetal,vLi(1),vRi(1));
    n20 = bssolver(1,vatom);
    Natom = sum((nm-n10).*n20)/sum(n20.^2);
    Natom = max(min(Natom,1),0);
    
    n20 = bssolver(Natom,vatom);
    vp = zeros(size(nm));
    dvp = zeros(size(nm));

    dN = 0;
    eval = 0;
    diffmu = 1;
    iter = 0;
    tol = eps;
    done = false;
    maxiter = 10;
    while (~done && iter <= maxiter)
        vp_old = vp;
        eval_old = eval;
        
        vp = invert({solver,bssolver},nm,[mu,Natom],vi,vLi,vRi,eps,vp+dvp);

        [na,chia] = solver(mu,vmetal+vp,vLi(1),vRi(1));
        [nb,chib] = bssolver(Natom,vatom+vp);
        [eval,evec] = evsolver(1,vatom+vp);
        fb = evec.^2;
        nf = na+nb;
        
        % differential for how vp changes for small changes in Natom with
        % total density fixed
        dvp_dNatom = -(chib+chia)\fb;
        % corresponding change in atom's chemical potential
        dmuatom_dN = (sum(dvp_dNatom.*fb)*dx);
        
        diffmu = mu-eval;

        if (abs(diffmu) < tol)
            done = true;
        end
        
        dN = diffmu/dmuatom_dN;
        dN = min(max(-Natom,dN),1-Natom);  % bound dN so N stays between 0 and 1
        
        % if N is at the edge of the interval and doesn't go anywhere we are
        % done
        if ((Natom == 1 || Natom == 0) && dN == 0)
            done = true;
        else
            Natom = Natom + dN;
        end
        
        dvp = dvp_dNatom*dN; 
        vp_pred = vp + dvp;
        eval_pred = eval + dmuatom_dN*dN;

        iter = iter + 1;
    end

    % changes in densities w.r.t changes in total chemical potential
    dnm_dmu = -2*real(ldos(shoot(mu,vmetal+vatom,sum(vLi),sum(vRi))));
    dna_dmu = -2*real(ldos(shoot(mu,vmetal+vp,sum(vLi),sum(vRi))));
    dNatom_dmu(i) = -(1-sum(fb.*((chia+chib)\(dnm_dmu-dna_dmu)))*dx)./(sum(fb.*((chia+chib)\fb))*dx);
    dNmetal_dmu(i) = sum(dna_dmu)*dx;
    dNtot_dmu(i) = sum(dnm_dmu)*dx;
    
    Nrec(i) = Natom;
    Hrec(i) = dmuatom_dN;
    mu_atom(i) = eval;
    vprec(:,i) = vp;
    
    narec(:,i) = na;
    nbrec(:,i) = nb;
    
end

