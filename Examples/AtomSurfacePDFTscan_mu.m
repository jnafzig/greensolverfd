clear;

L = 20; % Length of potential well
s = 5; % steepness parameter
V0 = 3.5; % Depth of potential well

R = 2;
Z = 3;

start = 7;
dx = 0.1; %spatial step size
xmin = -L-R;
xmax = start;

A = xmin;
B = xmax;
x = (xmin:dx:xmax)'; %initialize spatial grid
Nelem = max(size(x)); %Number of grid points

%Initialize potential of fragment 1
%va=-V0./(1+exp(-s*(x+R+L)))+V0./(1+exp(-s*(x+R)));
v_metal = -V0+V0./(1+exp(-s*(x+R)));

%initialize potential of fragment 2
v_atom = -Z*cosh(x).^-2;

%Create fragment potentials & left, right BCs on fragments
vi = [v_metal,v_atom];
vLi = [-V0,0];
vRi = [0,0];

solver = solver_fh(Nelem,dx);
shoot = shoot_fh(Nelem,dx);
bssolver = shootsolver_fh(Nelem,dx);
evsolver = shooteigsolver_fh(Nelem,dx);


% scan through a range of mu values
nummu = 100;
muspace = linspace(-2.5,-.35,nummu+1);
muspace = muspace(2:end);
dmu = muspace(2)-muspace(1);


Nrec = zeros(nummu,1);
Hrec = zeros(nummu,1);
mu_atom = zeros(nummu,1);
vprec = zeros(Nelem,nummu);
n_atomrec = zeros(Nelem,nummu);
n_metalrec = zeros(Nelem,nummu);
nm_rec = zeros(Nelem,nummu);

for i = 1:nummu

    mu = muspace(i);

    nm = solver(mu,v_metal+v_atom,sum(vLi),sum(vRi));

    n10 = solver(mu,v_metal,vLi(1),vRi(1));
    MaxEval = nodecount(shoot(0,v_atom,0,0));
    [evals,evecs] = evsolver(MaxEval,v_atom);
    diffs = evals - mu;
    [~,index]=min(abs(diffs));
    if index == 1
        n_core = zeros(Nelem,1);
    else
        n_core = sum(evecs(:,1:index-1).^2,2);
    end

    n_fukui = evecs(:,index).^2;
    N_atom0 = sum((nm-n10-n_core).*n_fukui)/sum(n_fukui.^2);
    N_atom0 = (index-1)+max(min(N_atom0,1),0);
    N_atom = N_atom0;

    n20 = bssolver(N_atom,v_atom);
    vp = zeros(size(nm));
    vp_pred = zeros(size(nm));
    dvp = zeros(size(nm));

    msk1 = nm<1e-8;
    msk2 = nm>1e-15;

    dN = 0;
    eval = 0;
    diffmu = 1;
    N_lowbound = -inf;
    N_upbound = inf;
    iter = 0;
    maxiter = 20;
    tol = eps;
    done = false;
    while (~done && iter <= maxiter)
        vp_old = vp;
        eval_old = eval;

        vp_pred(msk1) = vp_pred(find(~msk1,1,'last'));
    %     vp = invert({solver,bssolver},nm,[mu,Natom],vi,vLi,vRi,eps,vp_pred);

        vp_part = partialinvert({solver,bssolver},nm,[mu,N_atom],vi,vLi,vRi,msk2,eps,vp_pred);
        vp(msk2) = vp_part;
        vp(~msk2) = vp_pred(find(~msk1,1,'last'));

        [n_metal,chi_metal] = solver(mu,v_metal+vp,vLi(1),vRi(1));
        [n_atom,chi_atom] = bssolver(N_atom,v_atom+vp);
        nf = n_metal+n_atom;
        
        if mod(N_atom,1) ~= 0
            [eval,evec] = evsolver(N_atom,v_atom+vp);
            f_atom = evec(:,end).^2;
            mu_atom = eval(ceil(N_atom));
            mu_atom_plus = mu_atom;
            
            Epdft = sum(eval(1:floor(N_atom))) + eval(end)*mod(N_atom,1);
            E = Epdft - sum(vp.*n_atom)*dx;
        else
            [eval,evec] = evsolver(N_atom+1,v_atom+vp);
            
            mu_atom = eval(N_atom);
            mu_atom_plus = eval(N_atom+1);
            f_atom = evec(:,end-1).^2;
            
            if mu_atom < mu && mu < mu_atom_plus
                done = true;
            end
            
            Epdft = sum(eval(1:floor(N_atom)));
            E = Epdft - sum(vp.*n_atom)*dx;
        end
            
        % differential for how vp changes for small changes in Natom with
        % total density fixed
        dvp_dNatom = -(chi_atom+real(chi_metal))\f_atom; 
        % corresponding change in atom's chemical potential
        dmuatom_dN = (sum(dvp_dNatom.*f_atom)*dx);

        diffmu = mu-mu_atom;

        if (abs(diffmu) < tol)
            done = true;
        end

        dN = diffmu/dmuatom_dN; 

        % update bounds on N
        if dN > 0
            N_lowbound = N_atom;

            % should we worry about crossing integer or bound on N?
            if N_upbound <= floor(N_atom) + 1
               if dN + N_atom >= N_upbound
                   % then we bisect:
                   dN = (N_upbound + N_atom)/2 - N_atom;
               end
            else
                if dN + N_atom >= floor(N_atom) + 1
                   % then we go to upperbounding integer
                   dN = floor(N_atom) + 1 - N_atom;
                end
            end
        else
            N_upbound = N_atom;

            % should we worry about crossing integer or bound on N?
            if N_lowbound >= ceil(N_atom) - 1
               if dN + N_atom <= N_lowbound
                   % then we bisect:
                   dN = (N_lowbound + N_atom)/2 - N_atom;
               end
            else
                if dN + N_atom <= ceil(N_atom) - 1
                   % then we go to lowbounding integer
                   dN = ceil(N_atom) - 1 - N_atom;
                end
            end
        end


        % if upper and lower bounds are equal we are done
        if (N_lowbound == N_upbound) || done
            done = true;
        else
            N_atom = N_atom + dN;
        end

        dvp = dvp_dNatom*dN;
        vp_pred = vp + dvp;
        eval_pred = eval + dmuatom_dN*dN;

        iter = iter + 1;
    end
    
    Nrec(i) = N_atom;
    Hrec(i) = dmuatom_dN;
    mu_atomrec(i) = mu_atom;
    mu_atomrec_plus(i) = mu_atom_plus;
    vprec(:,i) = vp;
    
    n_metalrec(:,i) = n_metal;
    n_atomrec(:,i) = n_atom;
    nm_rec(:,i) = nm;
    
    Erec(i) = E;
    Epdftrec(i) = Epdft;

end