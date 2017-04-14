% Use the data from a run of AtomSurfacePDFTscan.m
AtomSurfacePDFTscan_mu

pdf = dNatom_dmu;
cdf = Nrec;

% initial guesses
mu = -0.48;
sigma =  5e-5;
nu = 1;
skew = 0;

% get the handle for the fit function
logtfit = logtfit_fcn();

% fit the curve
x = lsqcurvefit(logtfit,[mu,sigma,nu,skew],...
    muspace',log(pdf),[-1,0,0,-1],[0,1,5,1]);

plot(muspace,[log(pdf),logtfit(x,muspace')]);

mu = x(1)
sigma =  x(2)
nu = x(3)
skew = x(4)
