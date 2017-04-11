load scan_center.mat;

pdf = dNatom_dmu;
cdf = Nrec;

muvals = muspace;

load scan_tails.mat;


msk = muspace < -0.46 & muspace > -.6;

muspace  = [muvals,muspace(msk)];

pdf = [pdf;dNatom_dmu(msk)];
cdf = [cdf;Nrec(msk)];

[muspace,ix] = sort(muspace);
pdf = pdf(ix);
cdf = cdf(ix);




logit = @(x) log(x)-log(1-x);

% mu = muspace(find(cdf == 1, 1, 'first'))-sum(cdf(cdf<1))*dmu;
% mu = sum(diff([cdf]).*(muspace(1:end-1)'+muspace(2:end)')/2);
% var = sum(diff([0;cdf]).*(muspace'-mu).^2);


mu = -0.484350161528228;
sigma =  0.000058395539667;
nu = 1.004848813108375;

% student = truncate(makedist('tlocationscale'),lb,ub);
student = makedist('tlocationscale');
student.mu = mu;
student.sigma =  sigma;
student.nu = nu;

logtfit = logtfit_fcn();

x = lsqcurvefit(logtfit,[student.mu,student.sigma,student.nu,0],...
    muspace',log(pdf),[-1,0,0,-.1],[0,1,5,.1]);
% x = [-0.484350161528228   0.000058395539667   1.004848813108375];

% student.mu = 0;
% alpha = .001;
% 
% plot(muspace,[log(pdf),log(2*student.pdf(muspace'-mu).*student.cdf(alpha*(muspace'-mu)))]);

plot(muspace,[log(pdf),logtfit(x,muspace')]);
% xlim([-.58,-.42]);

% plot(muspace,logit([cdf,student.cdf(muspace')]));
% xlim([-.485,-.4838]);

% sqrt(var)
% student.std
