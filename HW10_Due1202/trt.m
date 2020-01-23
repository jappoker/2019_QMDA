clearvars;
% Read in variables: d (data values); xd (x-coordinates of data points);
% sdeve (standard deviation of the measurement errors)
filename = 'QMDA_HW_10.mat';
load(filename) 

N = length(d);
%sdeve = zeros(N,1); % set all sdeve_prior to be zero
%sdeve = diag(sdeve.^2);
corrlength = 0.1; % correlation length
functype = 'G'; % Gaussian autocovariance function

% Calculate the generalized least sqaures interpolation
% prior information
M = 500;
xm = linspace(0,1,M+1)'; % grid x coordinates
xm = xm(1:M);
mean_prior = 0;
covar_e = diag(sdeve.^2); % since mean_prior is zero
vari_prior = 1; % prior variance
% sigma_prior
covar_prior = zeros(M,M);
for i=(1:M)
    for j=(1:M)
        covar_prior(i,j) = acvf(xm(i)-xm(j), vari_prior, corrlength, functype);
    end 
end

% Construct matrix A (MxN)
A = zeros(M,N);
for i=(1:M)
    for j=(1:N)
        A(i,j) = acvf(xm(i)-xd(j), vari_prior, corrlength, functype);
    end
end

% Construct matrix B (NxN) B=G*A
B = zeros(N,N);
for i=(1:N)
    for j=(1:N)
        B(i,j) = acvf(xd(i)-xd(j), vari_prior, corrlength, functype);
    end 
end

% d_prior
%d_prior = G*mean_prior;

% Compute posterior 
covar_post = covar_prior-A*((covar_e+B)\(A'));  % Alberto's change
sdeve_post2 = diag(covar_post);
sdeve_post = sqrt(sdeve_post2);

mean_post = A*((covar_e + B)\d); % for sigma_e being all zeros
%mean_post = mean_prior+A*inv(covar_e+B)*(d);

% Plot interpolated posterior mean, its value ± one posterior standard 
% deviation, the data points and their error bars as a function of x
bluecolor = [0, 0.4470, 0.7410];
redcolor = [0.6350, 0.0780, 0.1840];

figure(5)
hold on;
% title(sprintf('Interpolation by generalized least squares with corrlength=%.2f, Function type is %s\n',corrlength,functype));
% prior
plot(xm,mean_prior*ones(length(xm),1),'k--',"LineWidth",1.5); % mean_prior
errorbar(xd,d,sdeve,'ko','MarkerSize',5,...
    'MarkerEdgeColor','black','MarkerFaceColor','black',"LineWidth",1)
% posterior
plot(xm,mean_post,'-',"LineWidth",1.5,'color',bluecolor);
plot(xm,mean_post+sdeve_post,'--',"LineWidth",1.5,'color',redcolor);
plot(xm,mean_post-sdeve_post,'--',"LineWidth",1.5,'color',bluecolor);
xlim([0 1]);
ylim([-1 1]);
hold off;