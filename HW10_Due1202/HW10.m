% initiate variables

clearvars;
load('QMDA_HW_10.mat');



% Define the matrix
M = 500;
N = length(d);
m = [linspace(0,1,M+1)]';
m = m(1:end-1);
mean_prior = 0;
s_prior = 1;
l_corr = 0.05;
type_corr = 'G';

% sdeve = zeros(N,1);

sig_e = diag(sdeve.^2);


G = zeros(N,M);

for i = 1:1:N
    [r,c] = find(m==xd(i));
    G (i,r)= 1;
end

sig_prior = zeros(M,M);

for i = 1:1:M
    for j = 1:1:M
        sig_prior(i,j) = acvf(m(i)-m(j),s_prior,l_corr,type_corr);
    end
end

u_prior = mean_prior * ones(M,1);
d_prior = G * u_prior;
A = sig_prior * G';
B = G * A;

u_post = u_prior + A /(sig_e + B)*(d-d_prior);
sig_post = sig_prior - A /(sig_e + B) * A';
sig_post_sd = sqrt(diag(sig_post));


figure(1);
clf;

% Plot the raw data
errorbar(xd,d,sdeve,'LineWidth',2,'LineStyle','none','Marker','o');
set(gca,'FontSize',14,'LineWidth',1);
hold on;

% Plot the u_prior
plot(m,u_prior,"k--","LineWidth",2);

% Plot the interpolation result

plot(m,u_post,"LineWidth",2);
plot(m,u_post + sig_post_sd ,"LineWidth",2,"LineStyle",":");
plot(m,u_post - sig_post_sd ,"LineWidth",2,"LineStyle",":");


% Some plot settings
ylim([-1.5 1.5]);
title(sprintf('Interpolation with correlation length=%.2f, Function type of %s\n',l_corr,type_corr));
legs = legend('Data with error bar','$\mu_{prior}$', '$\hat{p}(r)$','$\hat{p}(r) + \sigma$','$\hat{p}(r)- \sigma$');
set(legs,'Interpreter','latex');
set(legs,'FontSize',14,'Location','south');