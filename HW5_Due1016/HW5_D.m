% initiate variables

clearvars;
load('QMDA_HW_05.mat');


close all;

% Guess

%%% INITIAL VALUES:
% sig_e = 0.0142;
% sig_m = 0.4;

% sig_e = 1;
% sig_m = 0.4;

sig_e = 0.1;
sig_m = 5;


% Calculation 

I = eye(N);
sig2_e = sig_e ^2;
Cov_e = sig2_e * I;
Mu_prior = zeros(N,1);

Cov_prior = sig_m ^2 * I;

Mu_post = Mu_prior + inv (inv( Cov_prior) +G' * inv(Cov_e) * G) * G' * inv(Cov_e) * (dobs - G * Mu_prior);
Cov_post = inv (inv( Cov_prior) + G' * inv(Cov_e) * G);
sig2_post = diag(Cov_post);
sig_post = sqrt(sig2_post);

% 1 figure
figure(1);
set(gca,'LineWidth',1,'FontSize',14);

subplot(2,1,1);
hold on;
plot(x,Mu_post,"-",'LineWidth',2,'MarkerSize',5);
plot(x,Mu_post + sig_post*2 ,":",'LineWidth',2,'MarkerSize',5);
plot(x,Mu_post - sig_post*2,":",'LineWidth',2,'MarkerSize',5);

legend({ sprintf('\\mu_{post}'),  sprintf('\\mu_{post} + 2\\sigma')  , sprintf('\\mu_{post} - 2\\sigma') });
xlabel('x');
title(sprintf('\\mu_{post} with guess \\sigma_e = %g & \\sigma_m = %g', sig_e, sig_m));
hold off;

% 2 figure

subplot(2,1,2);
hold on;

plot(x,dobs,"o",'LineWidth',2,'MarkerSize',5);
plot(x,G*Mu_post, "-",'LineWidth',2 );

legend({ sprintf('d_{obs}'), sprintf('d_pred = G\\mu_{post}')})
title(sprintf('Data and predictions with guess \\sigma_e = %g & \\sigma_m = %g', sig_e, sig_m));
