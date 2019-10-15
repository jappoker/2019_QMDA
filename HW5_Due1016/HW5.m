% initiate variables

clearvars;
load('QMDA_HW_05.mat');

% A. plot data kernal G & observed data a a function of x
%    plot initialization
figure(1);
clf;
axis ij;
set(gca,'LineWidth',1,'FontSize',14);
colormap(flipud(gray)*0.93);
hold on;

%    plot kernal G
imagesc(G);
colorbar;

%    figure settings
title('Kernal G');
xlim([1 M]);
ylim([1 N]);
ylabel("N");
xlabel("M");

%    plot d obs
figure(2);
set(gca,'LineWidth',1,'FontSize',14);
hold on;

plot(x,dobs,"ko-",'LineWidth',2,'MarkerSize',5);

%    figure settings
title('Observed data vs x');
ylabel(sprintf('d_{obs}'));
xlabel("x");

% B.

mest = (G'*G)\(G'*dobs);
% plot(x,G*mest,'LineWidth',4);
e = dobs - G * mest;
E = e'*e;
sig2_e = E/(N-1);
sig_e = sqrt(sig2_e);

fprintf('The standard deviation of the data error \\sigma_e = %g\n',sig_e);

% C.

close all;

I = eye(N);
Cov_e = sig2_e * I;
Mu_prior = zeros(N,1);

% Guess

sig_m = 0.4;
Cov_prior = sig_m ^2 * I;

% Calculation 

Mu_post = Mu_prior + inv (inv( Cov_prior) +G' * inv(Cov_e) * G) * G' * inv(Cov_e) * (dobs - G * Mu_prior);
Cov_post = inv (inv( Cov_prior) + G' * inv(Cov_e) * G);
sig2_post = diag(Cov_post);
sig_post = sqrt(sig2_post);

% C.1 figure
figure(3);
set(gca,'LineWidth',1,'FontSize',14);
hold on;

plot(x,Mu_post,"-",'LineWidth',2,'MarkerSize',5);
plot(x,Mu_post + sig_post*2 ,":",'LineWidth',2,'MarkerSize',5);
plot(x,Mu_post - sig_post*2,":",'LineWidth',2,'MarkerSize',5);

legend({ sprintf('\\mu_{post}'),  sprintf('\\mu_{post} + 2\\sigma')  , sprintf('\\mu_{post} - 2\\sigma') });
xlabel('x');
title(sprintf('\\mu_{post} with guess \\sigma_m = %g', sig_m));

% C.2 figure
figure(4);
set(gca,'LineWidth',1,'FontSize',14);
hold on;

plot(x,dobs,"o",'LineWidth',2,'MarkerSize',5);
plot(x,G*Mu_post, "-",'LineWidth',2 );

legend({ sprintf('d_{obs}'), sprintf('d_pred = G\\mu_{post}')})
title(sprintf('Data and predictions with guess \\sigma_m = %g', sig_m));


