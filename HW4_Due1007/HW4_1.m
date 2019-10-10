n = 40;
l = [1 zeros(1,n-1)]; 
c = [1 1 zeros(1,n-2)];
G = toeplitz(c,l); 

sig_d = 1;
I = eye(n,n);

covar_d = sig_d ^ 2 * I;
M =  (G' * G) \ G';
covar_m = M * covar_d * M';

sig_m = sqrt(diag(covar_m));

figure(1);

plot(1:n,sig_m,"LineWidth",2);
xlabel("j");
ylabel(sprintf("\\sigma_{mj}"));
set(gca,'FontSize',16,'LineWidth',1);