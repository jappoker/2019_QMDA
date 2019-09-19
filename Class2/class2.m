% mcestimvar
%
% Monte Carlo calculation of variance estimates

clearvars;  % clean up workspace

n=10;  % number of samples
nmc=50000;  % number of Monte Carlo iterations
nbins=40;  % number of bins in histogram

% compute nmc estimates of variance for n smaples
varmc0=zeros(nmc,1);
varmc1=zeros(nmc,1);
for i=1:nmc
    x=randn(n,1);
    meanx=mean(x);
    % meanx=0;
    sumsq=sum((x-meanx).^2);
    varmc0(i)=sumsq/n;
    varmc1(i)=sumsq/(n-1);
end

% plot
figure(1);
subplot(1,2,1);
histogram(varmc0,nbins);
title(sprintf('Mean est. \\sigma^2 = %.3f (SSQ/N)',mean(varmc0)));
set(gca,'FontSize',14,'TickDir','out');
subplot(1,2,2);
histogram(varmc1,nbins);
title(sprintf('Mean est. \\sigma^2 = %.3f (SSQ/(N-1))',mean(varmc1)));
set(gca,'FontSize',14,'TickDir','out');
