clearvars; 

% plot initialization
figure(2);
clf;
set(gca,'LineWidth',1,'FontSize',14);
hold on;

% parameter settings
l_list = 1:1:10;

meanlist = zeros(10,1);
varlist = zeros(10,1);

for l=l_list
    p = @(d) l*exp(-l*d);
    
    % mean = E(x) = \int x * p(x) dx
    
    p_d = @(d)  p(d).*d;
    meanlist(l) = integral( p_d,0,Inf);
    
    % variance = E(x^2) - E(x)^2
    %          = \int x^2 * p(x) dx - mean ^2
    
    p_d2 = @(d)  p(d).*d.^2;
    varlist(l) = integral(p_d2,0,Inf) - meanlist(l)^2;
    
end

plot(l_list,meanlist,"-o",l_list,varlist,"-o",...
      'LineWidth',2, 'MarkerSize',10);
legend('Mean','Variance');

xlabel(sprintf('\\lambda'));
ylabel('Mean & Variance value');

