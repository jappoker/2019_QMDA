clearvars;  % clean up workspace


% plot part
% to better demonstrate the color gradient tendency, 
% d is chosen to be 0<= d<= 1, as p(1) = 5 exp(-5) or 
% p(1) = 10 exp (-10) are both approximate to 0

% parameter settings
d_min = 0;
d_max = 1;
n = 50; 
interval = d_max / n; 
d = [d_min:interval:d_max]';

% plot initialization
figure(2);
clf;
set(gca,'LineWidth',1,'FontSize',14);

axis equal;
set(gca,'XTick',[]); % turn off horizontal axis
set(gca,'YTick',[]); % turn off vertical axis

% ? borrowed from 'eda/ch03/eda03_06.m'

colormap(flipud(gray)*0.93);
hold on;

% ? color setting, with help of https://www.mathworks.com/help/matlab/ref/colormap.html


for l=[5 10]
    p = l*exp(-l*d);
    range = max(p)-min(p);
    imagesc( [3*l,3*l+3], [n , 0], (p-min(p))/range);
    text(3*l,-2,sprintf('\\lambda = %d',l));
    hold on;
end

colorbar;

clearvars; 

% calculation part
% to calculate it precious, as 0 <= d < infinity,
% d_max is set to be 10000

% parameter settings
d_min = 0;
d_max = 10000;
n = 1000000; 
interval = d_max / n; 
d = [d_min:interval:d_max]';



for l=[5 10]
    p = @(d) l*exp(-l*d);
    
    % mean = E(x) = \int x * p(x) dx
    
    p_d = @(d) l*exp(-l*d) .* d;
    meanval = integral( p_d,0,Inf);
    
    % variance = E(x^2) - E(x)^2
    %          = \int x^2 * p(x) dx - mean ^2
    
    p_d2 = @(d) l*exp(-l*d) .* d .*d;
    variance = integral(p_d2,0,Inf) - meanval^2;
    
    % median M: \int_0^M p(x) dx - 1/2 = 0
    
    intp = @(d) integral (p,0,d) - 0.5;
    medianval = fzero(intp,0);
    
    % mode is the d of p(d) = p_max
    
    [val,index] =  max(p(d));
    modeval = d(index);
    
    sprintf("when $\\lambda$ = %d, mean = %.3f, median = %.3f, variance = %.3f, mode = %.3f",l,meanval,medianval, variance, modeval)
end

