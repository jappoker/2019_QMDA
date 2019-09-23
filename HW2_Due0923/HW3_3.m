clearvars;  % clean up workspace


% plot part
% to better demonstrate the color gradient tendency, 
% m_max is set to be 2, cuz p(2) appximates 0

% parameter settings
m_min = 0;
m_max = 2;
n = 50; 
interval = m_max / n; 
m = [m_min:interval:m_max]';

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


p = 2* 2.*m/(sqrt(2*pi)) .* exp(- m.^4 /2) ;
range = max(p)-min(p);
imagesc( [0,3], [n , 0], (p-min(p))/range);
colorbar;
hold on;


% Calculate part

clearvars;
p = @(m) 2*2.*m/(sqrt(2*pi)) .* exp(- m.^4 /2) ;
    
% mean = E(x) = \int x * p(x) dx

p_m = @(m) 2*2.*m.^2/(sqrt(2*pi)) .* exp(- m.^4 /2);
meanval = integral( p_m,0,Inf);

% variance = E(x^2) - E(x)^2
%          = \int x^2 * p(x) dx - mean ^2

p_m2 =@(m) 2*2.*m.^3/(sqrt(2*pi)) .* exp(- m.^4 /2);
variance = integral(p_m2,0,Inf) - meanval^2;

% median M: \int_0^M p(x) dx - 1/2 * \int_0^Inf p(x) dx = 0

intp = @(m) 2*integral(p,0,m) - integral(p,0,Inf);
medianval = fzero(intp, [0,10]);

% mode is the d of p(d) = p_max

m_min = 0;
m_max = 2;
n = 5000; 
interval = m_max / n; 
m = [m_min:interval:m_max]';

[val,index] =  max(p(m));
modeval = m(index);

sprintf("mean = %.3f, median = %.3f, variance = %.3f, mode = %.3f",meanval,medianval, variance, modeval)
