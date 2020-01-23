% initiate variables

clearvars;
load('QMDA_HW_07.mat');


% A.

figure(1);
plot(ti,hi, 'LineWidth',2);
xlabel('Time (Day)');
ylabel('Sea surface elevation (m)');
set(gca,'LineWidth',1,'FontSize',14);

% B.

% To calculte the delta_t of sampling, we calculate the differenciate of ti

dtilist = diff(ti);

% to check if the time intervals are equal
dti = mean(dtilist);
fprintf('The mean of the time interval %g, and the max and min of the dti are %g, %g\n',dti, max(dtilist),min(dtilist));

% C.

% compute fourier coefficients
hft = fft(hi);

% compute vector of N +ve and -ve frequencies
N = length(ti);
Nf = N/2+1;
fNyq=1/(2*dti);  % Nyquist frequency 
fpos=linspace(0,fNyq,Nf)';
fneg=flipud(-fpos(2:N/2));
f=[fpos; fneg];


% D.
% compute periodogram
pg = abs(hft).^2;

N1 = 2;
N15 = round(15/(fNyq/Nf));
% plot the periodogram
figure(2);
% semilogy - https://www.mathworks.com/help/matlab/ref/semilogy.html
semilogy(f(N1:N15),pg(N1:N15),'LineWidth',2);

title('Periodogram');
xlabel('Frequency cycles/day');
ylabel('Squared magnitude of hFFT');
set(gca,'LineWidth',1,'FontSize',14);

% E.



