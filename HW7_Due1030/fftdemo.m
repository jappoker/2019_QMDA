% fftdemo
%
% Demonstrate how frequencies and coefficients are ordered in a fft

clearvars;  % clear variables

% parameters to play with
N=200;  % number of data points
deltat=1;  % sampling interval of d0 and d=d0+e
t=deltat*(0:N-1)';  % time coordinate
T0=deltat*N/4;  % period of sinusoidal d0
f0=1/T0;  % frequency of d0
phi=pi/4;  % phase of d0
sigmae=1;  % standard deviation of data errors e

% sinusoidal d0
d0=cos(2*pi*f0*t+phi);

% data d = d0 + data errors e, where e is normally distributed noise with a
% mean = zero and standard deviation = sigmae
e=sigmae*randn(N,1); 
d=d0+e;
d=d-mean(d);  % make sure d has zero mean

% plot d0 and d
figure(1);
plot(t,d0,'r--',t,d,'ko-','LineWidth',2);
title(sprintf('f_0 = %g, T_0 = %g',f0,T0));
legend('d_0','d');
xlabel('Time t');
ylabel('Data d_0 and d');
set(gca,'FontSize',18,'TickDir','out');

% compute vector of N +ve and -ve frequencies
fNyq=1/(2*deltat);  % Nyquist frequency
Npos=N/2+1;  % number of +ve frequencies
fpos=linspace(0,fNyq,Npos)';
fneg=flipud(-fpos(2:N/2));
f=[fpos; fneg];

% compute Fourier transform of d
ftd=fft(d);

% compute normalized periodogram pg (+ve and -ve frequencies)
pg=(1/N)*abs(ftd).^2;  % abs() returns the magnitude of complex numbers

% compute periodogram for positive frequencies only
pgpos=pg(1:N/2+1);  % take pg from f=0 to Nyquist
pgpos(2:N/2)=2*pgpos(2:N/2);  % add contribution from -ve frequencies

% plot periodogram
figure(2);
plot(fpos,pgpos,'bo-','LineWidth',2);
xlabel('Frequency f');
ylabel('Periodogram of d');
set(gca,'FontSize',18,'TickDir','out');


