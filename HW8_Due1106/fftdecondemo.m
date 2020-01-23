% fftdecondemo
%
% Demonstrate convolution and deconvolution with fft

clearvars;  % clear variables

% parameters to play with
N=1500;  % number of data points
deltat=1;  % sampling interval of data and error vectors
t=deltat*(0:N-1)';  % time coordinate
tg=100;  % characteristic time of instrument impulse response
sigmae=0.2;  % standard deviation of data errors e
t0ctrue=200;  % start and end time of pulse in ctrue
t1ctrue=700;

% plotting parameters
fontsize=18;
linewidth=3;

% filter impulse response g
g=exp(-t/tg);
g=g/sum(g);  % normalize so that sum(g)=1

% boxcar ctrue
ctrue=zeros(N,1);
ipulse=find((t>=t0ctrue)&(t<=t1ctrue));
ctrue(ipulse)=1;

% obtain "clean" measurements dclean by frequency domain convolution
fftg=fft(g);
fftctrue=fft(ctrue);
fftdclean=fftg.*fftctrue;
dclean=ifft(fftdclean);

% measured data dmeas = dclean + noise e
e=sigmae*randn(N,1);
dmeas=dclean+e;

% plot ctrue, filter, dclean
figure(1);
subplot(3,1,1);
plot(t,g,'r-','LineWidth',linewidth);
ylabel('Filter g');
set(gca,'FontSize',fontsize,'TickDir','out');
subplot(3,1,2);
plot(t,ctrue,'r-','LineWidth',linewidth);
ylabel('Concentration ctrue');
set(gca,'FontSize',fontsize,'TickDir','out');
subplot(3,1,3);
plot(t,dmeas,'r-',t,dclean,'b-','LineWidth',linewidth);
legend('dmeas','dclean');
ylabel('Measurements');
xlabel('Time t');
set(gca,'FontSize',fontsize,'TickDir','out');

% try naive deconvolution
fftdmeas=fft(dmeas);
fftcestnaive=fftdmeas./fftg;
cestnaive=ifft(fftcestnaive);
% plot
figure(2);
plot(t,cestnaive,'c-',t,ctrue,'r--','LineWidth',linewidth);
legend('cestnaive','ctrue');
ylabel('Estimated conc.');
xlabel('Time t');
title('Naive deconvolution');
set(gca,'FontSize',fontsize,'TickDir','out');

% water level deconvolution
eps2=0.01;  % water level
fftnumer=fftdmeas.*conj(fftg);
fftdenom=fftg.*conj(fftg)+eps2;
fftcestwl=fftnumer./fftdenom;
cestwl=ifft(fftcestwl);
% plot
figure(3);
plot(t,cestwl,'b-',t,ctrue,'r--','LineWidth',linewidth);
legend('cestwl','ctrue');
ylabel('Estimated conc.');
xlabel('Time t');
title('Water level deconvolution');
set(gca,'FontSize',fontsize,'TickDir','out');

% vector of +ve frequencies
fNyq=1/(2*deltat);  % Nyquist frequency
Npos=N/2+1;  % number of +ve frequencies
fpos=linspace(0,fNyq,Npos)';
fneg=flipud(-fpos(2:N/2));
f=[fpos; fneg];

% periodogram of measured data (+ve frequencies only!)
pgramdmeas=abs(fftdmeas(1:Npos)).^2;
figure(4);
subplot(2,1,1);
semilogy(fpos,pgramdmeas,'r-','LineWidth',linewidth);
ylabel('Periodogram of dmeas');
xlabel('Frequency f');
set(gca,'FontSize',fontsize,'TickDir','out');

% Wiener filter phi = Gaussian lowpass filter
f0=0.005;  % controls width of frequency response
phi=exp(-(abs(f)/f0).^2);
% plot
subplot(2,1,2);
plot(fpos,phi(1:Npos),'m-','LineWidth',linewidth);
ylabel('Wiener filter \phi');
xlabel('Frequency f');
set(gca,'FontSize',fontsize,'TickDir','out');

% Wiener deconvolution
fftcestwiener=(fftdmeas.*phi)./fftg;
cestwiener=ifft(fftcestwiener);
% plot
figure(5);
plot(t,cestwiener,'k-',t,ctrue,'r--','LineWidth',linewidth);
legend('cestwiener','ctrue');
ylabel('Estimated conc.');
xlabel('Time t');
title('Wiener deconvolution');
set(gca,'FontSize',fontsize,'TickDir','out');
