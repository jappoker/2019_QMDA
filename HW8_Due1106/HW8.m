% initiate variables

clearvars;
f=csvread('NYC_temp.csv',1,2);
t=f(:,1);
tmax=f(:,2);
tmin=f(:,3);

% Color setting:
Tmaxcolor = '#ef9a9a';
Tmincolor = '#81d4fa';

% A.

% creat new t
years=[1995:1:2014];
t_new=[];

leapyear=linspace(0,1,367);
leapyear=leapyear(1:366)';
normalyear=linspace(0,1,366);
normalyear=normalyear(1:365)';

for i = years
    if mod(i,4) == 0
        t_new = [t_new; i+leapyear];
    else 
        t_new = [t_new; i+normalyear];
    end
end

% plot
figure(1);
hold on;
plot(t_new,tmin, 'Marker','.','LineStyle','none','Color',Tmincolor)
plot(t_new,tmax,'Marker','.','LineStyle','none','Color',Tmaxcolor);
xlabel('Time (year)');
ylabel('degrees Celsius times 10');
legend('T_{min}','T_{max}')
set(gca,'LineWidth',1,'FontSize',14);


% B.

% compute fourier tranform of Tmax and Tmin
tmax_ft=fft(tmax);
tmax_pd=abs(tmax_ft).^2;
tmin_ft=fft(tmin);
tmin_pd=abs(tmin_ft).^2;

% compute vector of N +ve and -ve frequencies
N = length(t);
Nf = (N+1)/2;
dti = 1/365;
fNyq=1/(2*dti);  % Nyquist frequency 
fpos=linspace(0,fNyq,Nf)';
fneg=flipud(-fpos(2:Nf));
freq=[fpos; fneg];

% periodogram  (+ve frequencies only!)
figure(2);

subplot(3,1,1);
semilogy(fpos,tmax_pd(1:Nf),'LineWidth',2,'Color',Tmaxcolor);
ylabel('Periodogram of T_{max}');
xlabel('Frequency f');
set(gca,'FontSize',14);

subplot(3,1,2);
semilogy(fpos,tmin_pd(1:Nf),'LineWidth',2,'Color',Tmincolor);
ylabel('Periodogram of T_{min}');
xlabel('Frequency f');
set(gca,'FontSize',14);

% Wiener filter phi = Gaussian lowpass filter
f0=3;  % controls width of frequency response
phi=exp(-(abs(freq)/f0).^2);
% plot
subplot(3,1,3);
plot(fpos,phi(1:Nf),'m-','LineWidth',2);
ylabel('Wiener filter \phi');
xlabel('Frequency f');
set(gca,'FontSize',14);

% Wiener deconvolution to get smoothed version of Tmin and Tmax
tmax_filtered=(tmax_ft.*phi);
tmin_filtered=(tmin_ft.*phi);
tmax_filtered_ift=ifft(tmax_filtered);
tmin_filtered_ift=ifft(tmin_filtered);

% D.
% plot
figure(3);

subplot(2,1,1);
hold on;
plot(t_new,tmax,'Marker','.','LineStyle','none','Color',Tmaxcolor);
plot(t_new,tmax_filtered_ift,'k-','LineWidth',3);

legend('True T_{max}','smoothen T_{max}');
ylabel('T_{max}');
xlabel('Time t');
set(gca,'FontSize',14);

subplot(2,1,2);
hold on;
plot(t_new,tmin, 'Marker','.','LineStyle','none','Color',Tmincolor)
plot(t_new,tmin_filtered_ift,'k-','LineWidth',3);

legend('True T_{min}','smoothen T_{min}');
ylabel('T_{min}');
xlabel('Time t');
set(gca,'FontSize',14);

% E.
% obtain the difference between Tmin and Tmax
dt = tmax-tmin;
dt_filtered = tmax_filtered_ift - tmin_filtered_ift;


figure(4);

plot(t_new, dt,'.', t_new,dt_filtered,'k-','LineWidth',3);
for i = years
    xline(i,'k');
    xline(i+0.5,'r');
end

legend('True T_{difference}','smoothen T_{difference}');
ylabel('Temperature differrence');
xlabel('Time t');
set(gca,'FontSize',14);

