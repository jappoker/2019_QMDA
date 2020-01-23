% initiate variables

clearvars;
load('QMDA_HW_09.mat');


% A.

d_ft = fft(d);
N = length(t);
Npos = N/2+1;
dt = t(2)- t(1);

fNyq=1/(2*dt);  % Nyquist frequency 
fpos=linspace(0,fNyq,Npos)';
fneg=flipud(-fpos(2:N/2));
freq=[fpos; fneg];

d_ft_pd = abs(d_ft).^2 ;

% figure for periodogram

figure(1);
subplot(2,1,1);
plot(fpos,d_ft_pd(1:Npos),'LineWidth',2);

ylabel('Periodogram');
xlabel('Frequency f');
set(gca,'FontSize',14);

subplot(2,1,2);
semilogy(fpos,d_ft_pd(1:Npos),'LineWidth',2);

ylabel('Periodogram in log space');
xlabel('Frequency f');
set(gca,'FontSize',14);


% B.

Ntaper = 5;
deltataper =  round (N /  ((Ntaper+1)/2)/2);


% figure 2 for preview of tapers
figure(2);
clf;

subplot(3,1,1);

title(sprintf('Numer of tapers N_w = %g',Ntaper));
ylabel('Tapering function');
xlabel('Time t');
set(gca,'FontSize',14);

hold on;

h = zeros(N,Ntaper); % call an empty list for hanntapers;
d_h  = zeros(N,Ntaper);
d_ft_han = zeros(N,Ntaper); % call an empty list for fourier transform of the data applied hanntapers;
d_pd_han = zeros(N,Ntaper); % call an empty list for periodogram for d_ft_han

for i = 1:1:Ntaper
    h(:,i) = hanntaper(t,t(deltataper * i),deltataper);
    
    % calculate the d_ft_han per loop
    d_h(:,i) = h(:,i) .* d;
    d_ft_han(:,i) = fft(d_h(:,i));
    d_pd_han(:,i) = abs(d_ft_han(:,i) ).^2;
    
    % plot the hanntaper per loop
    plot(t,h(:,i),'LineWidth',2);
    hold on;
    
end


d_pd_han_avg = mean(d_pd_han,2);

subplot(3,1,2);

plot(fpos,d_pd_han_avg(1:Npos),'LineWidth',2);

xlabel('Frequency f');
ylabel('Periodogram');
set(gca,'FontSize',14);


subplot(3,1,3);
semilogy(fpos,d_pd_han_avg(1:Npos),'LineWidth',2);

xlabel('Frequency f');
ylabel('Periodogram');
set(gca,'FontSize',14);

% B. several N w's

figure(3);
clf;

Ntaper_list = [3, 5, 13, 25,51];
Ntaper_number = length(Ntaper_list);

for n = 1:1:Ntaper_number
    Nt = Ntaper_list(n);
    deltataper =  round (N /  ((Nt+1)/2)/2);

    % figure 2 for preview of tapers
    


    h = zeros(N,Nt); % call an empty list for hanntapers;
    d_h  = zeros(N,Nt);
    d_ft_han = zeros(N,Nt); % call an empty list for fourier transform of the data applied hanntapers;
    d_pd_han = zeros(N,Nt); % call an empty list for periodogram for d_ft_han

    for i = 1:1:Nt
        h(:,i) = hanntaper(t,t(deltataper * i),deltataper);

        % calculate the d_ft_han per loop
        d_h(:,i) = h(:,i) .* d;
        d_ft_han(:,i) = fft(d_h(:,i));
        d_pd_han(:,i) = abs(d_ft_han(:,i) ).^2;
    end

    d_pd_han_avg = mean(d_pd_han,2);

    figure(3);

    bx(n)= subplot( 'Position',[0.1 0.1+((1-0.15)/Ntaper_number)*(Ntaper_number-n) 0.3 (1-0.22)/Ntaper_number]);
    plot(t,h,'LineWidth',2);
    ylabel(sprintf('N_w = %g',Nt));
    
    
    ax(n)=subplot( 'Position',[0.43 0.1+((1-0.15)/Ntaper_number)*(Ntaper_number-n) 0.42 (1-0.22)/Ntaper_number]);
    semilogx(1./fpos,d_pd_han_avg(1:Npos),'LineWidth',2);
    xlabel('Period (kyr) in log space');

    legend(sprintf('Numer of tapers N_w = %g',Nt));
    set(gca,'FontSize',14);
    
end

set(ax,'LineWidth',1,'FontSize',14,'YAxisLocation', 'right');
set(ax(1:Ntaper_number-1),'XTickLabel','');
set(bx,'LineWidth',1,'FontSize',14,'YTickLabel','');
set(bx(1:Ntaper_number-1),'XTickLabel','');

