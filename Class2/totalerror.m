% totalerror
%
% Plot contour map of the total error in y = m1 + m2*x. x is LSAT and y GPA
% for the 15 law schools in the data of Diaconis and Efron (1983).

% clear all variables
clearvars;

% read the data
D=csvread('GPAvsLSAT.csv');
x=D(:,1);
y=D(:,2);

% plot the data
figure(1);
plot(x,y,'ko','MarkerFaceColor','m','MarkerSize',10);
set(gca,'FontSize',16,'TickDir','out');

% matrix G
N=length(x);
G=ones(N,2);
G(:,2)=x;

% guess for m
mguess=[0 mean(y./x)]';
ypred=G*mguess;

% plot the line for mguess
hold on;
plot(x,ypred,'b-');
hold off;

% error vector for mguess
e=y-ypred;
etotal=e'*e;
title(sprintf('Total error = %.4g at m_1 = %.4g, m_2 = %.4g',...
    etotal,mguess(1),mguess(2)));

% set up grid of total error
ngrid=100;
m1grid=linspace(-1.5,2,ngrid);
m2grid=linspace(0.2*mguess(2),1.5*mguess(2),ngrid);
etotalgrid=zeros(ngrid,ngrid);
for i=1:ngrid
    for j=1:ngrid
        m=[m1grid(j) m2grid(i)]';
        ypred=G*m;
        e=y-ypred;
        etotalgrid(i,j)=e'*e;
    end
end

% find minimum error and m values there
[etotalmin,jmin]=min(min(etotalgrid,[],1)); % min. of the min. of each column
[~,imin]=min(min(etotalgrid,[],2)); % min. of the min. of each row
m1min=m1grid(jmin);
m2min=m2grid(imin);

% plot contour map of total error and point with minimum error
figure(2);
levels=(0.5:0.5:5);
contour(m1grid,m2grid,etotalgrid,levels);  % contour map of etotal
hold on;
plot(m1min,m2min,'m*','MarkerSize',14);  % plot point where etotal is minimum
hold off;
axis xy;
colorbar;
title(sprintf('Min. total error = %.4g at m_1 = %.4g, m_2 = %.4g',....
    etotalmin,m1min,m2min));
xlabel('m_1');
ylabel('m_2');
ax=gca;
ax.YAxis.Exponent=0;  % avoid scientific notation on y axis
set(ax,'FontSize',16,'TickDir','out');

