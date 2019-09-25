clearvars;

% read the data
D=csvread('GPAvsLSAT.csv');
x = D(:,1);
y = D(:,2);

% plot the data
figure(1);
plot(x,y,'ko','MarkerFaceColor','m','MarkerSize',10);
set(gca,'FontSize',16,'TickDir','out');

% Matrix 6
N = length(x);
G = ones(N,2);
G(:,2) = x;

% Guess for m
mguess=[0 mean(y./x)]';
ypred = G*mguess; % my predicted y

% plot the line
hold on;
plot(x,ypred,'b-');

hold off;

% error vector;
e = y-ypred;
etotal = e'*e;
title(sprintf('Total error = %g',etotal));