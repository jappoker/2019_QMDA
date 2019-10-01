clearvars;

% HW 3 - A
% read the data
D=csvread('GPAvsLSAT.csv');
x = D(:,1);
y = D(:,2);

% plot the data
figure(1);
plot(x,y,'ko','MarkerFaceColor','b','MarkerSize',10);
set(gca,'FontSize',16,'TickDir','out');
xlabel('LSAT'); 
ylabel('GPA');


% HW 3 - B
% corrcoef: https://www.mathworks.com/help/matlab/ref/corrcoef.html

[R,P,RL,RU] = corrcoef(x,y);

r = R(1,2);
p = P(1,2);

% fprintf('r = %g, p = %g\n',r,p);
 
% HW 3 - C
% x = [1 2 3 4 5];
% y = [10 20 30 40 50];
% 
% [X,Y] = samplebootstrap(x,y);


% HW 3 - C & D

N = 1000;
rlist = zeros(N,1);
plist = rlist;
for i=1:N
    [X,Y] = samplebootstrap(x,y);
    [R,P,RL,RU] = corrcoef(X,Y);
    rlist(i) = R(1,2);
    plist(i) = P(1,2);
end

% plot r
figure(2);
histogram(rlist,40);
title("r");
set(gca,'FontSize',16,'TickDir','out');
xlabel("Correlation coefficient");
ylabel("# samples per interval of width");


% HW 3 - E

close all;
CI = 0.95;

% 2.5 % * 1000 = 25, so we will have the 26 to 975 values

new_r = sort(rlist);
a = round(N * (1-CI)/2) ;
b = N - a;
new_r = new_r(a + 1:b);

% fprintf("The 0.95 confidence interval is (%.3f, %.3f)\n",new_r(1), new_r(end));

% HW 3 - F

n_check = 5;
CI = 0.95;
N = 1000;       % initial N = 1000

syms RL_list RU_list;

while  1 
    RL_list = zeros(n_check,1);
    RU_list = RL_list;
    
    for c = 1:n_check
        rlist = zeros (N,1);
        
        for i=1:N
            [X,Y] = samplebootstrap(x,y);
            [R,P,RL,RU] = corrcoef(X,Y);
            rlist(i) = R(1,2);
        end
        
        new_r = sort(rlist);
        a = round(N * (1-CI)/2) ;
        b = N - a;
        new_r = new_r(a + 1:b);
        
        RL_list(c) = round(new_r(1),2);
        RU_list(c) = round(new_r(end),2);
    end
    
    if  max(RL_list) == min(RL_list) && max(RU_list) == min(RU_list)
        % make sure they give same two digits and jump out of the while
        break;
    else
        N = N + 1000;
    end
end

fprintf("The total sampling number to get the 2 digits same for repeating %g times: %g\n",n_check, N);
fprintf("The 0.95 confidence interval is (%g, %g)\n", RL_list(1), RU_list(1));