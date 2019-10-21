% initiate variables

clearvars;
load('MOR_EOF_228_Z.mat');


% A. 

[U,S,V] = svd(Z,'econ');
x=linspace(-40,40,81);

figure(1);

for i = 1:5
    ax(i)=subplot( 'Position',[0.1 0.1+0.17*(5-i) 0.8 0.17]);
    
    % To match the paper Small 1998, the Mode 1, 2, 5 are -V(:i), and the
    % Mode 3 and 4 are V(:i), as mentioned in the problem instructions.
    
    if i == 1 || i ==2 || i == 5
        plot(x,-V(:,i),"k-","LineWidth",2);
    else
        plot(x,V(:,i),"k-","LineWidth",2);
    end
    
    % Make the ylim looks better:
    
    lim = max(abs(V(:,i)))/0.8;
    ylim([-lim,lim])
    legend(sprintf('Mode %g',i));

end

set(ax,'LineWidth',1,'FontSize',14,'YTickLabel','');
set(ax(1:4),'XTickLabel','');
xlabel('Distance(km)');


% B.

var = diag(S);
var_ttl = sum(var);
var_percentage = var/var_ttl;


figure(2);
hold on;

plot(1:length(var_percentage), var_percentage, "ko-","MarkerSize",8,"LineWidth",2);
plot(1:5, var_percentage(1:5),"ko-","MarkerSize",8,"LineWidth",2,'MarkerFaceColor',"#000000");
text(20,0.1,sprintf('First 5 Sum is: %.3f',sum(var_percentage(1:5))),'FontSize',16 );
set(gca,'LineWidth',1,'FontSize',14);
xlabel('Mode');
ylabel('% Variance');
hold off;

% C.

S2 = S^2;

var2 = diag(S2);
var2_ttl = sum(var2);
var2_percentage = var2/var2_ttl;

figure(3);
hold on;

plot(1:length(var2_percentage), var2_percentage, "ko-","MarkerSize",8,"LineWidth",2);
plot(1:5, var2_percentage(1:5),"ko-","MarkerSize",8,"LineWidth",2,'MarkerFaceColor',"#000000");
text(20,0.2,sprintf('First 5 Sum is: %.3f',sum(var2_percentage(1:5))),'FontSize',16 );
set(gca,'LineWidth',1,'FontSize',14);
xlabel('Mode');
ylabel('% Variance');

hold off;



% D.
% close all;

Sp = zeros(length(S),length(S));
Sp(1:5,1:5) = S(1:5,1:5);

Zp = U*Sp*V';

rows = [1 80 120];

figure(4);
for r = rows
    i = find(rows==r);
    axx(i)=subplot('Position',[0.1 0.1+0.3*(3-i) 0.8 0.3]);
    plot(x, Zp(i,:),"r-" ,  ...
         x,Z(i,:), "k-",...
         'LineWidth',2);
    legend(sprintf('%g row of Z',r),sprintf('%g row of Z_p',r))
end

set(axx,'LineWidth',1,'FontSize',14,'YTickLabel','');
set(axx(1:2),'XTickLabel','');
xlabel('Distance(km)');
