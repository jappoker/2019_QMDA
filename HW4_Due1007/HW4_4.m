clearvars;
D = load('linedata01.txt');

x = D(:,1);
d = D(:,2);
N = length(x);

for i = 2 : 4
    G = ones(N,i+1);
    for n = 1 : i
        G(:,n+1) = x.^n;
    end
    
    figure(i-1);
    
    set(gca,'FontSize',16,'LineWidth',1);
    hold on;
    
    plot(x,d,"o",'LineWidth',2,"MarkerSize",10);
    hold on;
    
    mest = (G'*G)\(G'*d);
    
    nline = 100;
    xline=linspace(min(x),max(x),nline);
    Gline= ones(nline,i+1);
    for n = 1 : i
        Gline(:,n+1) = xline.^n;
    end
    
    plot(xline,Gline*mest,"-",'LineWidth',2,"MarkerSize",10);

    fprintf("The coefficients of polynomial of degree %g are\n",i);
    disp(mest);
    
    e = d - G*mest;
    E = e' * e;
    sigma2d= E / (length(d)- (i+1));
    Cm = sigma2d * inv(G'*G);
    
    sigma2m = diag(Cm);
    
    mu = mest + 2 * sqrt(sigma2m);
    ml = mest - 2 * sqrt(sigma2m);
    
    
    
    fprintf("The 0.95 confidence upper bound coefficients of polynomial of degree %g are\n",i);
    disp(mu);
    
    fprintf("The 0.95 confidence lower bound coefficients of polynomial of degree %g are\n",i);
    disp(ml);
    
    
    plot(xline,Gline*mu,"-",'LineWidth',2);
    plot(xline,Gline*ml,"-",'LineWidth',2);

    
    legend({"data", sprintf("poly: degree of %g",i), "upper bound of 0.95 confidence", "lower bound of 0.95 confidence" },'Location','northwest')
    title( sprintf("poly: degree of %g",i));
end
