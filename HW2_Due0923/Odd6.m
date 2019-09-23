clearvars;
s = [0:0.01:1]';
p = 0.03087.*s ./(0.03087.*s + (1-s)*0.00081);
plot(s,p,...
    'LineWidth', 1);

xlabel(sprintf('P(docter thinks the patient has Strep Throat) <--> the prior'));
ylabel('P(has the strep throat)');