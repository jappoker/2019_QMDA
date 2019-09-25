% samplebootstrap

function [xbs,ybs]=samplebootstrap(x,y)

n=length(x);
index=randi(n,n,1);
xbs=x(index);
ybs=y(index);
