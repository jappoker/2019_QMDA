% hanntaper
%
% Return a Hann(ing) cosine bell taper centered at xc
% and of half-width deltax as a function of x

function tap=hanntaper(x,xc,deltax)

xstart=xc-deltax;
xend=xc+deltax;
istart=find(x>=xstart,1,'first');
iend=find(x<=xend,1,'last');
tap=zeros(length(x),1);
for i=istart:iend
    tap(i)=0.5*(1-cos(pi*(x(i)-xstart)/deltax));
end