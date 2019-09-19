% eda03_06
% Normal curves with the same sigma but different means

% the all these column-vectors have 40 rows
L=40;
sigma=5;
d=[1:L]';

% First normal curve
dbar=10;
norm = 1/(sqrt(2*pi)*sigma);
a=norm*exp( -(d-dbar).^2 / (2*sigma^2) ) / (sigma*sqrt(2*pi));

% Second normal curve
dbar=15;
norm = 1/(sqrt(2*pi)*sigma);
b=norm*exp( -(d-dbar).^2 / (2*sigma^2) ) / (sigma*sqrt(2*pi));

% Third normal curve
dbar=20;
norm = 1/(sqrt(2*pi)*sigma);
c=norm*exp( -(d-dbar).^2 / (2*sigma^2) ) / (sigma*sqrt(2*pi));

% Fourth normal curve
dbar=25;
norm = 1/(sqrt(2*pi)*sigma);
dd=norm*exp( -(d-dbar).^2 / (2*sigma^2) ) / (sigma*sqrt(2*pi));

% Fifth normal curve
dbar=30;
norm = 1/(sqrt(2*pi)*sigma);
e=norm*exp( -(d-dbar).^2 / (2*sigma^2) ) / (sigma*sqrt(2*pi));

% Set up graphics
figure(2);
clf;
set(gca,'LineWidth',2);
hold on;
axis([-L/8, 15*L/8, -L/8, 15*L/8]);
axis equal;
set(gca,'XTick',[]); % turn off horizontal axis
set(gca,'YTick',[]); % turn off vertical axis

% black and white color map
bw=0.9*(256-linspace(0,255,256)')/256;
colormap([bw,bw,bw]);
hold on;

% plot the vecror a
range=max(a)-min(a);
if( range == 0.0 )
    range=1;
end
imagesc( [1, L/16], [L-1, 0], (a-min(a))/range);
text(L/16,-L/16,'a');
hold on;

% plot the vecror b
range=max(b)-min(b);
if( range == 0.0 )
    range=1;
end
imagesc( [1+L/4, L/16+L/4], [L-1, 0], (b-min(b))/range);
text(L/32+L/4,-L/16,'b');
hold on;

% plot the vecror c
range=max(c)-min(c);
if( range == 0.0 )
    range=1;
end
imagesc( [1+2*L/4, L/16+2*L/4], [L-1, 0], (c-min(c))/range);
text(L/32+2*L/4,-L/16,'c');
hold on;

% plot the vecror d
range=max(dd)-min(dd);
if( range == 0.0 )
    range=1;
end
imagesc( [1+3*L/4, L/16+3*L/4], [L-1, 0], (dd-min(dd))/range);
text(L/32+3*L/4,-L/16,'d');
hold on;

% plot the vecror e
range=max(e)-min(e);
if( range == 0.0 )
    range=1;
end
imagesc( [1+4*L/4, L/16+4*L/4], [L-1, 0], (e-min(e))/range);
text(L/32+4*L/4,-L/16,'e');
hold on;

% use simple drawing function that encapsulates all the graphics
eda_draw(' ',a,'caption a',b,'caption b',c,'caption c',dd,'caption d',e,'caption e');