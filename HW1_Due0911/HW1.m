n = 50;
l = [1 -1 zeros(1,n-2)];
ll = [1 zeros(1,n-1)];
M = toeplitz(ll,l);
M