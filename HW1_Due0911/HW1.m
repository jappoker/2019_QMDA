D=load('neuse.txt');
t=D(:,1);
d=D(:,2);
allPositive = true;
for i = 1:length(d)
    
    if d(i) < 0
        allPositive = false;
    end

end
allPositive
    