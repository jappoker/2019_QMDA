### Student Info

Jingyi Zhuang

Uni: jz2907

Email: jz2907@columbia.edu

### Problems & Solutions

1.1 (A)

Answer: 

```
46.5000
```

Code:

```matlab
a=2;
b=4;
c=8;
x=3.5;
y=a*x^2+b*x+c;
y
```

---

1.1 (B)

Answer: 

```
1.3304e-06
```

Script:

```matlab
p0 = 1.6;
c = 4;
x = 3.5;
p = p0 * exp(-c*x);
p
```

---

1.1 (C)

Answer:

```
2.0602
```

Script:

```matlab
h = 4;
theta = 31;
z = h * sin(theta/180 * pi);
z
```

---

1.1 (D)

Answer:

```
296.7580
```

Script:

```julia
h = 6.9;
r = 3.7;
v = pi * h * r ^2;
v
```

---

1.2

Answer:

```
C =

    31    31
    28    29
    31    31
    30    30
    31    31
    30    30
    31    31
    31    31
    30    30
    31    31
    30    30
    31    31
```

Script:

```matlab
a = [31;28;31;30;31;30;31;31;30;31;30;31];
b = [31;29;31;30;31;30;31;31;30;31;30;31];
C = [a b];
C
```

---

1.3

Answer:

```
x =

    11
    10
     8
     5
```

Script:

```matlab
y = [1;2;3;5];
M = [1 -1 0 0; 0 1 -1 0; 0 0 1 -1; 0 0 0 1];
x = M\y;
x
```

---

1.4

Script:

```matlab
n = 50;
l = [1 -1 zeros(1,n-2)];
ll = [1 zeros(1,n-1)];
M = toeplitz(ll,l);
M
```

---

1.5

Answer:

```
allPositive =

  logical

   1
```

None of the Neuse River discharge data is negative.

Script:

```matlab
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
```

