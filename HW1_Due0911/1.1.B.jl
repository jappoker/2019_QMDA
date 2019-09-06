function f(x)
    p0 = 1.6
    c = 4
    return p0 * exp(-c * x)
end

x = 3.5
y = f(x)
print(y)