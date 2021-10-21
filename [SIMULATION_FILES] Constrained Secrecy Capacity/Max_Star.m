function z = Max_Star(a,b)

if a>b
    z = a +log(1+exp(b-a));
elseif b> a
    z = b + log(1+exp(a-b));
else
    z = a;
end
