function f=fac2(n)
%Function fac2
%    It takes a number and return its double factorial value.

if (n==1||n==0)
    f=1;
else
    if (mod(n,2)==0) % even number
        f=prod(2:2:n);
    else % odd number
        f=prod(1:2:n);
    end
end
end