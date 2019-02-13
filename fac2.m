function [output] = fact2(n)
%Double Factorial function
%    ?? = ? ? (? ? 2) ? … ? 1 for ? > 0, and 1 otherwise
if n > 0
    if rem(n,2) == 0
        output = prod([2:2:n]);
    else
        output = prod([1:2:n]);
    end
else
    output = 1;
end
end

