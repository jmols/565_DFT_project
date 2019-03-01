function Func = buildPsi(C,x,y,z,basis)
% Builds an MO function given an MO coeff and gridpoints
contractions = numel(basis.d);
Func = 0;
%calculates the angular bit of the function
ang = C.*(x-basis.A(1)).^basis.a(1).*(y-basis.A(2)).^basis.a(2).*...
    (z-basis.A(3)).^basis.a(3);

%calulates the MO function with the contracted radial part
for k = 1:contractions
    Func = Func + ang.*basis.d(k).*basis.N(k).*exp(-basis.alpha(k).*...
        (x-basis.A(1)).^2).*exp(-basis.alpha(k).*(y-basis.A(2)).^2).*...
        exp(-basis.alpha(k).*(z-basis.A(3)).^2);
end
end

