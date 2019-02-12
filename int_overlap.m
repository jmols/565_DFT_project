function S = int_overlap(basis)
%Function int_overlap
%   It takes the information of the basis functions from basis (M
%   elements), and calculates the overlap integrals (MxM matrix).

M=numel(basis); % number of basis functions
for u=1:M
    for v=u:M
        s=0;
        for k=1:numel(basis(u).d)
            for l=1:numel(basis(v).d)
                s=s+basis(u).d(k)*basis(v).d(l)*basis(u).N(k)*basis(v).N(l)*overlap_primitive(basis(u).a,basis(v).a,basis(u).alpha(k),basis(v).alpha(l),basis(u).A,basis(v).A);
            end
        end
        S(u,v)=s;
        S(v,u)=s;
    end
end
end

