function S = int_overlap(basis)
%Function int_overlap
%   It takes the information of the basis functions from basis (M
%   elements), and calculates the overlap integrals (MxM matrix).

M=numel(basis); % number of basis functions
S = zeros(M,M); % to pre-allocate the MxM grid
for u=1:M
    for v=u:M %to fill the MxM grid
        for k=1:numel(basis(u).d)
            for l=1:numel(basis(v).d) % to account for the contractions
                % add up the overlaps to account for the contracted basis
                S(u,v)=S(u,v)+basis(u).d(k)*basis(v).d(l)*basis(u).N(k)...
                    *basis(v).N(l)*overlap_primitive(basis(u).a,basis(v).a,...
                    basis(u).alpha(k),basis(v).alpha(l),basis(u).A,basis(v).A);
            end
        end
        S(v,u) = S(u,v);
    end
end
end

