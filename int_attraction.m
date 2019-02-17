function Vne = int_attraction(atoms,xyz_a0,basis)

M=numel(basis); % number of basis functions
Vne = zeros(M,M);

for u=1:M
    for v=u:M %to fill half the MxM grid
        for c = 1:numel(atoms) %to iterate over atom positions
            for k=1:numel(basis(u).d) %to iterate over contracted basis values
                for l=1:numel(basis(v).d)
                    %using eq 37 along with eq 38 to do the various sums
                    Vne(u,v)=Vne(u,v)- atoms(c)*basis(u).d(k)*basis(v).d(l)...
                        *basis(u).N(k)*basis(v).N(l)*aux_integral(basis(u).a,...
                        basis(v).a,basis(u).alpha(k),basis(v).alpha(l),...
                        basis(u).A,basis(v).A,xyz_a0(c,:),0);
                end
            end
        Vne(v,u) = Vne(u,v); %to fill the remaining half
        end
    end    
end
end