function [basis] = buildbasis(atoms,xyz_a0,basissetdef)
%Builds a basisset given atoms, their cartesian coordinates, and the method
num = length(atoms);
p = 1;

for i = 1:num %iterate over the list of atoms in atoms
    %iterate over the number of basis functions for each atom 
    for n = 1:length(basissetdef{atoms(i)})
        if basissetdef{atoms(i)}(n).shelltype == 'S'
            a=[0 0 0];
            d_row=1; % factor to help determine the contraction coefficients
        elseif basissetdef{atoms(i)}(n).shelltype=='P'
            a=[1 0 0;0 1 0;0 0 1];
            d_row=0;
        elseif basissetdef{atoms(i)}(n).shelltype=='SP'
            a=[0 0 0;1 0 0;0 1 0;0 0 1];
            d_row=1;
        elseif basissetdef{atoms(i)}(n).shelltype=='D'
            a=[2 0 0;1 1 0;1 0 1;0 2 0;0 1 1;0 0 2];
            d_row=-1;
        end
        for q = 1:length(a(:,1))
            basis(p).atom = atoms(i); %puts the atom number in .atoms
            basis(p).A = xyz_a0(i, :); %puts the nuclear coordinates in .xyz_a0
            basis(p).a = a(q,:); %adds the cartesian exponents for the specific shell
            basis(p).alpha = basissetdef{atoms(i)}(n).exponents; %pulls the alphas from basissetdef
            basis(p).d = basissetdef{atoms(i)}(n).coeffs(sum(basis(p).a)+d_row,:); 
       %pulls the d_ks from basissetdef by using that d_row from earlier
            basis(p).N = Norm(basis(p).a(1),basis(p).a(2),basis(p).a(3),basis(p).alpha);
            %the above calulates the N_ks from my Norm function and this line iterates the count
            p = p+1;
        end
    end
end
end

