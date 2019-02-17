clear, clc

basissetdef = basisread('STO-3G');

atoms = [6 8];

xyz_a0 = [0 0 0;0.5 0 0];

basis = buildbasis(atoms, xyz_a0, basissetdef);


Z = atoms(1)*basis(1).d(1)*basis(6).d(3)*basis(1).N(1)*basis(6).N(3)*...
    aux_integral(basis(1).a,basis(6).a, basis(1).alpha(1),basis(6).alpha(3),...
    basis(1).A, basis(6).A,xyz_a0(1,:),0);

Vne = int_attraction(atoms,xyz_a0,basis);
%length(basis)

%basissetdef{atoms(6)};

%basissetdef{atoms(6)}(6).shelltype

