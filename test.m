clear, clc

basissetdef = basisread('STO-3G');

atoms = [6 8];

xyz_a0 = [0 0 0;1 0 0];

basis = buildbasis(atoms, xyz_a0, basissetdef);
a = basis(1).a;
b = basis(2).a;
alpha = basis(1).alpha(1);
beta = basis(2).alpha(1);
A = basis(1).A;
B = basis(2).A;
Tprim = int_kinprim(a,b,alpha,beta,A,B);

T = int_kinenergy(basis);
%length(basis)

%basissetdef{atoms(2)};

%basissetdef{atoms(2)}(2).shelltype

