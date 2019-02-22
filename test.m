clear, clc
load('motest.mat')
basissetdef = basisread('STO-3G');

atoms = [1 2];

xyz_a0 = [0 0 0;0 0 1.4632];

basis = buildbasis(atoms, xyz_a0, basissetdef);
settings.basisset = 'STO-3G';
settings.tolEnergy = 1e-8;
settings.tolDensity = 1e-8;
totalcharge = 1;

output = mocalc(atoms,xyz_a0,totalcharge,settings);

F = output.T + output.Vne + output.J - output.K;
err = F*output.C-output.S*output.C*output.epsilon;


%basissetdef{atoms(6)};

%basissetdef{atoms(6)}(6).shelltype

