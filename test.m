clear, clc
load('motest.mat')
basissetdef = basisread('STO-3G');

atoms = [8 1 1];

xyz_a0 = [0,0,0.1272;0,0.7581,-0.5086;0,-0.7581,-0.5086];

basis = buildbasis(atoms, xyz_a0, basissetdef);
settings.basisset = 'STO-3G';
settings.tolEnergy = 1e-8;
settings.tolDensity = 1e-8;
settings.method = 'RKS';
settings.ExchFunctional = 'Slater';
settings.CorrFunctional = 'VWN5';
settings.nRadialPoints = 100;
settings.nAngularPoints = 302;
totalcharge = 0;

out = mocalc(atoms,xyz_a0,totalcharge,settings);

%output = mocalc(atoms,xyz_a0,totalcharge,settings);

%basissetdef{atoms(6)};

%basissetdef{atoms(6)}(6).shelltype

