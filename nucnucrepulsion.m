function Vnn = nucnucrepulsion(atoms,xyz_a0)
%Function nucnumrepulsion
%   It calculates the totalnuclear repulsion repulsion energy in the
%   molecule based on equation (17).
%   atoms: list of the elementnumbers
%   xyz_a0: Cartesian coordinates of nuclei, in bohr

Vnn=0
for i=1:numel(atoms)
    for j=i+1:numel(atoms)
        Vnn=Vnn+(atoms(i)*atoms(j))/(norm(xyz_a0(i,:)-xyz_a0(j,:)));
    end
end
end

