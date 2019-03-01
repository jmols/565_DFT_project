function Plot = moplot(atoms,xyz_a0,out,iMO,level) 
% Function that plots a specific MO for a given molecule at a given
% isolevel
[x,y,z] = meshgrid(-5:0.1:5);
S = atoms.*50; %sets size of atoms for visualization's sake
hold on
for N = 1:numel(atoms) %plots the nuclei at their give coordinates
    scatter3(xyz_a0(N,1),xyz_a0(N,2),xyz_a0(N,3),S(N),'k','filled');
end

C =out.C(:,iMO);
M = numel(C);

% Build the MO given the MO coeff matrix and the basis functions
Psi = zeros(size(x));
for i = 1:M
   Psi = Psi + buildPsi(C(i),x,y,z,out.basis(i));
end

% Plot the positive and negative isosurfaces
isosurface(x,y,z,Psi,level)
isosurface(x,y,z,Psi,-level)

xlabel('r (a0)')
ylabel('r (a0)')
zlabel('r (a0)')
title(['MO ',num2str(iMO), ' with Energy = ',num2str(out.Etot),...
    ' Hatrees at isolevel = ',num2str(level),' (a_0^-^3^/^2)'])
hold off

%Visualization settings
view(3);
camlight
lighting gouraud
colorbar
end

