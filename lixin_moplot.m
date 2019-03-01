function Plot = moplot(atoms,xyz_a0,out,iMO,level) 
% Function that plots a specific MO for a given molecule at a given
% isolevel
[x,y,z] = meshgrid(-5:0.05:5); % generate 3D gird
% S = atoms.*50; %sets size of atoms for visualization's sake
hold on
for N = 1:numel(atoms) %plots the nuclei at their give coordinates
    scatter3(xyz_a0(N,1),xyz_a0(N,2),xyz_a0(N,3),atoms(N)*50,'k','filled');
end

% Build the MO given the MO coeff matrix and the basis functions
% The MO wavefunction Psi
for u=1:numel(out.C(:,iMO))
    dNe=0;
    for k=1:numel(out.basis(u).d)
        dNe=dNe+out.basis(u).d(k)*out.basis(u).N(k)* ...
        exp(-out.basis(u).alpha(k)*(norm([x,y,z]-out.basis(u).A)^2));
    end
    kai=((x-out.basis(u).A(0))^out.basis(u).a(0))*((y-out.basis(u).A(1)) ...
    ^out.basis(u).a(1))*((z-out.basis(u).A(2))^out.basis(u).a(2))*dNe;
    Psi=Psi+out.C(u,iMO)*kai;
end

% Plot the positive and negative isosurfaces
isosurface(x,y,z,Psi,+level)
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

