function out = mocalc(atoms,xyz_a0,totalcharge,settings)
% A function that calls all appropiate functions for an SCF procedure
% Input:
%atoms list of element numbers (array with K elements); e.g. [6 8] for CO 
%xyz_a0 K�3 array of Cartesian coordinates of nuclei, in bohr 
%totalcharge total charge of the molecule, in units of elementary charge 
%settings a structure that contains several fields    
%   .basisset string specifying the basis set, e.g. '6-31G', 'cc-pVDZ', 'STO-3G'   
%   .tolEnergy SCF convergence tolerance for the energy (hartrees)   
%   .tolDensity SCF convergence tolerance for the density (??0 ?3)

%Output:
%out a structure that contains several fields: 
%   .basis list of basis functions, as generated by buildbasis 
%   .S overlap matrix (M�M) .T kinetic energy matrix (M�M) 
%   .Vne electron-nuclear attraction matrix (M�M) 
%   .J matrix of Coulomb integrals (M�M) 
%   .K matrix of exchange integrals (M�M)
%   .ERI 4D array of electron-electron repulsion integrals (M�M�M�M) 
%   .epsilon MO energies (1�M), in hartrees, in ascending order, occupied and virtual orbitals 
%   .C MO coefficient matrix (M�M), of occupied and virtual orbitals,  
%       sorted in ascending order of orbital energy 
%   .P density matrix (M�M) 
%   .E0 electronic ground-state energy of the molecule, in hartrees 
%   .Etot total ground-state energy (including nuclear-nuclear repulsion;  
%       but without the vibrational zero-point energy), in hartrees 


basis = buildbasis(atoms, xyz_a0, basisread(settings.basisset));
S = int_overlap(basis);
T = int_kinenergy(basis);
Vne = int_attraction(atoms, xyz_a0, basis);
ERI = int_repulsion(basis);
N = sum(atoms) - totalcharge;
M = numel(basis);


%initial density matrix construction
%initial guess for F without two electron terms
F_initial = T + Vne;
%solve the Roothan Hall Equation
[C,epsilon] = eig(F_initial,S);
%normalize the MO coeffs
norms = zeros(M);
for k = 1:M
    norms(k) = sqrt(C(:,k)'*S*C(:,k));
end
%norms = C_sort'*S*C_sort;
C_norm = C./norms;
%sort the MO coeffs and epsilon
[epsi_sort, I] = sort(diag(epsilon));
epsi_sort = diag(epsi_sort);
C_sort = C_norm(:, I);
%Build the density matrix from the occupied MO orbitals
C_occupied = C_norm(:,1:N/2);
P = 2*(C_occupied*C_occupied');

%Start of the SCF procedure
maxcycles = 100; %one end condition
E = zeros(maxcycles); %Preallocate
converged = 0;
for n = 1:maxcycles
    %calculate the P for use in the two electron terms and for use 
    %in one of the end conditions
    P_past = 2*(C_occupied*C_occupied');
    
    %calculation of the two electron terms
    J = zeros(M,M);
    K = zeros(M,M);
    for mu = 1:M
        for nu = 1:M
            for kap = 1:M
                for lam = 1:M
                    J(mu,nu) = J(mu,nu) + P_past(mu,nu)*ERI(mu,nu,kap,lam);
                    K(mu,nu) = K(mu,nu) + (1/2).*P_past(mu,nu)*ERI(mu,kap,nu,lam);
                end
            end    
        end
    end
    Vee = J-K;
    
    %building the Fock matrix
    F = F_initial + Vee;
    %solve the Roothan Hall Equation
    [Cprime,epsilonprime] = eig(F,S);
    %normalize the MO coeffs
    for k = 1:M
        norms(k) = sqrt(Cprime(:,k)'*S*Cprime(:,k));
    end
    %norms = C_sort'*S*C_sort;
    C_norm = (Cprime)./norms;
    %sort the MO coeffs and epsilon
    [epsi_sort, I] = sort(diag(epsilonprime));
    C_sort = C_norm(:, I);
    epsilon_occupied = epsi_sort(:,1:N/2);
    C_occupied = C_norm(:,1:N/2);
    P = 2*(C_occupied*C_occupied');
    %Caluclate the F' matrix for the energy calculation
    A = F - (1/2).*Vee;
    E(n) = trace(P*A);
    E0 = E(n);
    %Checks to see if the energy is converged
    if (n>10) && abs(E(n)-E(n-1))<settings.tolEnergy
        for k = 1:M
            for l = 1:M
                %Checks to see if the density is converged
                if abs(P(k,l)-P_past(k,l))<settings.tolDensity
                    converged = 1;
                    
                else
                    continue
                end
            end
        end
    else
        continue
    end
    if converged == 1
        break
    end
end
      

out.basis = basis;
out.S = S;
out.T = T;
out.Vne = Vne;
out.Vee = Vee;
out.J = J;
out.K = K;
out.ERI = ERI;
out.C = C_norm;
out.P = P;
out.epsilon = epsilon_occupied;
out.E0 = E0;

out.Etot = out.E0 + nucnucrepulsion(atoms,xyz_a0);
end

