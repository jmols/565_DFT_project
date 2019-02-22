function out=mocalc(atoms,xyz_a0,totalcharge,settings)
%Function mocalc
%   It contains the SCF algorithm. It reads in the basis set definition and
%   calculate all the integrals by calling the corresponding functions, and
%   then returns to the SCF results.

% Read the basis function
basissetdef=basisread(settings.basisset);
out.basis=buildbasis(atoms,xyz_a0,basissetdef);

N=sum(atoms)-totalcharge; % Calculate the number of electrons
M=numel(out.basis); % Calculate the number of basis functions

out.S=int_overlap(out.basis); % overlap matrix
out.T=int_kinenergy(out.basis); % kinetic energy matrix
out.Vne=int_attraction(atoms,xyz_a0,out.basis); % electron-nuclear attraction matrix
out.ERI = int_repulsion(out.basis); % two-electron integrals
Vnn=nucnucrepulsion(atoms,xyz_a0); % nuclear-nuclear repulsion

% Estimates a starting density matrix by neglecting all electron-electron
% repulsion.
guess=SCF_oneloop(N,M,out.T,out.Vne,out.S,out.ERI,Vnn,zeros(M),0);
guess.e_Etot=1;
guess.e_P=1;
while guess.e_Etot>=settings.tolEnergy && guess.e_P>=settings.tolDensity
    guess=SCF_oneloop(N,M,out.T,out.Vne,out.S,out.ERI,Vnn,guess.P,guess.Etot);
end
% fprintf('SCF converged!');
out.J=guess.J;
out.K=guess.K;
out.epsilon=guess.epsilon;
out.C=guess.C;
out.P=guess.P;
out.E0=guess.E0;
out.Etot=guess.Etot;
end

function norm=normalization(C,S)
%Function normalization
%   It takes an unnormalized matrix C and returns to a normalized matrix.

for k = 1:numel(C(1,:))
  norms(k) = sqrt(C(:,k)'*S*C(:,k));
end
norm = C./norms;
end

function P=density_matrix(C,N)
%Function density_matrix
%   It calculates the density matrix from the MO coefficient matrix C,
%   based on the RHF. n is the number of atoms.
%   The density matrix can be obtained by picking the N/2 lowest-energy MOs
%   (occupied MOs, MxN/2 matrix) and calculating 2*C_occupied*C_occupied',
%   resulting in a matrix of size MxM. N is the number of electrons, i.e.
%   the number of MOs.

C_occupied=C(:,1:N/2); % pick occupied MOs
P=2*(C_occupied*C_occupied');
% Assure that the trace of the density matrix equals the number of electrons
%trace(P)
end

function [C_s,e_s]=sort_C_e(C,e)
%Function sort
%   It first sorts the energies in increasing order and then sort the MO
%   matrix accordingly.
%   C, e are MxM matrices.

[e_s,index]=sort(diag(e));
C_s=C(:,index);
end

function SCF=SCF_oneloop(N,M,T,Vne,S,ERI,Vnn,P_0,Etot_0)
% Function SCF_oneloop
%   It takes the necessary parameters for RHF calculation as well as the
%   density matrix and total energy value from the previous RHF
%   calculation.

% a) Build the Fock matrix
F=T+Vne; % one-electron part
SCF.J=zeros(M); % preallocate
SCF.K=zeros(M);
for q=1:M
    for p=1:M
        for s=1:M
            for r=1:M
                SCF.J(q,p)=SCF.J(q,p)+ERI(q,p,s,r)*P_0(r,s);
                SCF.K(q,p)=SCF.K(q,p)+(1/2)*ERI(q,r,s,p)*P_0(r,s);
            end
        end
        F(q,p)=F(q,p)+SCF.J(q,p)-SCF.K(q,p);        
    end
end

% b) Solve for the MO coefficients C
[C,e]=eig(F,S);
[C_s,epsilon]=sort_C_e(C,e); % Sort C and e
SCF.C=normalization(C_s,S); % normalize C_0
SCF.epsilon=epsilon';

% c) Build the new density matrix P
SCF.P=density_matrix(SCF.C,N); % RHF density matrix

% d) Calculate errors between two results for convergence judgement
SCF.E0=0; % electronic ground-state energy of the molecule
for u=1:M
    for v=1:M
        SCF.E0=SCF.E0+SCF.P(u,v)*(T(u,v)+Vne(u,v)+(1/2)*(SCF.J(u,v)-SCF.K(u,v)));
    end
end
SCF.Etot=SCF.E0+Vnn; % total ground-state energy

SCF.e_Etot=abs(SCF.Etot-Etot_0); % error between new Etot and the previous Etot
error_P=SCF.P-P_0;
SCF.e_P=max(abs(error_P(:))); % error between new P and the previous P

end



