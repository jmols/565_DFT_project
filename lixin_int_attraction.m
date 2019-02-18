function Vne = int_attraction(atoms,xyz_a0,basis)
%Function int_attraction
%   It calculates the MxM matrix of electron-nuclear attraction energy
%   integrals.

M=numel(basis); % number of basis functions
Vne=zeros(M,M); % preallocate

for u=1:M
    for v=u:M
        for k=1:numel(basis(u).d)
            for l=1:numel(basis(v).d)
                for C=1:numel(atoms)
                    Vne(u,v)=Vne(u,v)-atoms(C)*basis(u).d(k)*basis(v).d(l) ...
                        *basis(u).N(k)*basis(v).N(l)*Vne_primitive(basis(u).a, ...
                        basis(v).a,basis(u).alpha(k),basis(v).alpha(l),basis(u).A, ...
                        basis(v).A,xyz_a0(C,:));
                end
            end
        end
        Vne(v,u)=Vne(u,v);
    end
end
end

function Vp=Vne_primitive(a,b,alpha,beta,A,B,C)
%Function Vne_primitive
%   It calculates the integral for electron-nuclear attranction energy
%   element Vne(u,v) in the function "int_attraction".
%   a,b: vectors of Cartesian coordinates of the nucleus, in bohr.
%   alpha,beta: arrays of radial exponents of primitives, in inverse bohr.
%   A,B,C: vectors of Cartesian coordinates of the nucleus, in bohr.

K_AB=exp(-(alpha*beta/(alpha+beta))*((norm(A-B))^2)); % K_AB from equation (26).
p=alpha+beta; % p from equation (27).
P=(A.*alpha+B.*beta)/(alpha+beta); % P: the weighted midpoint between the two centers A and B, from equation (28).

% Precompute and store the s-type auxiliary integrals in BOYS matrix.
BOYS=0; % boys function will be different for each P{A(u),B(v),alpha(k),beta(l)},C and m values.
m_max=sum(a)+sum(b);
T=p*((norm(P-C))^2);
disp(T);
for i=1:m_max+1
    BOYS(i)=boysF(i-1,T);
end

Vp=auxiliary_int(0,BOYS,P,p,K_AB,a,b,A,B,C);

end

function a_rc_b=auxiliary_int(m,BOYS,P,p,K_AB,a,b,A,B,C)
%Function auxiliary_int
%    Based on the recursive relationships, it calculates the auxiliary
%    integrals for electron-nuclear attraction energy elements in function
%    "Vne_primitive".

check=[0 0];
for i=1:3
    if a(i)==0 && b(i)==0
        continue
    elseif a(i)<0 || b(i)<0
        a_rc_b=0;
        check=-1;
        break
    elseif a(i)>0
        check(1)=1;
    elseif b(i)>0
        check(2)=1;
    end
end

ang_mom=[1,0,0;0,1,0;0,0,1]; % The matrix for the change of angular momentum.

if isequal(check,[0,0]) % only involve s-type functions, i.e. a=[0 0 0] and b=[0 0 0].
    a_rc_b=(2*pi/p)*K_AB*BOYS(m+1);
% Start from reducing a
elseif check(1)==1
    for i=1:3
        if a(i)>0
            a_rc_b=(P(i)-A(i))*auxiliary_int(m,BOYS,P,p,K_AB,a-ang_mom(i,:),b,A,B,C) ...
                +(C(i)-P(i))*auxiliary_int(m+1,BOYS,P,p,K_AB,a-ang_mom(i,:),b,A,B,C) ...
                +((a(i)-1)/(2*p))*(auxiliary_int(m,BOYS,P,p,K_AB,a-2.*ang_mom(i,:),b,A,B,C) ...
                -auxiliary_int(m+1,BOYS,P,p,K_AB,a-2.*ang_mom(i,:),b,A,B,C)) ...
                +(b(i)/(2*p))*(auxiliary_int(m,BOYS,P,p,K_AB,a-ang_mom(i,:),b-ang_mom(i,:),A,B,C) ...
                -auxiliary_int(m+1,BOYS,P,p,K_AB,a-ang_mom(i,:),b-ang_mom(i,:),A,B,C));
            break
        end
    end
elseif check(1)==0 && check(2)==1 % a=[0 0 0], so reduce b
    for i=1:3
        if b(i)>0
            a_rc_b=(P(i)-B(i))*auxiliary_int(m,BOYS,P,p,K_AB,a,b-ang_mom(i,:),A,B,C) ...
                +(C(i)-P(i))*auxiliary_int(m+1,BOYS,P,p,K_AB,a,b-ang_mom(i,:),A,B,C) ...
                +((b(i)-1)/(2*p))*(auxiliary_int(m,BOYS,P,p,K_AB,a,b-2.*ang_mom(i,:),A,B,C) ...
                -auxiliary_int(m+1,BOYS,P,p,K_AB,a,b-2.*ang_mom(i,:),A,B,C));
            break
        end
    end
end

end

