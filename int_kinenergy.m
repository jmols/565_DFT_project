function T = int_kinenergy(basis)
%Function int_kinenergy
%   It calculates the MxM matrix of kinetic-energy integrals, based on the
%   equation (34)
M=numel(basis); % number of basis functions
for u=1:M
    for v=u:M
        t=0;
        for k=1:numel(basis(u).d)
            for l=1:numel(basis(v).d)
                t=t+basis(u).d(k)*basis(v).d(l)*basis(u).N(k)*basis(v).N(l)*I_kinenergy(basis(u).a,basis(v).a,basis(u).alpha(k),basis(v).alpha(l),basis(u).A,basis(v).A);
            end
        end
        T(u,v)=t;
        T(v,u)=t;
    end
end
end

function t_3d=I_kinenergy(a,b,alpha,beta,A,B)
%Function I_kinenery
%   It calculates the 3D integrals of the kinetic-energy integrals, based
%   on the equation (35) & (36). (for each k,l,a,b)

t_3d=0;
ang_mo=[1 0 0;0 1 0;0 0 1]; % used for 1 or 2 units of ?angular momentum? added/subtracted to/from the w-th component of the function |a] respectively.
for i=1:3
    if b(i)>=2
        t_3d=t_3d+beta*(2*b(i) +1)*overlap_primitive(a,b,alpha,beta,A,B)-2*beta^2*overlap_primitive(a,b+2*ang_mo(i,:),alpha,beta,A,B)-(1/2)*b(i)*(b(i)-1)*overlap_primitive(a,b-2*ang_mo(i,:),alpha,beta,A,B);
    else
        t_3d=t_3d+beta*(2*b(i) +1)*overlap_primitive(a,b,alpha,beta,A,B)-2*beta^2*overlap_primitive(a,b+2*ang_mo(i,:),alpha,beta,A,B);
    end
end
end
