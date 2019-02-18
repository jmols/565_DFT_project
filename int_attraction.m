function Vne = int_attraction(atoms,xyz_a0,basis)
M=numel(basis);
Vne=zeros(M,M);
for mu=1:M
    for nu=1:mu
        for c = 1:numel(atoms)
            for k=1:numel(basis(mu).d)
                for l=1:numel(basis(nu).d)
                    Vne=Vne+atoms(c).*basis(mu).d(k)*basis(nu).d(l)*basis(mu).N(k)*basis(nu).N(l)...
                        *int_recursion(xyz_a0,a,b,m,A,B,alpha,beta)
                end
            end
            Vne(mu,nu)=Vne(nu,mu);
        end
    end
end
end


%This is supposed to be equation (40) and (41)
function ab = int_recursion(xyz_a0,a,b,m,A,B,alpha,beta)
ab=0
w=[1 0 0;0 1 0;0 0 1];
for i=1:3
    if (a(i)>=2) && (b(i)>=1)
        ab = ab+((P(i)-A(i))*int_auxiliary(a-w(i,:),b,alpha,beta,A,B)^m...
            +(xyz_a0(i)*w(i,:)-P(i))*int_auxiliary(a-w(i,:),b,alpha,beta,A,B))^(m+1)...
            +(a(i)*w-1)/(2*p)*(int_auxiliary(a-2*w(i,:),b,alpha,beta,A,B)^m...
            -int_auxiliary(a-2*w(i,:),b,alpha,beta,A,B)^(m+1))...
            +b(i)*w(i,:)/(2*p)*(int_auxiliary(a-ang_mo(i,:),b-w(i,:),alpha,beta,A,B)^m...
            -int_auxiliary(a-w(i,:),b-w(i,:),alpha,beta,A,B)^(m+1))
    else if (a(i)>=1) || (b(i)==0)
            ab =r_3d+((P(i)-A(i))*w(i,:))*int_auxiliary(a-w(i,:),b,alpha,beta,A,B)^m...
                +(xyz_a0(i)*w(i,:)-P(i)*int_auxiliary(a-w(i,:),b,alpha,beta,A,B)^(m+1)... +b(i)*ang_mo/(2*p)*(overlap_whatever(a-ang_mo(i,:),b-ang_mo(i,:),alpha,beta,A,B)^m...
                -int_auxiliary(a-w(i,:),b-ang_mo(i,:),alpha,beta,A,B)^(m+1))
        else a(i)==0
            ab = int_auxiliary(a,b,alpha,beta,A,B)
        end
    end
end
for i=1:3
    if (b(i)>=2) || (a(i)>=1)
        ab=ab+((P(i)-B(i))*w(i,:))*int_auxiliary(a,b-w(i,:),alpha,beta,A,B)^m...
            +(xyz_a0(i)*w(i,:)-P(i)*int_auxiliary(a,b-w(i,:),alpha,beta,A,B))^(m+1)...
            +(b(i)*w(i,:)-1)/(2*p)*(int_auxiliary(a,b-2*w(i,:),alpha,beta,A,B)^m...
            -int_auxiliary(a,b-2*w(i,:),alpha,beta,A,B)^(m+1))...
            +a(i)*w(i,:)/(2*p)*(int_auxiliary(a-w(i,:),b-w(i,:),alpha,beta,A,B)^m...
            -int_auxiliary(a-w(i,:),b-w(i,:),alpha,beta,A,B)^(m+1))
    else if (b(i)>=1) || (a(i)==0)
            ab=ab+((P-B(i))*w(i,:)*int_auxiliary(a,b-w(i,:),alpha,beta,A,B))^m...
                +(xyz_a0(i)*w(i,:))-P(i)*int_auxiliary(a,b-w(i,:),alpha,beta,A,B)^(m+1)... +a(i)*ang_mo/(2*p)*(overlap_whatever(a-ang_mo(i,:),b-ang_mo(i,:),alpha,beta,A,B)^m...
                -int_auxiliary(a-w(i,:),b-w(i,:),alpha,beta,A,B)^(m+1)
        else b(i)==0
            ab=int_auxiliary(a,b,alpha,beta,A,B)
        end
    end
end
end

%equation (42) and boys function
function ss=int_auxiliary(a,b,alpha,beta,A,B)
K_AB=exp(-(alpha*beta/(alpha+beta))*((norm(A-B))^2));
p=alpha+beta;
P=(A.*alpha+B.*beta)/(alpha+beta);
for j = 1:numel(atoms)
    T = p * norm(P-xyz_a0(j))^2;
    w=[1 0 0;0 1 0;0 0 1];
    for i=1:3
        if a(0)<0
            ss=0
        elseif b(i)<0
            ss=0
        else a(i)==0&&b(i)==0
            ss=((2*pi)/p)*K_AB*boysF(m,T)
        end
    end
end
end

function y = boysF(m,T)

% Evaluate in terms of the lower incomplete gamma function.
% Warning: Compared to the common definition of the lower
% incomplete gamma function (Wikipedia, Digital Library of
% Mathematical Functions, Wolfram Functions), Matlab's
% gammainc() includes a normalization/regularization factor,
% and the order of arguments is reversed.
mp = m+1/2;
y = (gammainc(T,mp,'lower').*gamma(mp))./(2*T.^mp);

% Limit for T -> 0
threshold = 1e-13;
if numel(T)>1
    idx = abs(T)<threshold;
    y(idx) = 1/(2*m+1);
else
    if abs(T)<threshold
        y = 1./(2*m+1);
    end
end
end



