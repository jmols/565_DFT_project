function Sp = overlap_primitive(a,b,alpha,beta,A,B)
%Function overlap_primitive
%   It calculates the overlap integral between two unnormalized Gaussian
%   primitives, based on equation (31) in the file
%   MolecularOrbitalsFromScratch.pdf.
%   a,b: Cartesian exponent vectors, [ax ay az] and [bx by bz].
%   alpha,beta: radial exponents for the two primitives, in 1/bohr^2.
%   A,B: center vectors,[Ax Ay Az] and [Bx By Bz], in bohr.

K_AB=exp(-(alpha*beta/(alpha+beta))*((norm(A-B))^2)); % K_AB from equation (26).
p=alpha+beta; % p from equation (27).
P=(A.*alpha+B.*beta)/(alpha+beta); % P: the weighted midpoint between the two centers A and B, from equation (28).

Sp=((pi/p)^(3/2))*K_AB;
for i=1:3
    Sp=Sp*I_overl_prim(a(i),b(i),P(i),A(i),B(i),p);
end
end

function i_1d=I_overl_prim(a,b,P,A,B,p)
%Function I_overl_prim
%   It calculates the 1D integrals for each coordinates, based on the
%   equation (32) in the file MolecularOrbitalsFromScratch.pdf.

i_1d=0;
for i=0:1:(a+b)/2
    i_1d=i_1d+f_overl_prim(2*i,a,b,P,A,B)*fac2(2*i-1)/((2*p)^i);
end
end

function f=f_overl_prim(k,a,b,P,A,B)
%Function f_overl_prim
%   It subsitutes one part of the calculation for the 1D integral for the
%   function overlap_primitive, based on the equation (33) in the file
%   MolecularOrbitalsFromScratch.pdf.

f=0;
for j=max(0,k-a):min(k,b)
    f=f+nchoosek(a,k-j)*nchoosek(b,j)*((P-A)^(a-k+j))*((P-B)^(b-j));
end
end