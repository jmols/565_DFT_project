function basis = buildbasis(atoms,xyz_a0,basissetdef)
%Function buildbasis 
%   It takes a list of atoms and their coordinates (xyz_a0), and the basis
%   set information read in by basisread, and generates a full list of
%   basis functions for a given molecule.

count=0;
for i=1:numel(atoms)
    for j=1:numel(basissetdef{atoms(i)})
        if basissetdef{atoms(i)}(j).shelltype=='S'
            a=[0 0 0];
            d_row=1; % help to determine the contraction coefficients
        elseif basissetdef{atoms(i)}(j).shelltype=='P'
            a=[1 0 0;0 1 0;0 0 1];
            d_row=0;
        elseif basissetdef{atoms(i)}(j).shelltype=='SP'
            a=[0 0 0;1 0 0;0 1 0;0 0 1];
            d_row=1;
        elseif basissetdef{atoms(i)}(j).shelltype=='D'
            a=[2 0 0;1 1 0;1 0 1;0 2 0;0 1 1;0 0 2];
            d_row=-1;
        end
        for p=1:numel(a(:,1))
            count=count+1;
            nbf.atom=atoms(i); % element number
            nbf.A=xyz_a0(i,:); % vector of Cartesian coordinates of the nucleus, in bohr
            nbf.a=a(p,:); % vector of Cartesian exponents ([0 0 0] for s, [1 0 0] for px, etc.)
            nbf.alpha=basissetdef{atoms(i)}(j).exponents; % array of radial exponents of primitives, in inverse bohr
            nbf.d=basissetdef{atoms(i)}(j).coeffs(sum(nbf.a)+d_row,:); % array of contraction coefficients
            nbf.N=0;
            for k=1:numel(basissetdef{atoms(i)}(j).exponents)
                nbf.N(k)=((2/pi)^(3/4)*(2^sum(a(p,:)))*(basissetdef{atoms(i)}(j).exponents(k)^((2*sum(a(p,:))+3)/4)))/sqrt(fac2(2*a(p,1)-1)*fac2(2*a(p,2)-1)*fac2(2*a(p,3)-1));
            end% normalization constants
            basis(count)=nbf;
        end
    end 
end
end