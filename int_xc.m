function [Vxc,Exc,rhoInt]=int_xc(basis,P,grid,ExchFunctional,CorrFunctional)
%Function int_xc
%   It returns the matrix of exchange-correlation integrals (MxM) in Vxc,
%   the exchange-correlation energy in Exc, and the integral over the
%   electron density in rhoInt.

m=numel(grid.xyz(:,1)); % number of molecular grid points
M=numel(basis); % number of basis functions

%1) Use 'eval_bf' to evaluate all basis functions at each grid points r_k.
val_bf=zeros(M,m); % preallocate
for u=1:M
    val_bf(u,:)=eval_bf(basis(u),grid.xyz);
end

%2) Calculate the density rho at each grid point k based on eq.(69).
rho=zeros(1,m); % preallocate
for k=1:m
    for r=1:M
        for s=1:M
            rho(k)=rho(k)+val_bf(r,k)*val_bf(s,k)*P(r,s);
        end
    end
end

%3) the exchange and correlation potential, energy density, density
%integral, i.e. Vxc, Exc, and rhoInt.
if strcmp(ExchFunctional,'Slater')
    Cx = (3/4)*(3/pi)^(1/3);
    ex=-(rho.^(1/3)).*Cx; % Slater(S) exchange energy density
    Vx=(4/3)*ex; % Slater exchange potential
else
    fprintf('Please reset the exchange functional.');
end

if strcmp(CorrFunctional,'VWN3')
    b=13.0720;
    c=42.7198;
    x0=-0.409286;
elseif strcmp(CorrFunctional,'VWN5')
    b=3.72744;
    c=12.9352;
    x0=-0.10498;
else
    fprintf('Please reset the correlation functional.');
end
ec=zeros(1,m); % preallocate
Vc=zeros(1,m); 
for k=1:m
    [ec(k),Vc(k)]=VWN(rho(k),b,c,x0);
end
exc=ex+ec;
V_xc=Vx+Vc;
Vxc=zeros(M,M); % preallocate
for u=1:M
    for v=1:M
        for k=1:m
            Vxc(u,v)=Vxc(u,v)+grid.weights(k)*val_bf(u,k)*V_xc(k)*val_bf(v,k);
        end
    end
end
Exc=0;
rhoInt=0;
for k=1:m
    Exc=Exc+grid.weights(k)*exc(k)*rho(k);
    rhoInt=rhoInt+grid.weights(k)*rho(k);
end
end

function [ec,Vc]=VWN(rho,b,c,x0)
%Function VWN
%   It calculates the Vosko-Wilk-Nusair (VWN) correlation and returns the
%   energy density and potential for each rho value.

A=0.0310907; % in hartrees
Q=sqrt((4*c)-(b^2));
if rho>=1e-16
    x=(3/(4*pi*rho))^(1/6);
    eta=atan(Q/(2*x+b));
    ec=A*(log(x^2/xi(x,b,c))+2*b*eta/Q-(b*x0/xi(x0,b,c))* ...
        (log((x-x0)^2/xi(x,b,c))+2*(2*x0+b)*eta/Q));
    Vc=ec-(A/3)*(c*(x-x0)-b*x*x0)/(xi(x,b,c)*(x-x0));
else
    ec=0;
    Vc=0;
end
end

function v=xi(zeta,b,c)
%Function xi
%   It calculate the function for xi value v from the input zeta.
v=zeta^2+b*zeta+c;
end

