function [N] = Norm(ax,ay,az,alpha)
%Calculates the normalization constant Nk
a = ax + ay + az;
denom = sqrt(fac2(2*ax-1)*fac2(2*ay-1)*fac2(2*az-1));
N = (2/pi)^(3/4)*(2^(a).*alpha.^((2*(a)+3)/4))/denom;
end

