function val = eval_bf(basisfun,xyz_a0)
%Returns an m element array of the basis function evaluated at m positions
%in space given a basisfunctiona and an mx3 array of positions

points = numel(xyz_a0(:,1)); %determine the amount of positions
contractions = numel(basisfun.d); %determine the amount of primitives in
                                  %the basis function

val = zeros(points,1); %preallocate
for m = 1:points
    for k = 1:contractions
        %Use equation 18 from the math write up
        val(m) = val(m) + basisfun.d(k)*basisfun.N(k)*(xyz_a0(m,1)-...
            basisfun.A(1))^basisfun.a(1)*(xyz_a0(m,2)-basisfun.A(2))...
            ^basisfun.a(2)*(xyz_a0(m,3)-basisfun.A(3))^basisfun.a(3)...
            *exp(-basisfun.alpha(k)*norm(xyz_a0(m,1:3)-basisfun.A)^2);
    end
end
end

