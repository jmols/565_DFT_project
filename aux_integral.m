function Z = aux_integral(a, b, alpha, beta, A, B, C, m)
Kab=exp(-(alpha*beta/(alpha+beta))*((norm(A-B))^2)); % K_AB from equation (26).
p=alpha+beta; % p from equation (27).
% P: the weighted midpoint between the two centers A and B, from equation (28).
P=(A.*alpha+B.*beta)/(alpha+beta);
T = p * norm(P-C)^2; % T from equation 42

w = [1 0 0; 0 1 0; 0 0 1]; %to subtract from the ang mo

for i=1:3 % to check to see if any elements are below zero
    if a(i)<0
        Z=0;
    elseif b(i)<0
        Z=0;
    end
end

if sum(a)==0 && sum(b)==0 %implementing eq 42
    Z = (2*pi)/p * Kab * boysF(m,T);
    
else
    for i = 1:3 % using equations 40 and 41 if any ang mo elements are above zero
        if a(i)>0
            Z = (P(i)-A(i))*aux_integral(a-w(i,:),b,alpha,beta,A,B,C,m) + ...
                (C(i)-P(i))*aux_integral(a-w(i,:),b,alpha,beta,A,B,C,m+1) + ...
                (a(i)-1)/(2*p)*(aux_integral(a-2*w(i,:),b,alpha,beta,A,B,C,m) - ...
                aux_integral(a-2*w(i,:),b,alpha,beta,A,B,C,m+1)) + ...
                b/(2*p)*(aux_integral(a-w(i,:),b-w(i,:),alpha,beta,A,B,C,m) - ...
                aux_integral(a-w(i,:),b-w(i,:),alpha,beta,A,B,C,m+1));
        elseif b(i)>0
            Z = (P(i)-B(i))*aux_integral(a,b-w(i,:),alpha,beta,A,B,C,m) + ...
                (C(i)-P(i))*aux_integral(a,b-w(i,:),alpha,beta,A,B,C,m+1) + ...
                (b(i)-1)/(2*p)*(aux_integral(a,b-2*w(i,:),alpha,beta,A,B,C,m) - ...
                aux_integral(a,b-2*w(i,:),alpha,beta,A,B,C,m+1));
        end
    end
end
end

