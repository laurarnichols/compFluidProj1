function [ T, count, a ] = gaussSeidel( n, T, a, b )
%------------------------------------------------------
% Laura Nichols
% 13 February 2018
%
% Uses the Gauss-Seidel algorithm to solve a set of 
% linear equations
%------------------------------------------------------

count = 1;
R = 1;
epsilon = 0.001;

while R > epsilon
    Told = T;
    
    T = zeros(n+2,1);
    for i = 1:n+2
        for j = 1:i-1
            T(i) = T(i) + (-a(i,j)/a(i,i))*T(j);
        end
        
        for j = i+1:n+2
            T(i) = T(i) + (-a(i,j)/a(i,i))*Told(j);
        end
        
        T(i) = T(i) + b(i)/a(i,i);
    end
    
    if count == 1
        R0 = getResidual(Told, T);
    else
        R = getResidual(Told, T)/R0;
    end
    
    count = count + 1;
end
end

