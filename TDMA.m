function [ T, count ] = TDMA( n, T, A, Cp )
%------------------------------------------------------
% Laura Nichols
% 13 February 2018
%
% Uses the TDMA algorithm to solve a set of linear
% equations
%------------------------------------------------------

count = 1;
R = 1;
epsilon = 0.001;

% Loop until converged
while R > epsilon
    Told = T;
    
    % Actual solution equation
    for i = n+1:-1:2
        T(i) = A(i)*T(i+1) + Cp(i);
    end
    
    % Get residual to test
    if count == 1
        R0 = getResidual(Told, T);
    else
        R = getResidual(Told, T)/R0;
    end
    
    count = count + 1;
end
end




