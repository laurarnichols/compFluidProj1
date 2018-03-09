function [ T, count ] = solve( algo, n, T0, BCtype, BCs, k, A, dx, dy )

% Get the number of dimensions
if isvector(T0)
    dim = 1;
else
    dim = length(size(T0));
end

%------------------------------------------------------
% Solve for T using different methods

if algo == 1 && dim == 1 
    [Ac, Cp, ~, ~] = setMatrix(T0, 1, algo, n, BCtype, BCs, k, A, dx);
    [T, count] = TDMA(n, T0, Ac, Cp);
elseif algo == 2 && dim == 1
    [~, ~, a, b] = setMatrix(T0, 1, algo, n, BCtype, BCs, k, A, dx);
    [T, count] = gaussSeidel(n, T0, a, b);
elseif dim == 2
    error('Shouldn''t get here');
    count = 1;
    R = 1;
    epsilon = 0.01;
    T = T0;
    
    while R > epsilon
        for j = 2:n(2)+1
            [Ac, Cp, ~, ~] = setMatrix(T, j, algo, n, BCtype, BCs, k, A, dx, dy);
            [T(j,:), count] = TDMA(n, T0(j,:), Ac, Cp);
        end
        
        T = fixBounds(T, BCtype, BCs, k, A, dx, dy);
        
        Told = T;
        
        if count == 1
            R0 = getResidual(Told, T);
        else
            R = getResidual(Told, T)/R0;
        end

        count = count + 1;
    end
else
    error(['Sorry, solve is not able to handle'... 
             ' that option at this time.']);
end

end

