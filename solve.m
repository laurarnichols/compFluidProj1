function [ T, count ] = solve( algo, n, T0, BCtype, BCs, k, A, dx, dy, u, rho, epsilon )

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
    % Use TDMA if dimension is 2
    loopCount = 1;
    R = 1;
    %epsilon = 5/n(1);
    T = T0;
    count = 0;
    
    % Loop until entire solution converges
    while R > epsilon || loopCount < 100
        if mod(count,2) == 0
            % Loop over rows
            for j = 2:n(2)+1
                [Ac, Cp, ~, ~] = setMatrix(T, j, algo, n, BCtype, BCs, k, A, dx, dy, u, rho, 1);
                % Solve single row
                [T(j,:), newCount] = TDMA(n(1), T(j,:), Ac, Cp);
                count = count + newCount;
            end
        else
            % Loop over columns
            for j = 2:n(2)+1
                [Ac, Cp, ~, ~] = setMatrix(T, j, algo, n, BCtype, BCs, k, A, dx, dy, u, rho, 2);
                % Solve single row
                [T(:,j), newCount] = TDMA(n(1), T(:,j), Ac, Cp);
                count = count + newCount;
            end
        end
        
        % Set boundary temps
        %T = fixBounds(T, BCtype, BCs, k, A, dx, dy);
                
        if loopCount < 100
            %R0 = getResidual(T0, T);
            R = 1;
        else
            R = getResidual(Told, T);
        end

        Told = T;
        loopCount = loopCount + 1;
    end
else
    error(['Sorry, solve is not able to handle'... 
             ' that option at this time.']);
end

end

