function [ Ac, Cp, a, b ] = setMatrix( T, row, algo, n, BCtype, BCs, gamma, A, dx, dy )
%--------------------------------------------------------
% This formula generates the coefficients in the 
% equation
%-beta(i)*phi(i-1)+D(i)*phi(i)-alpha(i)*(phi(i+1))=C(i)
% for the purpose of solving
%                  d/dx(-gamma*dT/dx) = 0
%
% Functions called:
%           getCoefs - gets coefficients for current,
%                      eastern, and western points
%
% Called by functions:
%           
%
% Input data:
%           algo - choice of algorithm
%                   1 = TDMA
%                   2 = Gauss Seidel
%           n - length of single dimension passed
%           TA and TB - temperature of endpoints
%           gamma - some scalar related to problem
%           A - area
%           dx - point spacing
%
% Output data:
%           Ac and Cp - coefficients for use in TDMA
%           a and b - matrix/vector for use in Gauss
%                     Seidel
%
% Laura Nichols
%--------------------------------------------------------

if isvector(BCs)
    dim = 1;
else
    dim = length(size(BCs));
end

% Initialize return variables, so don't return empty
Ac = 0;
Cp = 0;
a = 0;
b = 0;

% Initialize for speed
C = zeros(n(1)+2,1);
beta = zeros(n(1)+2,1);
D = ones(n(1)+2,1);
alpha = zeros(n(1)+2,1);

% Populate coefficients
for i = 1:n(1)+2
    if dim == 1
        [aP, aE, aW, const] = getCoefs(i, BCtype, BCs, gamma, A, dx);

        C(i) = const;
        beta(i) = aW;
        D(i) = aP;
        alpha(i) = aE;
    elseif dim == 2
        [aP1, aE, aW, const1] = getCoefs(i, BCtype(1,:), BCs(2,:), gamma(row,:), A(row,:), dx(row,:));
        [aP2, aN, aS, const2] = getCoefs(i, BCtype(2,:), BCs(2,:), gamma(:,i), A(:,i), dy(:,i));
        
        C(i) = const1 + const2 + aN*T(row+1,i) + aS*T(row-1,i);
        beta(i) = aW;
        D(i) = aP1 + aP2;
        alpha(i) = aE;
    end
end

if algo == 1

    % Use back solving and forward substitution
    % to make TDMA more accurate

    % Initialize for speed
    Ac = zeros(n(1)+2,1);
    Cp = zeros(n(1)+2,1);

    % Set first values
    Ac(1) = 0;
    Cp(1) = T(row,1);

    % Populate cofficients
    for i = 2:n(1)+1
        Ac(i) = alpha(i)/(D(i) - beta(i)*Ac(i-1));
        Cp(i) = (beta(i)*Cp(i-1) + C(i))/(D(i) - beta(i)*Ac(i-1));
    end


    % Set last values
    Ac(n(1)+2) = 0;
    Cp(n(1)+2) = T(row,end);
elseif algo == 2
    % Put original coefficients in matrix for 
    % Gauss-Seidel
    
    b = C;
    
    % Initialize for speed
    a = zeros(n(1)+2);

    % Set first and last values
    a(1,1) = 1;
    a(n(1)+2,n(1)+2) = 1;

    % Set loop variable
    startCol = 1;
    
    % Populate matrix
    for row = 2:n(1)+1
        a(row,startCol) = -beta(row);
        a(row,startCol+1) = D(row);
        a(row,startCol+2) = -alpha(row);

        startCol = startCol + 1;
    end
end

end

