function [ aP, aE, aW, const ] = getCoefs( i, BCtype, BCs, gamma, A, dx, u, rho)
%--------------------------------------------------------
% This function assumes you are only calculating inner
% points so that you don't go outside the bounds.
% 
% This function also assumes a 1D problem, so extras 
% will have to be done outside function.
%
% Functions called:
%           
%
% Called by functions:
%           setMatrix - puts all coefficients into 
%                       matrix for solving
%
% Input data:
%           i - index
%           gamma - some scalar related to problem
%           A - area
%           dx - point spacing
%
% Output data:
%           aP - coefficient for current point
%           aE - coefficient for eastern point
%           aW - coefficient for western point
%
% Laura Nichols
%--------------------------------------------------------

% Check if input has correct size
if ~isvector(BCs)
    error(['Error: getCoefs expects to only deal' ...
          ' with one dimension.']);
end

const = 0;

% Treat inner and outer points seperately
if i == 1 %#ok<*IJCL>
    % Boundary
    if BCtype(1) == 1
        aW = 0;
        aP = 1;
        aE = 0;
        const = BCs(1);
    elseif BCtype(1) == 2
        aW = 0;
        aP = 0;
        aE = 0;
    else 
        error('That BC type isn''t allowed.');
    end
elseif i == length(A)
    % Boundary
    if BCtype(2) == 1
        aW = 0;
        aP = 1;
        aE = 0;
        const = BCs(2);
    elseif BCtype(2) == 2
        aW = 0;
        aP = 0;
        aE = 0;
    else
        error('That BC type isn''t allowed.');
    end
else 
    % Next to boundary
    if i == 2
        % Take average of neighbors to get gamma and area at
        % eastern face of volume
        gamma_e = (gamma(i+1) + gamma(i)) / 2;
        A_e = (A(i+1) + A(i)) / 2;
        
        if BCtype(1) == 1
            gamma_w = gamma(1);
            A_w = A(1);
        elseif BCtype(1) == 2
            gamma_w = 0;
            A_w = 0;
            const = BCs(1);
        else
            error('That BC type isn''t allowed.');
        end
    % Next to boundary
    elseif i == length(A) - 1
        % Take average of neighbors to get gamma and area at
        % western face of volume
        gamma_w = (gamma(i-1) + gamma(i)) / 2;
        A_w = (A(i-1) + A(i)) / 2;
        
        if BCtype(2) == 1
            gamma_e = gamma(end);
            A_e = A(end);
        elseif BCtype(2) == 2
            gamma_e = 0;
            A_e = 0;
            const = -BCs(2);
        else
            error('That BC type isn''t allowed.');
        end
    else
        % Take average of neighbors to get gamma and area at
        % eastern and western faces of volume
        gamma_e = (gamma(i+1) + gamma(i)) / 2;
        A_e = (A(i+1) + A(i)) / 2;
        gamma_w = (gamma(i-1) + gamma(i)) / 2;
        A_w = (A(i-1) + A(i)) / 2;
    end
    
    D_e = gamma_e*A_e/dx(i);
    D_w = gamma_w*A_w/dx(i-1);
    
    F_e = rho(i)*u(i)*A_e;
    F_w = rho(i-1)*u(i-1)*A_w;
    %input('Waiting...');
%     F_e = 0;
%     F_w = 0;
    
    % Define coefficients
    aP = D_w + D_e + F_e;
    aE = D_e + max(0,-F_e);
    aW = D_w + max(F_w, 0);
end

end

