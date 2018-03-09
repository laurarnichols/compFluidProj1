function [ aP, aE, aW, const ] = getCoefs( i, BCtype, BCs, gamma, A, dx)
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

if ~isvector(BCs)
    error(['Error: getCoefs expects to only deal' ...
          ' with one dimension.']);
end

const = 0;

% Take average of neighbors to get gamma and area at
% eastern and western faces of volume
if i == 1 %#ok<*IJCL>
    if BCtype(1) == 1
        aW = 0;
        aP = 1;
        aE = 0;
    elseif BCtype(1) == 2
        aW = 0;
        aP = 0;
        aE = 0;
    else 
        error('That BC type isn''t allowed.');
    end
elseif i == length(A)
    if BCtype(2) == 1
        aW = 0;
        aP = 1;
        aE = 0;
    elseif BCtype(2) == 2
        aW = 0;
        aP = 0;
        aE = 0;
    else
        error('That BC type isn''t allowed.');
    end
else 
    if i == 2
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
    elseif i == length(A) - 1
        gamma_w = (gamma(i-1) + gamma(i)) / 2;
        A_w = (A(i-1) + A(i)) / 2;
        
        if BCtype(2) == 1
            gamma_e = gamma(end);
            A_e = A(end);
        elseif BCtype(2) == 2
            gamma_e = gamma(end);
            A_e = A(end);
            const = -BCs(2);
        else
            error('That BC type isn''t allowed.');
        end
    else
        gamma_e = (gamma(i+1) + gamma(i)) / 2;
        A_e = (A(i+1) + A(i)) / 2;
        gamma_w = (gamma(i-1) + gamma(i)) / 2;
        A_w = (A(i-1) + A(i)) / 2;
    end
    
    % Define coefficients
    aP = gamma_e*A_e/dx(i) + gamma_w*A_w/dx(i-1);
    aE = gamma_e*A_e/dx(i);
    aW = gamma_w*A_w/dx(i-1);
end

end

