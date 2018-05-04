function [ R ] = getResidual(phi1, phi2)
%------------------------------------------------------
% Laura Nichols
% 14 February 2018
%
% Calculates the residual as the square root of the 
% sum of squares of the difference between 2 solutions
%------------------------------------------------------

R = sqrt(sum(sum((phi2 - phi1).^2)));

end

