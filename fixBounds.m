function [ T ] = fixBounds( T, BCtype, BCs, k, A, dx, dy )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if BCtype(1,1) == 2
    T(:,1) = dx(:,1)*BCs(1,1)./(k(:,1).*A(:,1)) - T(:,2);
elseif BCtype(1,2) == 2
    T(:,end) = -dx(:,end)*BCs(1,2)./(k(:,end).*A(:,end)) + T(:,end-1);
elseif BCtype(2,2) == 2
    T(1,:) = dy(1,:)*BCs(2,2)./(k(1,:).*A(1,:)) - T(2,:);
elseif BCtype(2,1) == 2
    T(end,:) = -dy(end,:)*BCs(2,1)./(k(end,:).*A(end,:)) + T(end-1,:);
end

end

