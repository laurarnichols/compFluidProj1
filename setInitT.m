function [ T0 ] = setInitT( n, BCtype, BCs )
% Set initial temperature guess based on given 
% boundary conditions

% T0 = (BCs(1,1) + BCs(2,1) + BCs(1,2) + BCs(2,2))/4*ones(n(1)+2,n(2)+2);
% 
% for i = 1:n(2)+2
%     T0(:,i) = (min(min(BCs)) - max(max(BCs)))/(n(2)+2)*i + max(max(BCs));
% end

T0 = zeros(n(1)+2,n(2)+2);


if BCtype(1,1) == 1
    T0(:,1) = BCs(1,1);
end

if BCtype(2,1) == 1
    T0(1,:) = BCs(2,1);
end

if BCtype(1,2) == 1
    T0(:,end) = BCs(1,2);
end

if BCtype(2,2) == 1
    T0(end,:) = BCs(2,2);
end

end

