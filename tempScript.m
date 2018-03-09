%------------------------------------------------------
% Define parameters of problem
l = [0.5 0.5];
start = [0 0];
n = [10 10];

temp = 1;
flux = 2;
BCtype = [temp temp; temp temp];
BCs = [-900 100; -900 500];

% Make k and A a function of space in general
k = 100*ones(n(2)+2, n(1)+2);
A = 1e-3*ones(n(2)+2, n(1)+2);

%------------------------------------------------------
% Define grid spacing for each dimension over the
% whole space
dx = l(1)/n(1)*ones(n(2)+2,n(1)+1);

dx(1,:) = l(1)/(2*n(1));
dx(:,1) = l(1)/(2*n(1));
dx(n(1)+1,:) = l(1)/(2*n(1));
dx(:,n(1)+1) = l(1)/(2*n(1));
%---------------------------------
dy = l(2)/n(2)*ones(n(2)+1,n(1)+2);

dy(2,:) = l(2)/(2*n(2));
dy(:,2) = l(2)/(2*n(2));
dy(n(2)+1,:) = l(2)/(2*n(2));
dy(:,n(2)+1) = l(2)/(2*n(2));

%------------------------------------------------------
% Generate grid points
X = zeros(n(1)+2,n(1)+2);
X(:,1) = start(1);

for i = 2:n(1)+2
    for j = 1:n(2)+2
        X(j,i) = X(j,i-1)+ dx(j,i-1);
    end
end

Y = zeros(n(2)+2,n(2)+2);
Y(end,:) = start(2);

for i = 1:n(1)+2
    for j = 2:n(2)+2
        Y(j,i) = Y(j-1,i)+ dy(j-1,i);
    end
end

%------------------------------------------------------
% Set initial conditions for T
T0 = zeros(n(1)+2,n(2)+2);
T0(:,end) = BCs(1,2);
T0(end,:) = BCs(2,2);

%======================================================
%------------------------------------------------------
% Solve problem
% Inputs: algorithm, initial temp., thermal
% conductivity, area, spacing
% For algorithm: 1 = TDMA, 2 = gaussSeidel
[T, count1] = solve(1, n, T0, BCtype, BCs, k, A, dx, dy);

%------------------------------------------------------
% Set actual solution for comparison


%------------------------------------------------------
% Plot
%plot(x,Tact)
%hold
surf(X,Y,T)
